# N-body simulation

### Programmazione Concorrente, Parallela e su Cloud

***Università degli Studi di Salerno***

*Laurea Magistrale in Informatica*

*Anno accademico 2017/2018*

**Studente:** Pasquale Settembre

**Matricola:** 0522500554

------

### Problem statement

In un problema n-body, vogliamo trovare la posizione e le velocità di una collezione di particelle che interagiscono tra loro per un periodo di tempo.  Un n-body solver è un programma che trova la soluzione ad un problema n-body attraverso la simulazione del comportamento delle particelle.

### Soluzione proposta

Gli ***input*** all'n-body solver sono, il numero di particelle e il numero di iterazioni che deve simulare, mentre la   posizione e la velocità di ogni particella viene generata in maniera casuale ogni qualvolta viene eseguito l'n-body solver. L'***output*** sarà costituito dalla posizione e dalla velocità di ogni particella alla fine del numero di iterazioni che il programma deve simulare.
La soluzione proposta considera una soluzione N^2, ogni particella interagisce con tutte le altre particelle e conosce la posizione e la velocità di ognuna di esse. La comunicazione tra le particelle è avvenuta usando la Comunicazione Collettiva, attraverso le funzioni MPI ***MPI_Bcast*** e ***MPI_Allgatherv***.  Attraverso l'utiizzo della libreria MPI, è stato possibile parallelizzare il calcolo delle particelle su più processori.  
L'n-body solver sviluppato, è in grado di distribuire in maniera equa ai vari processori una porzione di particelle da calcolare, anche nel caso il numero di processori e il numero di particelle non siano divisibili tra loro.

#### Implementazione

***Variabili utilizzate***:

- Struttura **Body**, corrisponde alla struttura di ogni particella. Dato che stiamo considerando lo spazio 3D, le variabili x,y e z, corrispondono alla ***posizione*** della particella sugli assi x, y e z; mentre le variabli vx, vy e vz, corrispondono alla ***velocità*** della particella sugli assi x, y e z

  ```c
  typedef struct { 				
  	float x, y, z, vx, vy, vz;  							  
  } Body;
  ```

-  `float *buf` array di float con i valori casuali delle particelle

- `int nBodies = atoi(argv[1]);` numero di particelle, passato da linea di comando

- `int nIters = atoi(argv[2]); ` numero iterazioni di simulazione, passato da linea di comando

- `const float dt = 0.01f;` rappresenta l'istante di tempo per ogni iterazione

- `int *startRange;`  array che indica qual è l'indice di inizio della porzione di particelle per ogni processore

- `int *countSend;` array che indica quante particelle ogni processore deve calcolare

- `MPI_Datatype myStruct;`  datatype per la struttura Body

- `Body *p = (Body*)buf;`  variabile di tipo Body  di ogni processore che punta agli elementi dell'array buf
  



##### ALGORITMO:

L'algoritmo per il calcolo dei valori delle particelle da parte dei processori, è diviso in diverse fasi:

**Avvio MPI**

Viene avviato MPI e vengono presi i parametri passati da linea di comando

```c
/* start up MPI */
MPI_Init(&argc, &argv);
/* find out process rank */
MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
/* find out number of processes */
MPI_Comm_size(MPI_COMM_WORLD, &nproc);

int nBodies = atoi(argv[1]);      //numero body
const float dt = 0.01f;  		  //tempo iterazione
const int nIters = atoi(argv[2]); //numero iterazioni
```



**Allocazione array**

- Allocazione dell'array dei valori delle particelle con il numero di particelle fornito in input 
- Allocazione degli array countSend e startRange con il numero di processori fornito in input

```c
bytes = nBodies*sizeof(Body);            
buf = (float*)malloc(sizeof(int)*bytes);			
countSend = malloc (sizeof(int)*nproc);				 
startRange = malloc (sizeof(int)*nproc); 
```



**Creazione datatype**

C'è stata la necessità della creazione di un Datatype poiché tale struttura deve essere trasferita a tutti i processori attraverso MPI, quindi è stato definito il layout del Datatype contenente i campi per la struttura Body. 

```c
const int numitem=6;
int blocklen[6] = {1,1,1,1,1,1};
MPI_Datatype types[6] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};   
MPI_Datatype myStruct;
MPI_Aint     displ[6];
displ[0] = offsetof(Body, x);
displ[1] = offsetof(Body, y);
displ[2] = offsetof(Body, z);
displ[3] = offsetof(Body, vx);
displ[4] = offsetof(Body, vy);
displ[5] = offsetof(Body, vz);
MPI_Type_create_struct(numitem, blocklen, displ, types, &myStruct);
MPI_Type_commit(&myStruct);
```



**Calcolo della porzione delle particelle**

Inizialmente viene prima calcolato il resto e la porzione da assegnare ad ogni processore

```c
resto = nBodies % (nproc);
portion = nBodies/(nproc);
```

Fatto ciò, ogni processore calcola qual è la propria porzione di particelle che deve calcolare, e la porzione che deve essere assegnata ad ogni processore, poiché ogni processore attraverso l'utilizzo della ***MPI_Allgatherv*** , deve conoscere qual è la porzione assegnata di tutti i processori. Il calcolo della porzione viene effettuato considerando il resto della divisione tra il numero di particelle e il numero di processori:

- se il resto > 0, allora l' i-esimo processore avrà come porzione la "porzione assegnata + 1" e decrementa di 1 il valore del resto

In questo modo i primi processori avranno da calcolare una particella in più rispetto alla porzione di particelle assegnata. Successivamente, ogni processore calcola l'indice di inizio della porzione di particelle da calcolare. 

```c
for(int i = 0; i < nproc; i++){
    countSend[i] = portion;
    if(resto > 0){						
        countSend[i] = countSend[i] + 1;
        resto--;
	}
    startRange[i] = count;
    count = count + countSend[i];
}
```



**Inizializzazione particelle**

Il processore MASTER inizializza in maniera casuale la posizione e la velocità di ogni particella. Attraverso la funzione ***MPI_Bcast***, comunica a tutti i processori i valori delle particelle 

```c
if(my_rank == 0){					
    randomizeBodies(buf, 6*nBodies); 
}
MPI_Bcast(buf,bytes,MPI_FLOAT,0,MPI_COMM_WORLD); 
```

**Funzione randomizeBodies**: genera casualmente la posizione e la velocità per ogni particella, memorizzando i valori nell'array di float 

```c
void randomizeBodies(float *data, int n) {  
	for (int i = 0; i < n; i++) {							
		data[i] = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;   
	}
}
```



**Simulazione**

Prima di iniziare la simulazione, ogni processore ha un puntatore di tipo Body che punta all'array dei valori delle particelle conosce i valori di tutte le particelle, in questo modo tutti processori conoscono tutti i valori delle particelle
`Body *p = (Body*)buf;	 `

Dopodiché, la simulazione inizia con il calcolo delle velocità della porzione delle particelle assegnata ad ogni processore. Ogni processore chiama la funzione *bodyForce*, passando le particelle(attraverso il puntatore Body), il tempo per l'iterazione (dt = 0.01), il numero di particelle che deve considerare e l'indice di inizio all'interno dell'insieme delle particelle

```c
void bodyForce(Body *p, float dt, int n, int startRange, int numBodies) {
    for (int i = startRange; i < startRange + n; i++) {         
        float Fx = 0.0f;				   
        float Fy = 0.0f;
        float Fz = 0.0f;
        for (int j = 0; j < numBodies; j++) {      
			float dx = p[j].x - p[i].x;
			float dy = p[j].y - p[i].y;
			float dz = p[j].z - p[i].z;
			float distSqr = dx*dx + dy*dy + dz*dz + SOFTENING;  
			float invDist = 1.0f / sqrtf(distSqr);         
			float invDist3 = invDist * invDist * invDist;   
			Fx += dx * invDist3;   
			Fy += dy * invDist3;
			Fz += dz * invDist3;
		}
		p[i].vx += dt*Fx;
		p[i].vy += dt*Fy;
		p[i].vz += dt*Fz;
	}
}
```

Successivamente, ogni processore calcola la posizione della porzione delle particelle:

```c
for (int i = startRange[my_rank] ; i < startRange[my_rank] + countSend[my_rank]; i++) {
    p[i].x += p[i].vx*dt;
    p[i].y += p[i].vy*dt;
    p[i].z += p[i].vz*dt;
}
```

**Comunicazione tra i processori**: dato che ogni processore ha bisogno della posizione e della velocità di ogni particella, è stata utilizzata la funzione *MPI_Allgatherv* che permette ad ogni processore di ricevere le porzioni di particelle calcolate dagli altri processori, e inviare la porzione di particelle calcolate dall'i-esimo processore a tutti gli altri processori. In questo modo ogni processore avrà i valori di tutte le particelle calcolate dai vari processori. 

```c
	MPI_Allgatherv(MPI_IN_PLACE,0,myStruct,p,countSend,startRange,myStruct,MPI_COMM_WORLD);
```

La funzione Allgatherv prende come parametri:

- *MPI_IN_PLACE* : in questo caso viene utilizzata questa variabile poiché il buffer di input è lo stesso di quello di output; ogni processore inserisce i valori delle particelle calcolate nel buffer di ricezione
- *myStruct*: rappresenta il datatype creato per la struttura Body
- *p*: rappresenta il buffer di ricezione dove verranno memorizzati valori delle particelle di ogni processore
- *countSend:* rappresenta la porzione di particelle di ogni processore, indica quante particelle ci sono nel buffer di ricezione per ogni processore
- *startRange*: rappresenta la posizione, all'interno del buffer di ricezione, delle particelle per ogni processore

Una volta eseguita la Allgatherv, il passo di simulazione termina e i processori rieseguono un nuovo passo di simulazione.

**Stampa dei valori finali delle particelle**

Una volta terminati tutti i passi di simulazione, attraverso l'utilizzo della funzione MPI_Barrier si assicura che tutti i processori terminino la computazione. Il processore MASTER stampa le velocità e le posizioni finali di tutte le particelle e il tempo totale di computazione dei valori delle particelle all'interno della simulazione

```c
if(my_rank == 0){
		for(int i=0; i < nBodies; i++){
			printf("BODY n: %d X:%0.3f Y:%0.3f Z:%0.3f VX:%0.3f VY:%0.3f VZ:%0.3f \n",i,p[i].x,p[i].y,p[i].z,p[i].vx,p[i].vy,p[i].vz);
		}
    printf("\nTEMPO TRASCORSO %f \n",endTime - startTime);
}
/* shut down MPI */
MPI_Type_free(&myStruct);    //Viene liberata la memoria allocata per il datatype creato
//Vengono deallocati i puntatori
free(buf);
free(countSend);
free(startRange);
MPI_Finalize();
return 0;
```

