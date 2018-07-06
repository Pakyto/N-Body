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
  

##### Algoritmo:

L'algoritmo per il calcolo dei valori delle particelle da parte dei processori, è diviso in diverse fasi:

**AVVIO MPI**
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



**ALLOCAZIONE ARRAY**

- Allocazione dell'array dei valori delle particelle con il numero di particelle fornito in input 
- Allocazione degli array countSend e startRange con il numero di processori fornito in input

```c
bytes = nBodies*sizeof(Body);            
buf = (float*)malloc(sizeof(int)*bytes);			
countSend = malloc (sizeof(int)*nproc);				 
startRange = malloc (sizeof(int)*nproc); 
```



**CREAZIONE DATATYPE**
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



**CALCOLO PORZIONE PER OGNI PROCESSORE**
Inizialmente viene prima calcolato il resto e la porzione da assegnare ad ogni processore

```c
resto = nBodies % (nproc);
portion = nBodies/(nproc);
```

Fatto ciò, ogni processore calcola qual è la propria porzione di particelle che deve calcolare e la porzione che deve essere assegnata ad ogni processore, poiché ogni processore attraverso l'utilizzo della ***MPI_Allgatherv*** , deve conoscere qual è la porzione assegnata ad ogni processore, ciò viene effettuato attraverso gli array di puntatori countSend e startRange. Il calcolo della porzione viene effettuato considerando il resto della divisione tra il numero di particelle e il numero di processori:

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



- Il processo MASTER inizializza in maniera casuale la posizione e la velocità di ogni particella. Attraverso la funzione ***MPI_Bcast***, comunica a tutti i processori i valori delle particelle

