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
La soluzione proposta considera una soluzione N^2, ogni particella interagisce con tutte le altre particelle e conosce la posizione e la velocità di ognuna di esse. La comunicazione tra le particelle è avvenuta usando la Comunicazione Collettiva, attraverso le funzioni MPI ***Mpi_Bcast*** e ***Mpi_Allgatherv***.  Attraverso l'utiizzo della libreria MPI, è stato possibile parallelizzare il calcolo delle particelle su più processori.  
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

- `int *startRange;`  array che indica qual è l'indice di inizio della porzione di particelle per ogni processore

- `int *countSend;` array che indica quante particelle ogni processore deve calcolare

- `MPI_Datatype myStruct;`  datatype per la struttura Body

- `Body *p = (Body*)buf;`  variabile di tipo Body che punta agli elementi dell'array buf







L'idea di base dell'algoritmo è:

- Il processo MASTER inizializza in maniera casuale la posizione e la velocità di ogni particella. Attraverso la funzione ***Mpi_Bcast***, comunica a tutti i processori i valori delle particelle

