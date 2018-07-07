/*
 ============================================================================
 Name        : N-Body.c
 Author      : Pasquale Settembre
 Version     :
 Copyright   : Your copyright notice
 Description : N-body solver in C
 ============================================================================
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#define SOFTENING 1e-9f

typedef struct { 				//STRUTTURA BODY
	float x, y, z, vx, vy, vz;  //x, y, z sono le posizioni nello spazio 3D;
								//vx, vy, vz sono le velocità sulle cordinate x,y e z
} Body;

//Calcolo dei valori dei bodies in maniera casuale
void randomizeBodies(float *data, int n) {  //La funzione prende in input un array di float e il numero di bodies
	for (int i = 0; i < n; i++) {							//Per n body, si calcolano i valori dei bodies
		data[i] = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;   //Allocazione degli elementi dell'array con valori casuali
	}
}

/**Funzione bodyForce, prende come paramentri:
 * i bodies da considerare
 * il tempo dell'iterazione
 * la porzione del numero di body assegnata al processore
 * l'indice di inizio di un processore all'interno della struttura body
 * numero totale i bodies
 */
void bodyForce(Body *p, float dt, int n, int startRange, int numBodies) {

	for (int i = startRange; i < startRange + n; i++) {          //Calcolo delle forze per la porzione di bodies di un processore

		float Fx = 0.0f;				   //Si inizializza a 0 la forza su x,y e z
		float Fy = 0.0f;
		float Fz = 0.0f;

		for (int j = 0; j < numBodies; j++) {      //Calcolo della distanza da un determinato body verso tutti gli altri


			/*CALCOLO DISTANZA su x,y e z (si calcola prendendo la posizione x, y, z di un altro body
			 * 									-
												la posizione x, y, z del body corrente)*/
			float dx = p[j].x - p[i].x;
			float dy = p[j].y - p[i].y;
			float dz = p[j].z - p[i].z;

			float distSqr = dx*dx + dy*dy + dz*dz + SOFTENING;  //Calcolo distanza totale da un body
			float invDist = 1.0f / sqrtf(distSqr);          //Calcolo dell'inversa della distanza
			float invDist3 = invDist * invDist * invDist;   //L'inversa si moltiplica per 3, perché è uno spazio a 3 dimensioni

			Fx += dx * invDist3;   //Calcola la forza su x come, distanza di x * l'inverso della distanza
			Fy += dy * invDist3;
			Fz += dz * invDist3;
		}
		//CALCOLO VELOCITÀ SU x,y,z
		p[i].vx += dt*Fx;
		p[i].vy += dt*Fy;
		p[i].vz += dt*Fz;

	}

}

int main(int argc, char* argv[]){
	int  my_rank; /* rank of process */
	int  nproc;       /* number of processes */

	int resto;
	int *startRange;
	int *countSend;
	int portion;
	double startTime,endTime;
	int bytes;
	float *buf;

	/* start up MPI */
	MPI_Init(&argc, &argv);

	/* find out process rank */
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	/* find out number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	int nBodies = atoi(argv[1]);      //numero body
	const float dt = 0.01f;  		  //tempo iterazione
	const int nIters = atoi(argv[2]);  //numero iterazioni
	int count = 0;


	bytes = nBodies*sizeof(Body);             //Allocazione numero di byte in base alla dimensione di Body
	buf = (float*)malloc(sizeof(int)*bytes);			  //Allocazione array buffer di puntatori con il numero di byte dei body


	countSend = malloc (sizeof(int)*nproc);				//Array che indica quanti elementi ha ogni processore
	startRange = malloc (sizeof(int)*nproc);		    //Array che indica quale l'indice di inizio di un processore

	/*---DATATYPE CREATION---*/
	/*Viene effettuata la creazione di un Datatype per la struttura Body, perché la struttura deve essere serializzata
	e dovrà essere inviata agli altri processori attraverso MPI*/

	const int numitem=6;
	int blocklen[6] = {1,1,1,1,1,1};
	MPI_Datatype types[6] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};   //Creazione datatype con i campi della struttura Body
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
	/*---END DATATYPE CREATION---*/

	Body *p = (Body*)buf;					      //Puntatore elemento Body che contiene l'array buf

	resto = nBodies % (nproc);
	portion = nBodies/(nproc);

	/**Calcolo del resto da parte di tutti i processori:
	 * ogni processore conosce qual è la porzione degli altri processori
	 */

	for(int i = 0; i < nproc; i++){
		countSend[i] = portion;
		//printf("RESTO FOR %d \n",resto);
		if(resto > 0){						 //Se il resto è maggiore di 0, i primi processori avranno un Body in più
			countSend[i] = countSend[i] + 1;
			resto--;
		}

		startRange[i] = count;
		count = count + countSend[i];

	}


	if(my_rank == 0){					 //Il master effettua l'inizializzazione della struttura Body
		randomizeBodies(buf, 6*nBodies); //Inizializzazione attraverso l'array di float buf
	}

	MPI_Bcast(buf,bytes,MPI_FLOAT,0,MPI_COMM_WORLD);    //Viene inviato a tutti i processori l'array inizializzato

	startTime = MPI_Wtime();

	for (int iter = 1; iter <= nIters; iter++) {              //INIZIO SIMULAZIONE

		printf("ITER %d PROCESSORE %d  Indice di inizio %d Num.particelle da calcolare %d \n",iter,my_rank,startRange[my_rank],countSend[my_rank]);

		//Calcolo della forza dei Body da parte di un processore
		bodyForce(p,dt,countSend[my_rank],startRange[my_rank],nBodies);


		//CALCOLO DELLE POSIZIONI x,y,z della porzione di bodies di un processore
		for (int i = startRange[my_rank] ; i < startRange[my_rank] + countSend[my_rank]; i++) {
			p[i].x += p[i].vx*dt;
			p[i].y += p[i].vy*dt;
			p[i].z += p[i].vz*dt;
		}

		//Viene inviata a tutti i processori la porzione di bodies modificati da ogni processore, ogni processore conosce le posizioni e la forze dei bodies calcolate dagli altri proce
		/**
		 * Vengono utilizzati i puntatori countSend e startRange, perché ogni processore conosce qual è la porzione assegnata agli altri processori
		 * Utilizzando la Allgatherv, viene inviata ad ogni procesore la porzione calcolata da un processore
		 * In questo modo ogni processore conosce le posizioni e le forze dei bodies calcolati dagli altri processori
		 */
		MPI_Allgatherv(MPI_IN_PLACE,0,myStruct,p,countSend,startRange,myStruct,MPI_COMM_WORLD);

	}



	MPI_Barrier(MPI_COMM_WORLD); //Si attende che tutti processori abbiano completato la simulazione

	endTime = MPI_Wtime();

	//Il master stampa il valore finale dei bodies e il tempo della simulazione
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
}

