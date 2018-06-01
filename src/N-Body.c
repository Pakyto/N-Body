/*
 ============================================================================
 Name        : N-Body.c
 Author      : Pasquale Settembre
 Version     :
 Copyright   : Your copyright notice
 Description : Hello MPI World in C
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
		//data[i] = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;   //Allocazione degli elementi dell'array con valori casuali
		data[i] = 1;
	}
	printf("USCITA randomizeBodies \n");

}

//Calcolo della forza di un body
void bodyForce(Body *p, float dt, int n, int startRange, int numBodies) {       //Prende in input un body, il tempo trascorso e il numero di body
	printf("FUNZIONE bodyForce COUNTSEND: %d  STARTRANGE: %d \n",n,startRange);

	for (int i = startRange; i < startRange + n; i++) {          //Per ogni body

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
		printf("USCITA bodyForce %0.3f %0.3f %0.3f \n",p[i].vx,p[i].vy,p[i].vz);

	}
	//QUESTA FASE DEVE ESSERE PARALLELIZZATA CON MPI
}

int main(int argc, char* argv[]){
	int  my_rank; /* rank of process */
	int  nproc;       /* number of processes */
	MPI_Status status ;   /* return status for receive */

	int resto;
	int *startRange;
	int *countSend;
	int portion;
	double startTime,endTime;

	/* start up MPI */
	MPI_Init(&argc, &argv);

	/* find out process rank */
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	/* find out number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	int nBodies = 50;      //numero body
	const float dt = 0.01f; // time step (conteggia il tempo trascorso)
	const int nIters = 5;  // simulation iterations
	int count = 0;

	if (argc > 1)
		nBodies = atoi(argv[1]);

	int bytes = nBodies*sizeof(Body);             //Allocazione numero di byte in base alla dimensione di Body
	float *buf = (float*)malloc(bytes);			  //Allocazione array buffer di puntatori con il numero di byte dei body

	countSend = malloc (sizeof(int)*nproc);
	startRange = malloc (sizeof(int)*nproc);

	/*---DATATYPE CREATION---*/
	const int numitem=6;
	int blocklen[6] = {1,1,1,1,1,1};
	MPI_Datatype types[6] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
	MPI_Datatype body_type;
	MPI_Aint     displ[6];

	displ[0] = offsetof(Body, x);
	displ[1] = offsetof(Body, y);
	displ[2] = offsetof(Body, z);
	displ[3] = offsetof(Body, vx);
	displ[4] = offsetof(Body, vy);
	displ[5] = offsetof(Body, vz);

	MPI_Type_create_struct(numitem, blocklen, displ, types, &body_type);
	MPI_Type_commit(&body_type);
	/*-------------------*/

	Body *p = (Body*)buf;					      //Puntatore ad un elemento Body che contiene l'array buf

	resto = nBodies % (nproc);
	printf("RESTO %d \n",resto);
	portion = nBodies/(nproc);
	printf("PORTION %d \n\n",portion);

	for(int i = 0; i < nproc; i++){
		countSend[i] = portion;
		printf("RESTO FOR %d \n",resto);
		if(resto > 0){
			countSend[i] = countSend[i] + 1;   //Alcuni processori hanno un elemento in più
			resto--;
		}
		printf("COUNTSEND pos %d %d \n",i,countSend[i]);

		startRange[i] = count;
		printf("startRange %d \n",startRange[i]);

		count = count + countSend[i];
		printf("COUNT %d \n\n",count);

	}

	MPI_Barrier(MPI_COMM_WORLD);
	printf("DOPO BARRIERA \n");

	float recvbuf[countSend[my_rank]];

	if(my_rank == 0){
		randomizeBodies(buf, 6*nBodies); // Init pos / vel data  (calcolo dei valori dei bodies)
	}

	printf("PRIMA SCATTER \n");

	MPI_Scatterv(buf,countSend,startRange,MPI_FLOAT,recvbuf,countSend[my_rank],MPI_FLOAT,0,MPI_COMM_WORLD);

	//MPI_Scatterv(&buf,countSend,startRange,MPI_FLOAT,&recvbuf,countSend[my_rank],MPI_FLOAT,0,MPI_COMM_WORLD);
	printf("dopo SCATTER \n");


	for (int i = 0; i < countSend[my_rank]; i++) {
		printf("RANK %d BUFFF %f \n",my_rank, recvbuf[i]);
	}

	printf("\n");


	startTime = MPI_Wtime();

	for (int iter = 1; iter <= nIters; iter++) {              //INIZIO SIMULAZIONE

		printf("PROCESSOR %d \n",my_rank);
		bodyForce(p,dt,countSend[my_rank],startRange[my_rank],nBodies);


		/*
		//CALCOLO DELLE POSIZIONI x,y,z dei bodies
		for (int i = 0 ; i < nBodies; i++) { // integrate position
			//Si calcola, velocità su x per dt(tempo trascorso) e il risultato viene sommato alla posizione x
			p[i].x += p[i].vx*dt;
			p[i].y += p[i].vy*dt;
			p[i].z += p[i].vz*dt;
		}
		*/
	}
	endTime = MPI_Wtime();


	MPI_Barrier(MPI_COMM_WORLD);
	printf("DOPO BARRIERA \n");

	if(my_rank == 0){
		printf("\nTEMPO TRASCORSO %f \n",endTime - startTime);
		/*for(int i=0; i < nBodies; i++){
			printf("BODY n: %d %0.3f %0.3f %0.3f %0.3f \n",i,p[i].x,p[i].y,p[i].z,p[i].vx);
		}*/

	}



	/* shut down MPI */
	free(buf);
	free(countSend);
	free(startRange);
	MPI_Finalize();

	return 0;
}

