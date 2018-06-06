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
void bodyForce(Body *p, float dt, int n, int startRange, int numBodies,MPI_Datatype body_type,int rank) {       //Prende in input un body, il tempo trascorso e il numero di body
	printf("FUNZIONE bodyForce COUNTSEND: %d  STARTRANGE: %d \n",n,startRange);

	printf("STRUTTURA BODY bodyForce \n");
	for (int i = startRange; i < startRange+n; i++) {
		printf("RANK %d bodyForce BODIES %0.3f %0.3f %0.3f  \n",rank, p[i].x, p[i].y, p[i].z);
	}

	MPI_Barrier(MPI_COMM_WORLD);


	printf("J BODY \n");
	for (int i = 0; i < numBodies; i++) {
		printf("RANK %d j BODIES %0.3f %0.3f %0.3f  \n",rank, p[i].x, p[i].y, p[i].z);
	}

	MPI_Barrier(MPI_COMM_WORLD);


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

		printf("RANK %d USCITA bodyForce %0.3f %0.3f %0.3f \n",rank,p[i].vx,p[i].vy,p[i].vz);
	}
	//MPI_Allgather(p,sizeof(*p),body_type,p,sizeof(*p),body_type,MPI_COMM_WORLD);


	//QUESTA FASE DEVE ESSERE PARALLELIZZATA CON MPI
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

	int nBodies = 30;      //numero body
	const float dt = 0.01f; // time step (conteggia il tempo trascorso)
	const int nIters = 5;  // simulation iterations
	int count = 0;

	if (argc > 1)
		nBodies = atoi(argv[1]);


	bytes = nBodies*sizeof(Body);             //Allocazione numero di byte in base alla dimensione di Body
	buf = (float*)malloc(sizeof(int)*bytes);			  //Allocazione array buffer di puntatori con il numero di byte dei body


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

	MPI_Datatype myStruct;
	MPI_Type_create_resized( body_type, 0, sizeof( Body ), &myStruct );

	MPI_Type_commit(&myStruct);
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

	//printf("DOPO BARRIERA \n");

	//float *recvbuf = (float*)malloc(sizeof(int)*bytes);

	//float recvbuf[countSend[my_rank]];                   //Buffer con il numero di elementi assegnato ad ogni proccesso

	if(my_rank == 0){
		randomizeBodies(buf, 6*nBodies); // Init pos / vel data  (calcolo dei valori dei bodies)
	}

	printf("BARRIER PRIMA BCAST \n");
	MPI_Bcast(buf,bytes,MPI_FLOAT,0,MPI_COMM_WORLD);

	printf("PRIMA SCATTER \n");


	//MPI_Scatterv(buf,countSend,startRange,MPI_FLOAT,recvbuf,countSend[my_rank],MPI_FLOAT,0,MPI_COMM_WORLD);


	printf("dopo SCATTER \n");


	/*
	for (int i = startRange[my_rank]; i < startRange[my_rank]+countSend[my_rank]; i++) {
		printf("RANK %d BUFFF %0.3f \n",my_rank, buf[i]);
	}*/

	/*
	for (int i = startRange[my_rank]; i < startRange[my_rank]+countSend[my_rank]; i++) {
		printf("RANK %d BODIES %0.3f %0.3f %0.3f  \n",my_rank, p[i].x, p[i].y, p[i].z);
	}*/

	printf("\n");

	MPI_Barrier(MPI_COMM_WORLD);


	startTime = MPI_Wtime();

	Body *bodies = (Body*)buf;

	for (int i = startRange[my_rank]; i < startRange[my_rank]+countSend[my_rank]; i++) {
		printf("RANK %d BODIES %0.3f %0.3f %0.3f  \n",my_rank, bodies[i].x, bodies[i].y, bodies[i].z);
	}

	for (int iter = 1; iter <= nIters; iter++) {              //INIZIO SIMULAZIONE

		printf("ITER %d PROCESSOR %d \n",iter,my_rank);

		bodyForce(p,dt,countSend[my_rank],startRange[my_rank],nBodies,body_type,my_rank);


		//CALCOLO DELLE POSIZIONI x,y,z dei bodies
		for (int i = startRange[my_rank] ; i < startRange[my_rank] + countSend[my_rank]; i++) { // integrate position
			//Si calcola, velocità su x per dt(tempo trascorso) e il risultato viene sommato alla posizione x
			p[i].x += p[i].vx*dt;
			p[i].y += p[i].vy*dt;
			p[i].z += p[i].vz*dt;
		}

		printf("PRIMA ALLGATHER \n");
		for (int i = 0; i < nBodies; i++) {
			printf("RANK %d BODIES n %d %0.3f %0.3f %0.3f  \n",my_rank,i, p[i].x, p[i].y, p[i].z);
		}


		MPI_Allgatherv(MPI_IN_PLACE,0,myStruct,p,countSend,startRange,myStruct,MPI_COMM_WORLD);

		printf("DOPO ALLGATHER \n");
		for (int i = 0; i < nBodies; i++) {
			printf("RANK %d BODIES n.%d %0.3f %0.3f %0.3f  \n",my_rank,i, p[i].x, p[i].y, p[i].z);
		}

		MPI_Barrier(MPI_COMM_WORLD);

	}



	MPI_Barrier(MPI_COMM_WORLD);
	//printf("DOPO BARRIERA \n");

	endTime = MPI_Wtime();



	if(my_rank == 0){
		printf("\nTEMPO TRASCORSO %f \n",endTime - startTime);
		for(int i=0; i < nBodies; i++){
			printf("BODY n: %d X:%0.3f Y:%0.3f Z:%0.3f VX:%0.3f VY:%0.3f VZ:%0.3f \n",i,p[i].x,p[i].y,p[i].z,p[i].vx,p[i].vy,p[i].vz);
		}

	}


	/* shut down MPI */
	MPI_Type_free(&body_type);

	free(buf);
	free(countSend);
	free(startRange);
	MPI_Finalize();

	return 0;
}

