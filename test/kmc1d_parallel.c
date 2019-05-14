/*******************************************************************************
Run 1-D KMC in parallel with p processors.
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <string.h>
#include "mt64.h"

#ifndef uniform64
#define uniform64() genrand64_real3()
#endif

const double PI = 3.14159265358979323846264338327950288;

typedef struct {
	int height;
	double Lrate, Rrate;
	int Lnbr, Rnbr;
} crystal_site;


int my_itoa(int val, char* buf);

double getRate(crystal_site *h, int i, double K);

void myprint(long N, char* name, crystal_site* h, MPI_Comm comm);

double initialize_lattice(crystal_site* h, int L_loc, int rank, int L, double K, int P, MPI_Comm comm);
double KMC_section(int section_num, crystal_site *h, double t_stop, double R_loc, double K, int L_loc, int rank);
int drawSite(crystal_site* h, double R_loc);






void myprint(long N, char* name, crystal_site* h, MPI_Comm comm) {
	FILE *fid;
	fid = fopen(name, "w");
  int rank, np;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &np);

  for (int i = 0; i < np; i++) {
    MPI_Barrier(comm);
    if (rank == i) {
      for (long k = 0; k < N; k++) {
        fprintf(fid, "%d ", h[k].height);
      }
    }
    MPI_Barrier(comm);
  }

	fclose(fid);
}

/*
Initialize the lattice for the block of size L_loc = L/P
*/
double initialize_lattice(crystal_site* h, int L_loc, int rank, int L, double K, int P, MPI_Comm comm) {

	MPI_Status status1;

	//initialize heights
	for (int i = 1; i < L_loc + 1; i++) {
		double x = (double) L_loc * rank + i;
		double heightVal = L*sin(2 * PI* x / ((double)L));
		h[i].height = (int)floor(heightVal);
		if (uniform64() < heightVal - h[i].height)
			h[i].height++;
	}


	// Send to the left first.
	MPI_Send(&(h[1].height), 1, MPI_INT, (rank - 1 + P) % P, 0, comm);
	MPI_Recv(&(h[L_loc + 1].height), 1, MPI_INT, (rank + 1) % P, 0, comm, &status1);

	MPI_Send(&(h[L_loc].height), 1, MPI_INT, (rank + 1) % P, 0, comm);
	MPI_Recv(&(h[0].height), 1, MPI_INT, (rank - 1 + P) % P, 0, comm, &status1);
printf("Rank = %d\n", rank);
	// Now send to the right.



	double R_loc = 0;
	for (int i = 1; i < L_loc + 1; i++) {
		h[i].Rrate = getRate(h, i, K);
		h[i].Lrate = h[i].Rrate;
		R_loc += (h[i].Lrate + h[i].Rrate);
	}

	return R_loc;
}


double getRate(crystal_site *h, int i, double K) {
	double rate;

	int hL, hR, h0, z1, z2;
	double cord;
	hL = h[i - 1].height;
	hR = h[i + 1].height;
	h0 = h[i].height;
	z1 = h0 - hL; z2 = hR - h0;
	cord = (double)(z2 - z1 + 1);
	rate = 0.5*exp(-2 * K*cord);

	return rate;
}

int my_itoa(int val, char* buf)
{
	const unsigned int radix = 10;

	char* p;
	unsigned int a;        //every digit
	int len;
	char* b;            //start of the digit char
	char temp;
	unsigned int u;

	p = buf;

	if (val < 0)
	{
		*p++ = '-';
		val = 0 - val;
	}
	u = (unsigned int)val;

	b = p;

	do
	{
		a = u % radix;
		u /= radix;

		*p++ = a + '0';

	} while (u > 0);

	len = (int)(p - buf);

	*p-- = 0;

	//swap
	do
	{
		temp = *p;
		*p = *b;
		*b = temp;
		--p;
		++b;

	} while (b < p);

	return len;
}

double KMC_section(int section_num, crystal_site *h, double t_stop, double R_loc, double K, int L_loc, int rank) {
	double t = 0;

	double timetonext;
	int i, j, i2, whichNbr, siteupdatelist[4];
	int numEvents = 0;
	while (t < t_stop) {
		if (rank ==0)
			printf("section = %d, t = %f, t_stop = %f\n", section_num, t, t_stop);
		//numEvents = 0;
		timetonext = -log(uniform64()) / R_loc; //draw from Exp(R_local)

		//if (t + timetonext > t_stop)
			//break;
		if (rank ==0)
			printf("drew a time\n");
		//draw from distribution{ j w.p.r_j / R_loc, j = 0,...,L_loc - 1 }

		i = drawSite(h, R_loc);
		if(rank==0)
   	printf("section = %d, selected index i = %d\n", section_num, i);
		whichNbr = 0;
		if (uniform64() > 0.5)
			whichNbr = 1;
		//	printf("Section = %d,  chose i = %d\n", section_num, i);
		if ((int)(2*i/L_loc) == section_num && i != 0 && i != L_loc + 1) {
			printf("index i = %d was chosen\n", i);
			//printf("i=%d, 2*i/L_loc = %d, section_num = %d\n", i, 2*i/L_loc, section_num);
			numEvents++;
			siteupdatelist[0] = i;
			if (whichNbr == 1) {
				siteupdatelist[1] = h[siteupdatelist[0]].Lnbr;
				siteupdatelist[2] = h[siteupdatelist[1]].Lnbr;
				siteupdatelist[3] = h[siteupdatelist[0]].Rnbr;
			}
			else {
				siteupdatelist[1] = h[siteupdatelist[0]].Rnbr;
				siteupdatelist[2] = h[siteupdatelist[0]].Lnbr;
				siteupdatelist[3] = h[siteupdatelist[1]].Rnbr;
			}

			h[siteupdatelist[0]].height--;
			h[siteupdatelist[1]].height++;
			if(rank==0)
      printf("got here\n");
			for (j = 0; j< 4; j++) {
				i2 = siteupdatelist[j];
				if (i2 != 0 && i2 != L_loc + 1) {
					R_loc -= h[i2].Rrate;
					R_loc -= h[i2].Lrate;
					h[i2].Rrate = getRate(h, i2, K);
					h[i2].Lrate = h[i2].Rrate;
					R_loc += h[i2].Rrate;
					R_loc += h[i2].Lrate;
				}
			}
			if(rank==0)
			printf("got here\n");

			/*R_loc = 0;
			for (i = 0; i < L; i++)
				R_loc += (h[i].Lrate + h[i].Rrate);*/
		}

		t += timetonext;

	}
	//printf("num events in section %d = %d\n", section_num, numEvents);
	return R_loc;

}


int drawSite(crystal_site* h, double R_loc) {
	double eta = R_loc*uniform64();
	double z2 = 0;
	int i = 0;

	while (z2<eta) {
		i++;
		z2 += 2 * h[i].Rrate;
	}
	int whichNbr = 0;
	if (uniform64() > 0.5)
		whichNbr = 1;
	return i;
}

int main(int argc, char * argv[]) {

	int rank; // Which processor is executing currently.
	int P;		// The number of blocks/processors.

	MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &P);

	// parse seed from command line
	long long unsigned int sd;
	sd = atoi(argv[1]);
	init_genrand64(sd);
	char seed[20];
	my_itoa(sd, seed);

	int i, j, i2, L, L_loc, numEvents;
	int siteupdatelist[4];
	double R_max, R_loc, c, n, t_stop, K, Tfinal, timetonext, z2, eta, cpuTime, start, end;
	char parameters[100] = "./parameters.txt";

	FILE *fid;
	fid = fopen(parameters, "r");
	fscanf(fid, "%d %lf %lf %lf",  &L, &K, &Tfinal, &c);
	printf("L=%d, K=%f, T=%e, c = %f\n", L, K, Tfinal, c);
	fclose(fid);


	Tfinal *= pow(L, 4);
	n = c * L;

	// Initialize all of the blocks.
	L_loc = L / P; // Get the size of the block.
	crystal_site* h = (crystal_site *) malloc((L_loc + 2) * sizeof(crystal_site));
	R_loc = initialize_lattice(h, L_loc, rank, L, K, P, MPI_COMM_WORLD);

	char str_hInit[100] = "./hInit.txt";
	myprint(L_loc, str_hInit, h, MPI_COMM_WORLD);

	MPI_Allreduce(&R_loc, &R_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);


	t_stop = n / R_max;


	double t = 0;  // Keep track of the time.
	while (t < Tfinal) {
		if (t + t_stop > Tfinal)
			break;
			if (rank==0)
		//printf("Starting section 0: rank = %d, t = %f, Tfinal = %f\n", rank, t, Tfinal);
		R_loc = KMC_section(0, h, t_stop, R_loc, K, L_loc, rank);
		//COMMUNICATE BETWEEN PROCESSORS
		MPI_Comm comm = MPI_COMM_WORLD;
		MPI_Status status1;
			if (rank==0)
    //printf("After section 0: rank = %d, t = %f, Tfinal = %f\n", rank, t, Tfinal);
		MPI_Send(&(h[0].height), 1, MPI_INT, (rank - 1 + P) % P, 123, comm);
		MPI_Send(&(h[1].height), 1, MPI_INT, (rank - 1 + P) % P, 124, comm);
			if (rank==0)
   //printf("After sending: rank = %d, t = %f, Tfinal = %f\n", rank, t, Tfinal);
		MPI_Recv(&(h[L_loc].height), 1, MPI_INT, (rank + 1) % P, 123, comm, &status1);
		MPI_Recv(&(h[L_loc + 1].height), 1, MPI_INT, (rank + 1) % P, 124, comm, &status1);
			if (rank==0)
  // printf("After receiving: rank = %d, t = %f, Tfinal = %f\n", rank, t, Tfinal);
		R_loc = KMC_section(1, h, t_stop, R_loc, K, L_loc, rank);
		t += t_stop;


		//COMMUNICATE BETWEEN PROCESSORS

		MPI_Send(&(h[L_loc].height), 1, MPI_INT, (rank + 1) % P, 123, comm);
		MPI_Send(&(h[L_loc + 1].height), 1, MPI_INT, (rank + 1) % P, 124, comm);

		MPI_Recv(&(h[0].height), 1, MPI_INT, (rank - 1 + P) % P, 123, comm, &status1);
		MPI_Recv(&(h[1].height), 1, MPI_INT, (rank - 1 + P) % P, 124, comm, &status1);

		//RECOMPUTE R_max AND t_stop
		MPI_Allreduce(&R_loc, &R_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		if (rank ==0)
		printf("rank = %d, t = %f, Tfinal = %f\n", rank, t, Tfinal);
		t_stop = n / R_max;
	}

	end = clock();
	cpuTime = (end - start) / (CLOCKS_PER_SEC);
	printf("\n");
	printf("KMC CPU time: %g minutes\n", cpuTime / 60.0);

	char str_hFinal[100] = "./hFinal.txt";
	myprint(L_loc, str_hFinal, h, MPI_COMM_WORLD);

}
