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
	double rate;

} crystal_site;


int my_itoa(int val, char* buf);
double getRate(crystal_site *h, int i, double K);
void myprint(long N, char* name, crystal_site* h);
double initialize_lattice(crystal_site* h, int L_loc, int rank, int L, double K, int P, MPI_Comm comm);
double KMC_section(int section_num, crystal_site *h, double t_stop, double R_loc, double K, int L_loc, int rank);
int drawSite(crystal_site* h, double R_loc);
void rateMax(crystal_site* h, int rank, int L_loc);

void myprint(long N, char* name, crystal_site* h) {
	FILE *fid;
/*	int* heights = (int*) malloc(N*sizeof(int));
  for(int i = 0; i < N; i++)
    heights[i] = h[i].height;
	int rank, np;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &np);*/
 
  fid = fopen(name, "w+");

  for (long k = 0; k < N; k++) 
		fprintf(fid, "%d ", h[k].height);
  fclose(fid);
	
}

/*
Initialize the lattice for the block of size L_loc = L/P
*/
double initialize_lattice(crystal_site* h, int L_loc, int rank, int L, double K, int P, MPI_Comm comm) {

	MPI_Status status1;

	//initialize heights
	for (int i = 1; i < L_loc + 1; i++) {
		double x = (double)L_loc * rank + i;
		double heightVal = L*sin(2 * PI* x / ((double)L));
		h[i].height = (int)floor(heightVal);
		if (uniform64() < heightVal - h[i].height)
			h[i].height++;
	}


  
	MPI_Send(&(h[1].height), 1, MPI_INT, (rank - 1 + P) % P, 0, comm); // send h[1] to the left
	MPI_Send(&(h[L_loc].height), 1, MPI_INT, (rank + 1) % P, 0, comm); // send h[L_loc] to the right
  
  MPI_Recv(&(h[L_loc + 1].height), 1, MPI_INT, (rank + 1) % P, 0, comm, &status1); //receive to L_loc + 1 from the right
	MPI_Recv(&(h[0].height), 1, MPI_INT, (rank - 1 + P) % P, 0, comm, &status1); //receive to 0 from the left
//	printf("Rank = %d\n", rank);
	// Now send to the right.



	double R_loc = 0;
	for (int i = 1; i < L_loc + 1; i++) {
		h[i].rate = getRate(h, i, K);

		R_loc += 2*h[i].rate;
	}
  
  
 	/*MPI_Send(&(h[1].rate), 1, MPI_DOUBLE, (rank - 1 + P) % P, 0, comm); // send h[1].rate to the left
	MPI_Send(&(h[L_loc].rate), 1, MPI_DOUBLE, (rank + 1) % P, 0, comm); // send h[L_loc].rate to the right
  
  MPI_Recv(&(h[L_loc + 1].rate), 1, MPI_DOUBLE, (rank + 1) % P, 0, comm, &status1); //receive to L_loc + 1 from the right
	MPI_Recv(&(h[0].rate), 1, MPI_DOUBLE, (rank - 1 + P) % P, 0, comm, &status1); //receive to 0 from the left*/

  
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
  double increment = 0;
	double timetonext;
	int i, j, i2, whichNbr, siteupdatelist[4];
	int numEvents = 0;
	while (t < t_stop) {
		//numEvents = 0;
		timetonext = -log(uniform64()) / R_loc; //draw from Exp(R_loc)
    if (t + timetonext > t_stop)
				break;
		
		//draw from distribution{ j w.p.r_j / R_loc, j = 0,...,L_loc - 1 }
    // if (rank==0)
    // printf("about to draw site\n");
  
		i = drawSite(h, R_loc);
	 // if (rank==0)
    // printf("drew site\n");
		whichNbr = 0;
		if (uniform64() > 0.5)
			whichNbr = 1;
		//	printf("Section = %d,  chose i = %d\n", section_num, i);
		if (((i < L_loc/2.0 && section_num==0) || (i >= L_loc/2.0 && section_num == 1))&& i != 0 && i != L_loc + 1) {
     
			numEvents++;
			siteupdatelist[0] = i;
			
      
      if (whichNbr == 1){
				siteupdatelist[1] = i-1;
				siteupdatelist[2] = i-2;
				siteupdatelist[3] = i+1;       
      }
      else{
 				siteupdatelist[1] = i+1;
				siteupdatelist[2] = i-1;
				siteupdatelist[3] = i+2;             
      }

			h[siteupdatelist[0]].height--;
			h[siteupdatelist[1]].height++;
		
			for (j = 0; j< 4; j++) {
				i2 = siteupdatelist[j];
    
				if (i2 > 0 && i2 < L_loc + 1) {
				//	R_loc -= 2*h[i2].rate;
          increment -= 2*h[i2].rate;
					h[i2].rate = getRate(h, i2, K);
				//	R_loc += 2*h[i2].rate;
         increment += 2*h[i2].rate;
				}
			}
		
	
		}

		t += timetonext;

	}

//	printf("Rank %d: num events in section %d = %d\n", rank, section_num, numEvents);
	//return R_loc;
   return increment;

}

void rateMax(crystal_site* h, int rank, int L_loc) {

     int max_idx = 1;
     double max_rate = 0;
     for(int k = 1; k <= L_loc; k++){
        if (h[k].rate >= max_rate){
          max_rate = h[k].rate;
          max_idx = k;
        } 
     }
  printf("Rank %d: new max rate = %f at index %d, heights = %d %d %d\n", rank, max_rate, max_idx, h[max_idx-1].height, h[max_idx].height, h[max_idx+1].height);
  
 
}

int drawSite(crystal_site* h, double R_loc) {
	double eta = R_loc*uniform64();
	double z2 = 0;
	int i = 0;

	while (z2<eta) {
		i++;
		z2 += 2 * h[i].rate;
	}

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

	int L, L_loc, numEvents;
	double R_max, R_loc, c, n, t_stop, K, Tfinal, cpuTime;
        clock_t start, end;
	char parameters[100] = "./parameters.txt";

	FILE *fid;
	fid = fopen(parameters, "r");
	fscanf(fid, "%d %lf %lf %lf", &L, &K, &Tfinal, &c);
	
	fclose(fid);


	Tfinal *= pow(L, 4);
  if (rank == 0)
  printf("L=%d, K=%f, T=%e, c = %f\n", L, K, Tfinal, c);
	n = c * L;

	// Initialize all of the blocks.
	L_loc = L / P; // Get the size of the block.
	crystal_site* h = (crystal_site *)malloc((L_loc + 2) * sizeof(crystal_site));
	R_loc = initialize_lattice(h, L_loc, rank, L, K, P, MPI_COMM_WORLD);


  char str_hInit[100];
    snprintf(str_hInit, 100, "hInit%02d.txt", rank);
  
	myprint(L_loc, str_hInit, h);

	MPI_Allreduce(&R_loc, &R_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);


	t_stop = n / R_max;

  start = clock();
	double t = 0;  // Keep track of the time.
  int num_itr = 0;
	while (t < Tfinal) {
		if (t + t_stop > Tfinal)
		  break;
     // MPI_Barrier(MPI_COMM_WORLD);  
     double increment = 0;
     
	  //if (rank==0)
    //  printf("before section 0 in rank 0, R_loc = %f, R_max = %f\n", R_loc, R_max);
     // printf("starting section 0 in rank with largest R_loc\n");
	 // R_loc = KMC_section(0, h, t_stop, R_loc, K, L_loc, rank);
      increment = KMC_section(0, h, t_stop, R_loc, K, L_loc, rank);
    
	  //if (rank==0)
     // printf("finished section 0 in rank 0, R_loc = %f, R_max = %f\n", R_loc, R_max);
      //printf("finished section 0 in rank with largest R_loc\n");
    
   // R_loc -= 2*h[L_loc].rate;
   // R_loc -= 2*h[L_loc -1].rate;
   increment -=2*h[L_loc].rate;
   increment -= 2*h[L_loc-1].rate;
		//COMMUNICATE BETWEEN PROCESSORS
		MPI_Comm comm = MPI_COMM_WORLD;
		MPI_Status status1;

		MPI_Send(&(h[0].height), 1, MPI_INT, (rank - 1 + P) % P, 123, comm);
		MPI_Send(&(h[1].height), 1, MPI_INT, (rank - 1 + P) % P, 124, comm);
	
		MPI_Recv(&(h[L_loc].height), 1, MPI_INT, (rank + 1) % P, 123, comm, &status1);
		MPI_Recv(&(h[L_loc + 1].height), 1, MPI_INT, (rank + 1) % P, 124, comm, &status1);
   
   
    h[L_loc].rate = getRate(h, L_loc, K);
    h[L_loc-1].rate = getRate(h, L_loc-1, K);
    
  //  R_loc += 2*h[L_loc].rate;
   // R_loc += 2*h[L_loc-1].rate;
    increment +=2*h[L_loc].rate;
   increment += 2*h[L_loc-1].rate;
    //MPI_Barrier(MPI_COMM_WORLD);

	//	R_loc = KMC_section(1, h, t_stop, R_loc, K, L_loc, rank);
    increment = increment + KMC_section(1, h, t_stop, R_loc, K, L_loc, rank);
    R_loc+= increment;
    
   // if (rank==0)
    //  printf("finished section 1 in rank 0, R_loc = %f, R_max = %f\n", R_loc, R_max);
  
    R_loc -= 2*h[1].rate;
    R_loc -= 2*h[2].rate;
  
		t += t_stop;

    //MPI_Barrier(MPI_COMM_WORLD);
		//COMMUNICATE BETWEEN PROCESSORS

		MPI_Send(&(h[L_loc].height), 1, MPI_INT, (rank + 1) % P, 123, comm);
		MPI_Send(&(h[L_loc + 1].height), 1, MPI_INT, (rank + 1) % P, 124, comm);

		MPI_Recv(&(h[0].height), 1, MPI_INT, (rank - 1 + P) % P, 123, comm, &status1);
		MPI_Recv(&(h[1].height), 1, MPI_INT, (rank - 1 + P) % P, 124, comm, &status1);
    
    
    h[1].rate = getRate(h, 1, K);
    h[2].rate = getRate(h,2,K);

    R_loc += 2*h[1].rate;
    R_loc += 2*h[2].rate;

		//RECOMPUTE R_max AND t_stop
		MPI_Allreduce(&R_loc, &R_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
		if (rank == 0)
	  	printf("num_itr = %d, R_max = %f, new t = %f\n", num_itr, R_max, t);
   // if (num_itr == 0)
   // rateMax(h, rank, L_loc);                 
		t_stop = n / R_max;

    num_itr++;
	}

	end = clock();
	cpuTime = (end - start) / (CLOCKS_PER_SEC);
  if (rank==0){
	printf("\n");
	printf("KMC CPU time: %g minutes\n", cpuTime / 60.0);
  }


  char str_hFinal[100];
  snprintf(str_hFinal, 100, "hFinal%02d.txt", rank);
	myprint(L_loc, str_hFinal, h);
 
  MPI_Finalize();

}
