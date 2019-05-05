#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "mt64.h"

#ifndef uniform64
#define uniform64() genrand64_real3()
#endif

typedef struct {
	int height;
	double Lrate, Rrate;
	int Lnbr, Rnbr;
} crystal_site;

int my_itoa(int val, char* buf);
double getRate(crystal_site *h, int i, int whichNbr, double K, int is_dH);
void write2file(char *name, crystal_site *h, int dim);

int main(int argc, char *argv[]) {

	// parse seed from command line
	long long unsigned int sd;
	sd = atoi(argv[1]);
	init_genrand64(sd);
	char seed[20];
	my_itoa(sd, seed);

	int i, j, i2, L, L_loc, whichNbr, numEvents; 
	int siteupdatelist[4];
	const double PI = 3.14159265358979323846264338327950288;
	double R_max, R_loc, c, n, t_stop, K, Tfinal, timetonext, z2, eta;
	int is_dH = 0;
	int nsmpls = 1;
	char parameters[100] = "./parameters.txt";

	FILE *fid;
	fid = fopen(parameters, "r");
	fscanf(fid, "%d %lf %lf %lf",  &L, &K, &Tfinal, &c);
	printf("L=%d, K=%f, T=%e, c = %f\n", L, K, Tfinal, c);
	fclose(fid);
		
	Tfinal *= pow(L, 4);
	n = c*L; // avg number of events we want to occur in t_stop seconds
	L_loc = L;

	crystal_site *h;
	h = (crystal_site *)malloc(L * sizeof(crystal_site));
	h[0].Lnbr = L - 1;
	h[0].Rnbr = 1;
	h[L - 1].Lnbr = L - 2;
	h[L - 1].Rnbr = 0;
	for (i = 1; i<L - 1; i++) {
		h[i].Lnbr = i - 1;
		h[i].Rnbr = i + 1;
	}

	
	//initialize heights
	for (i = 0; i < L; i++) {
		double heightVal = L*sin(2 * PI*((double)i) / ((double)L));
		h[i].height = (int)floor(heightVal);
		if (uniform64() < heightVal - h[i].height)
			h[i].height++;
	}

	char str_hInit[100] = "./hInit.txt";
	write2file(str_hInit, h, L);
	
	double t = 0;
	R_loc = 0;
	for (i = 0; i<L; i++) {
		h[i].Rrate = getRate(h, i, 0, K, is_dH);
		h[i].Lrate = getRate(h, i, 1, K, is_dH);
		R_loc += (h[i].Lrate + h[i].Rrate);
	}

	R_max = R_loc; // this is the case when there's only one block
	t_stop = n / R_max;
	printf("n=%f, t_stop = %f, 1/R_max = %f\n", n, t_stop, 1/R_max);
	clock_t start, end;
	double cpuTime;
	//while(t < 1){
	while (t < Tfinal) {
		//printf("R_loc = %f, current t = %f, t_stop = %f, Tfinal = %f\n", R_loc, t, t_stop, Tfinal);
		if (t + t_stop > Tfinal)
			break;
		double t_1 = 0;
		double t_2 = 0;
		numEvents = 0;
		while (t_1 < t_stop) {
			//numEvents = 0;
			timetonext = -log(uniform64()) / R_loc; //draw from Exp(R_local)
		//	printf("timetonext = %lf\n", timetonext);
			if (t_1 + timetonext > t_stop)
				break;

			//draw from distribution{ j w.p.r_j / R_loc, j = 0,...,L_loc - 1 }
			eta = R_loc*uniform64();
			z2 = 0;
			i = -1;

			while (z2<eta) {
				i++;
				z2 += 2 * h[i].Rrate;
			}
			whichNbr = 0;
			if (uniform64() > 0.5)
				whichNbr = 1;
		//	printf("i = %d, L_loc = %d\n", i, L_loc);
			if (i < L_loc / 2) {
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

				for (j = 0; j< 4; j++) {
					i2 = siteupdatelist[j];
					h[i2].Rrate = getRate(h, i2, 0, K, is_dH);
					h[i2].Lrate = h[i2].Rrate;
				}

				R_loc = 0;
				for (i = 0; i < L; i++)
					R_loc += (h[i].Lrate + h[i].Rrate);
			}

			t_1 += timetonext;

		}
		//printf("num events in left block=%d, expected num = %f\n", numEvents, n/2);
		//COMMUNICATE BETWEEN PROCESSORS
		numEvents = 0;
		while (t_2 < t_stop) {
			timetonext = -log(uniform64()) / R_loc; //draw from Exp(R_local)
			if (t_2 + timetonext > t_stop)
				break;

			//draw from distribution{ j w.p.r_j / R_loc, j = 0,...,L_loc - 1 }
			eta = R_loc*uniform64();
			z2 = 0;
			i = -1;

			while (z2<eta) {
				i++;
				z2 += 2 * h[i].Rrate;
			}
			whichNbr = 0;
			if (uniform64() > 0.5)
				whichNbr = 1;

			if (i >= L_loc / 2) {
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

				for (j = 0; j< 4; j++) {
					i2 = siteupdatelist[j];
					h[i2].Rrate = getRate(h, i2, 0, K, is_dH);
					h[i2].Lrate = h[i2].Rrate;
				}

				R_loc = 0;
				for (i = 0; i < L; i++)
					R_loc += (h[i].Lrate + h[i].Rrate);
			}

			t_2 += timetonext;
		}
		//printf("num events in right block: %d\n", numEvents);
		t += t_stop;
		//COMMUNICATE BETWEEN PROCESSORS 
		//RECOMPUTE R_max AND t_stop
		R_max = R_loc;
		t_stop = n / R_max;
	}

	end = clock();
	cpuTime = (end - start) / (CLOCKS_PER_SEC);
	printf("\n");
	printf("KMC CPU time: %g minutes\n", cpuTime / 60.0);

	char str_hFinal[100] = "./hFinal.txt";
	write2file(str_hFinal, h, L);
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

double getRate(crystal_site *h, int i, int whichNbr, double K, int is_dH) {
	double rate;
	if (is_dH) {
		int h0, hL, hLL, hR, hRR, z1, z2, z3;
		double dH;
		h0 = h[i].height;
		hL = h[h[i].Lnbr].height;
		hLL = h[h[h[i].Lnbr].Lnbr].height;
		hR = h[h[i].Rnbr].height;
		hRR = h[h[h[i].Rnbr].Rnbr].height;

		if (whichNbr == 0) { //jump right
			z1 = h0 - hL;
			z2 = hR - h0;
			z3 = hRR - hR;
			dH = 6 - 2 * z1 + 4 * z2 - 2 * z3;
		}
		else { //jump left
			z1 = hL - hLL;
			z2 = h0 - hL;
			z3 = hR - h0;
			dH = 6 + 2 * z1 - 4 * z2 + 2 * z3;
		}
		rate = exp(-0.5*K*dH);
	}
	else {
		int hL, hR, h0, z1, z2;
		double cord;
		hL = h[h[i].Lnbr].height;
		hR = h[h[i].Rnbr].height;
		h0 = h[i].height;
		z1 = h0 - hL; z2 = hR - h0;
		cord = (double)(z2 - z1 + 1);
		rate = 0.5*exp(-2 * K*cord);
	}
	return rate;
}

void write2file(char *name, crystal_site *h, int dim){
	FILE *fid;
	int i;
	fid = fopen(name, "w");
	for (i = 0; i<dim; i++) {
			fprintf(fid, "%d ", h[i].height);
		//fprintf(fid, "\n");
	}
	fclose(fid);
}
