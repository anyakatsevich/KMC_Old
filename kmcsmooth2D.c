#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "mt64.h"

#ifndef uniform64
#define uniform64() genrand64_real3()
#endif



typedef struct{
  int height;
  double rate;
  double cord;
  int Lnbr, Rnbr, Unbr, Dnbr;
} crystal_site;




double potential(int z)
{
//	double p = 1.0;
//	return pow(fabs((double)z), p);
	
	return (double)(z*z);
}



double getcord(crystal_site *h, int i)
{
  int h0, hL,hR,hU,hD;
  double cord;
	
  h0 = h[i].height;
  hL = h[h[i].Lnbr].height;
  hR = h[h[i].Rnbr].height;
  hD = h[h[i].Dnbr].height;
  hU = h[h[i].Unbr].height;
  
  /* cord = 0; */

  /* cord = 0.5*(  potential(hR-h0+1, p)-potential(hR-h0, p) */
  /* 		+ potential(h0-hL-1, p)-potential(h0-hL, p) ); */

  cord = 0.5*(  potential(hR-h0+1)-potential(hR-h0)
  		+ potential(h0-hL-1)-potential(h0-hL)
  		+ potential(hU-h0+1)-potential(hU-h0)
  		+ potential(h0-hD-1)-potential(h0-hD) );

  /* if (cord>=15.0){ */
  /*   printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",  */
  /* 	   i,h[i].Lnbr,h[i].Rnbr,h[i].Dnbr,h[i].Unbr,h0,hL,hR,hD,hU); */
  /*   exit(1); */
  /* } */
	
  return cord;
}





int main (void)
{
  int L = 100;
  long nsmpls = 1;
	
  double K = 1.5;
  double T = (1e-4)*(double)pow(L,4);
	
  //	printf("%g\n",T);
  //	exit(1);
	
  long m;
  int i, j, i1, j1, j2, i2, j3, k, nbr0, nbr1, nbr2, oldlistpos, updates;
  double z1,z2, eta, timetonext, t, ratesum, oldrate, injecttime, r1, oldcord;
  int siteupdatelist[8];
	
  const double PI = 3.14159265358979323846264338327950288;
	
  init_genrand64(4);
  //srandom(5);
	
  int **map, *map_mem;
  double *hini, *rates, *avh, *cords;
	
  crystal_site *h;

  map = (int **)malloc(L*sizeof(int*));
  map_mem = (int *)malloc(L*L*sizeof(int));
  for(i=0;i<L;i++)
    map[i] = &map_mem[i*L];

  hini = (double *)malloc(L*L*sizeof(double));
  rates = (double *)malloc(L*L*sizeof(double));
  avh = (double *)malloc(L*L*sizeof(double));
  cords = (double *)malloc(L*L*sizeof(double));

  h = (crystal_site *)malloc(L*L*sizeof(crystal_site));

  /* int map[L][L]; */
  /* double hini[L*L], rates[L*L], avh[L*L], cords[L*L]; */

  /* crystal_site h[L*L]; */
	
  FILE *fid;
	
  for(i=0;i<L;i++){
    for(j=0;j<L;j++){
      map[i][j] = i*L+j;
    }
  }
	
  for(i=0;i<L;i++){
    for(j=0;j<L;j++){
      //   hini[i] = 0;
      hini[map[i][j]] = 0*L*sin(2*PI*i/L)*sin(2*PI*j/L);
      //hini[map[i][j]] = 0;
      avh[map[i][j]] = 0;
      rates[map[i][j]] = 0;
      cords[map[i][j]] = 0;
      //    rates0[i] = 0;
    }
  }
	
  /////////////////////////////////////////////
  //////////// The neighbor maps ///////////////

  /////////// corners
  h[map[0][0]].Lnbr = map[L-1][0];
  h[map[0][0]].Rnbr = map[1][0];
  h[map[0][0]].Dnbr = map[0][L-1];
  h[map[0][0]].Unbr = map[0][1];

  h[map[0][L-1]].Lnbr = map[L-1][L-1];
  h[map[0][L-1]].Rnbr = map[1][L-1];
  h[map[0][L-1]].Dnbr = map[0][L-2];
  h[map[0][L-1]].Unbr = map[0][0];

  h[map[L-1][0]].Lnbr = map[L-2][0];
  h[map[L-1][0]].Rnbr = map[0][0];
  h[map[L-1][0]].Dnbr = map[L-1][L-1];
  h[map[L-1][0]].Unbr = map[L-1][1];

  h[map[L-1][L-1]].Lnbr = map[L-2][L-1];
  h[map[L-1][L-1]].Rnbr = map[0][L-1];
  h[map[L-1][L-1]].Dnbr = map[L-1][L-2];
  h[map[L-1][L-1]].Unbr = map[L-1][0];

  ////////////////// sides
  for(i=1;i<L-1;i++){
    h[map[0][i]].Lnbr = map[L-1][i];
    h[map[0][i]].Rnbr = map[1][i];
    h[map[0][i]].Dnbr = map[0][i-1];
    h[map[0][i]].Unbr = map[0][i+1];

    h[map[L-1][i]].Lnbr = map[L-2][i];
    h[map[L-1][i]].Rnbr = map[0][i];
    h[map[L-1][i]].Dnbr = map[L-1][i-1];
    h[map[L-1][i]].Unbr = map[L-1][i+1];

    h[map[i][0]].Lnbr = map[i-1][0];
    h[map[i][0]].Rnbr = map[i+1][0];
    h[map[i][0]].Dnbr = map[i][L-1];
    h[map[i][0]].Unbr = map[i][1];

    h[map[i][L-1]].Lnbr = map[i-1][L-1];
    h[map[i][L-1]].Rnbr = map[i+1][L-1];
    h[map[i][L-1]].Dnbr = map[i][L-2];
    h[map[i][L-1]].Unbr = map[i][0];
  }

  ///////////// interior
  for(i=1;i<L-1;i++){
    for(j=1;j<L-1;j++){
      h[map[i][j]].Lnbr = map[i-1][j];
      h[map[i][j]].Rnbr = map[i+1][j];
      h[map[i][j]].Dnbr = map[i][j-1];
      h[map[i][j]].Unbr = map[i][j+1];
    }
  }
  /////////////////////////////////////
  ////////////////////////////////////	
	

  clock_t start, end;
  double cpuTime;
	
  start = clock();
	
  for(m=0;m<nsmpls;m++){
    for(i=0;i<L*L;i++){
      h[i].height = (int)floor(hini[i]);
      if(uniform64() < hini[i] - h[i].height){
	h[i].height++;
      }
    }
		
    ratesum = 0;
    for(i=0;i<L*L;i++){
      h[i].cord = getcord(h,i);
      h[i].rate = exp(-2*K*h[i].cord);  //  Note:  Rate for an atom to break loose, not break loose and move left.
      ratesum += h[i].rate;
    }
		
		
		
    t = 0;
    while(t<T){
			
      ratesum=0;
      for(i=0;i<L*L;i++)
	ratesum += h[i].rate;
			
      timetonext = -log(uniform64())/ratesum;
			
      //printf("%g\n",eta);
			
      if(t+timetonext>T){
	t = T;
	break;
      }
			
//      for(i=0;i<L*L;i++){
//	ratesint[i] += h[i].rate*timetonext/((double) nsmpls)/((double) pow(L,4));
//      }
			
			
			
      /////////////////// This sampling scheme only requires 1 random variable but scales
      /////////////////// like L so that the whole algorithm scales like L^2
      eta = ratesum*uniform64();
      z2 = 0;
      i = -1;
      while(z2<eta){
	i++;
	z2+= h[i].rate;
      }
      if(i>L*L-1){
	printf("resample overflow\n");
	printf("iteration:  %ld\n", m);
	printf("physical time:  %20.14g\n", t);
	printf("cpu time:  %20.14g\n",(double)((clock()-start)/CLOCKS_PER_SEC));
	//	printf("r1-RAND_MAX:  %20.14g\n",(double)r1-(double)RAND_MAX);
	printf("ratesum:  %20.14g\n", ratesum);
	printf("eta:  %20.14g\n", eta);
	printf("z2:  %20.14g\n", z2);
	return 1;
      }
      if(i<0){
	printf("resample underflow\n");
	printf("iteration:  %ld\n", m);
	printf("physical time:  %20.14g\n", t);
	printf("cpu time:  %20.14g\n",(double)((clock()-start)/CLOCKS_PER_SEC));
	printf("ratesum:  %20.14g\n", ratesum);
	printf("eta:  %20.14g\n", eta);
	printf("z2:  %20.14g\n", z2);
	return 1;
      }
      siteupdatelist[0] = i;
			
			
      ////////////////////This sampling scheme requires lots of random variables but scales
      ////////////////////like 1 so that the whole algorithm scales like L
      //			eta = (double)random()/RAND_MAX;
      //			i = (int)floor(L*(double)random()/RAND_MAX);
      //			while( eta > h[i].rate/r1 ){
      //				eta = (double)random()/RAND_MAX;
      //				i = (int)floor(L*(double)random()/RAND_MAX);
      //			}
      //			siteupdatelist[0] = i;
			
      
      k = (int)floor(4*uniform64());
		
      if (k==0){//////////// jump to the left //////////////
	i1 = h[i].Lnbr;
	siteupdatelist[1] = i1;
	siteupdatelist[2] = h[i].Rnbr;
	siteupdatelist[3] = h[i].Dnbr;
	siteupdatelist[4] = h[i].Unbr;
	siteupdatelist[5] = h[i1].Lnbr;
	siteupdatelist[6] = h[i1].Dnbr;
	siteupdatelist[7] = h[i1].Unbr;
      }
      if (k==1){////////////// jump to the right ////////////
	i1 = h[i].Rnbr;
	siteupdatelist[1] = i1;
	siteupdatelist[2] = h[i].Lnbr;
	siteupdatelist[3] = h[i].Dnbr;
	siteupdatelist[4] = h[i].Unbr;
	siteupdatelist[5] = h[i1].Rnbr;
	siteupdatelist[6] = h[i1].Dnbr;
	siteupdatelist[7] = h[i1].Unbr;
      }
      if (k==2){///////////// jump down ////////////////
	i1 = h[i].Dnbr;
	siteupdatelist[1] = i1;
	siteupdatelist[2] = h[i].Lnbr;
	siteupdatelist[3] = h[i].Rnbr;
	siteupdatelist[4] = h[i].Unbr;
	siteupdatelist[5] = h[i1].Lnbr;
	siteupdatelist[6] = h[i1].Rnbr;
	siteupdatelist[7] = h[i1].Dnbr;
      }
      if (k==3){////////////////// jump up /////////////////
	i1 = h[i].Unbr;
	siteupdatelist[1] = i1;
	siteupdatelist[2] = h[i].Lnbr;
	siteupdatelist[3] = h[i].Rnbr;
	siteupdatelist[4] = h[i].Dnbr;
	siteupdatelist[5] = h[i1].Lnbr;
	siteupdatelist[6] = h[i1].Rnbr;
	siteupdatelist[7] = h[i1].Unbr;
      }
			
			
      h[siteupdatelist[0]].height--;
      h[siteupdatelist[1]].height++;
			
			
			
      for(k=0;k< 8;k++){
      	i2 = siteupdatelist[k];
				
      	oldcord = h[i2].cord;
      	h[i2].cord = getcord(h,i2);
				
      	if(h[i2].cord != oldcord){
      	  oldrate = h[i2].rate;
      	  ratesum -= oldrate;
      	  h[i2].rate = exp(-2*K*h[i2].cord);
      	  ratesum += h[i2].rate;
      	}
      }
			
      /* ratesum = 0; */
      /* for(i=0;i<L*L;i++){ */
      /* 	  h[i].cord = getcord(h,i,p); */
      /* 	  h[i].rate = exp(-2*K*h[i].cord); */
      /* 	  ratesum += h[i].rate; */
      /* } */
			
      t += timetonext;
    }
    for(i=0;i<L*L;i++){
      avh[i] += (double) h[i].height/((double) nsmpls)/((double) L);
      rates[i] += (double) h[i].rate/((double) nsmpls);
      cords[i] += (double) h[i].cord/((double)nsmpls);
    }
  }
	
	
  end = clock();
	
	
  cpuTime= (end-start)/ (CLOCKS_PER_SEC);
	
  printf("\n");
  printf("CPU time: %g minutes\n",cpuTime/60.0);
	
	
  fid = fopen("kmcsmooth2D_p2_K1.5_L150_T1em4_heights.txt","w");
  for(i = 0;i<L;i++){
    for(j=0;j<L;j++){
      fprintf(fid,"%5d\t %5d\t %20.14lf\n",i,j,avh[map[i][j]]);
    }
    fprintf(fid,"\n");
  }
  fclose(fid);
	
  fid = fopen("kmcsmooth2D_p2_K1.5_L150_T1em4_rates.txt","w");
  for(i = 0;i<L;i++){
    for(j=0;j<L;j++){
      fprintf(fid,"%5d\t %5d\t %20.14lf\n",i,j,rates[map[i][j]]);
    }
    fprintf(fid,"\n");
  }
  fclose(fid);
      
	
	
  fid = fopen("kmcsmooth2D_p2_K1.5_L150_T1em4_cords.txt","w");
  for(i = 0;i<L;i++){
    for(j=0;j<L;j++){
      fprintf(fid,"%5d\t %5d\t %20.14lf\n",i,j,cords[map[i][j]]);
    }
    fprintf(fid,"\n");
  }
  fclose(fid);	
	
	
  return 0;
}

