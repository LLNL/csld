#include <math.h>
#include <stddef.h>
#define NDIM 3
#define NDIMSQ 9
#define ERRTOL 1e-9
#define MAXNPT 10
//#define DEBUG
#define MAXSGP 48
#ifdef DEBUG
#include <stdio.h>
#include <stdlib.h>
#endif


void intDigits(int * digits, int nInput, int base, int len);

void fct_trans(
		      int npt,
		      int dim,
		      double *gamma,
		      int *pi,
		      double *Gm
		      )
{
  int dimTensor, idx1[MAXNPT], idx2[MAXNPT], i1, i2, ipt, matdim;
  double *GmThis, tmp;
  int dims[] = {0, 0};


  dimTensor= (int) pow(dim, npt);
  matdim=dimTensor*dimTensor;
  dims[0]=dimTensor;
  dims[1]=dimTensor;
//  Gm=new double[matdim];
  for (i1=0; i1< matdim; i1++) Gm[i1]=0.0;
#ifdef DEBUG
  printf("mat dim= %d\n pi= [", matdim);
  for (int i =0; i<npt; i++) printf("%d ", pi[i]);
  printf("]\n");
#endif

  GmThis=Gm;
  for (i1=0; i1 < dimTensor; i1++) {
    intDigits(idx1, i1, dim, npt);
    for (i2=0; i2 < dimTensor; i2++) {
      intDigits(idx2, i2, dim, npt);
      tmp=1.0;
#ifdef DEBUG
      printf("%d %d\n", i1, i2);
#endif
      for (ipt=0; ipt<npt; ipt++) tmp*= gamma[dim*idx1[ipt]+ idx2[pi[ipt]]];
#ifdef DEBUG
      for (ipt=0; ipt<npt; ipt++) printf("%d\n", dim*idx1[ipt]+ idx2[pi[ipt]]);
#endif
      *GmThis = tmp;
      GmThis++;
    }
  }


//  delete [] Gm;
}
