#include <math.h>
//#include <stddef.h>
#define NDIM 3
#define NDIMSQ 9
#define ERRTOL 1e-9
//#define MAXNPT 10
#define MAXSGP 48
// #define DEBUG
#ifdef DEBUG
#include <stdio.h>
#endif

void intDigits(int * digits, int nInput, int base, int len);

void ForceCorrML( 
		 double *dx, int dxN,
		 int *clusALLpts, int clusALLptsN,
		 int *clusALLunq, int nClus,
		 int *clusALLsym, int clusALLsymN,
		 int *uniqueNpt, int nUnique,
		 double *uniqueFac, int uniqueFacN,
		 double *sgopALLmat, int sgopALLmatN,
		 int maxnpt,
		 double *A
		  )
{
  int dimTensor, *dimTensorList, nAtom, npt;
  double tmp2, fac, *sgopThis, *sgopDot, *dxDot;
  int dim[] = {0, 0};
  int  iA, ix, iC, p, q, iclus, ifree, ii, *idx, *clusThis;
  double * gammadx; // pre-compute gamma.u for all sgp matrices
  int Nsgp, i1, i2, *startidx, ifree0, isym, icorr;

  nAtom= dxN/NDIM;
  maxnpt= clusALLptsN/nClus;
  dimTensorList=new int[nUnique];
  idx=new int[maxnpt];
  dim[0]= dxN;
  dim[1]=0;
  for(iC=0; iC<nUnique; iC++) {
    dimTensorList[iC]= (int) pow(NDIM, uniqueNpt[iC]);
    dim[1]+= dimTensorList[iC];
  }
  startidx = new int[nUnique];
  startidx[0]=0;
  for(iC=0; iC<nUnique-1; iC++) startidx[iC+1]=startidx[iC] + dimTensorList[iC];

  // append total energy correlation (one row) at the end
  dim[0]++;

  //  A=new double[dim[0]*dim[1]];
  for (ix=0; ix<dim[0]*dim[1]; ix++) A[ix]=0.0;

    Nsgp=sgopALLmatN/NDIMSQ;
  if(Nsgp>MAXSGP){
#ifdef DEBUG
    printf("No. of sgop %d > max sgp= %d\n",Nsgp, MAXSGP);
#endif
    return;
  }
  gammadx= new double[dxN*Nsgp];
  for (ix=0; ix<dxN*Nsgp; ix++) gammadx[ix]=0.0;
#ifdef DEBUG
  printf("gamma.dx allocated with Nsgp=%d size=%d\n",Nsgp,(int) dxN*Nsgp);
  printf("input u[0]=%lf, %lf, %lf\n",dx[0],dx[1],dx[2]);
  printf("input u[1]=%lf, %lf, %lf\n",dx[3+0],dx[3+1],dx[3+2]);
#endif
  for(iC=0; iC<Nsgp; iC++) {
    for (ix=0; ix<nAtom; ix++) {
      for (i1=0; i1<NDIM; i1++){
	for (i2=0; i2<NDIM; i2++){
	  gammadx[iC*dxN+ NDIM*ix+i1]+= sgopALLmat[iC*NDIMSQ+i1*NDIM+i2]*dx[NDIM*ix+i2];
	}
      }
    }
  }
#ifdef DEBUG
  printf("gamma.dx initialized\n");
#endif


  for(iclus=0;  iclus< nClus; iclus++) {
    iC=clusALLunq[iclus];
    isym=clusALLsym[iclus];
    npt = uniqueNpt[iC];
    dimTensor = dimTensorList[iC];
    fac = -1.0/uniqueFac[iC]; // note -1/factorial
    clusThis= clusALLpts+ iclus*maxnpt;
    sgopThis= sgopALLmat+ isym*NDIMSQ;
    dxDot= gammadx+ isym*dxN;
    for(ii=0; ii<dimTensor; ii++){
      intDigits(idx, ii, NDIM, npt);
      icorr= startidx[iC]+ii;

      // calculate force correlation
      for(q=0; q<npt; q++) {
	/* take derivative of the q-th atom in cluster */
	iA=clusThis[q];
	ifree0= iA*NDIM;
	
	sgopDot=sgopThis+idx[q]*NDIM;
	tmp2 = fac;
	for(p=0; p<npt; p++){
	  if(p == q) continue;
	  tmp2 *= dxDot[clusThis[p]*NDIM+idx[p]];
#ifdef DEBUG
	  //	  fout << "ifree0, icorr, tmp2,q,p= " << ifree0<<" "<<icorr<<" "<<tmp2<< " "<<q <<" "<<p << endl;
#endif
	}
	for(ix=0; ix<NDIM; ix++) {
	  ifree= ifree0+ix;
	  A[ifree*dim[1]+ icorr]+= tmp2* sgopDot[ix];
	}
      } // end for q
      
        // calculate total energy correlation
      tmp2= -fac;
      for (p=0; p<npt; p++) {
	tmp2 *= dxDot[clusThis[p]*NDIM+idx[p]];
      }

#ifdef DEBUG
      if (fabs(tmp2)> 1E-10)
	printf(" calc tot En npt=%d, p=%d, tmp2=%lf, icorr=%d, clus=%d %d, idx= %d %d\n", npt, p, tmp2, icorr, clusThis[0], clusThis[1], idx[0], idx[1]);
#endif

      A[((int) dxN)*dim[1]+ icorr]+= tmp2;
    }
  }
      
#ifdef DEBUG
      printf("done\n");
#endif

      delete [] dimTensorList;
      delete [] idx;
      delete [] gammadx;
      delete [] startidx;

  //  delete [] A;
}



