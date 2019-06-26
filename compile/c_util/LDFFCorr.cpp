#include <math.h>
#include "spline.h"
#include <stdio.h>
// #include <stddef.h>

#define NDIM 3
#define MAXNPT 4
#define MAXPAIR 6
#define ERRTOL 1e-9
#define MAXSGP 48
#define MAXL 50

// #define DEBUG
#ifdef DEBUG
#endif

/* {dx, _Real, 2}, {clusPts, _Integer, 2}, {clusOrb, _Integer, 1}, {clusSym, _Integer, 1}, 
    LDFFCorr::usage = "C++ code to calculate correlation of lattice dynamics force field given supercell";
*/

inline double vlen(double *r) {
  return sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
}


tk::spline *b2=NULL, *b3=NULL;

void init_basis(int tgt, int lmax, int npt, const double *x, const double *y)
{
  tk::spline * b;
  if (tgt == 2) {
    if (b2) delete [] b2;
    b2 = new tk::spline[lmax];
    b = b2;
  } else {
    if (b3) delete [] b3;
    b3 = new tk::spline[lmax];
    b = b3;
  }
   
   std::vector<double> X(x, x+npt);
   for (int l=0; l<lmax; l++) {
     std::vector<double> Y(y+l*npt, y+(l+1)*npt);
     b[l].set_points(X,  Y);    // currently it is required that X is already sorted
   }
}

void calc_val(double x, int lmax, double *val) {
  for (int il=0; il<lmax; il++)
    val[il] = b2[il](x);
}

void calc_der(double x, int lmax, double *val) {
  for (int il=0; il<lmax; il++)
    val[il] = b2[il].derivative(x);
}

/*void calc_bas_val(char *typ, double *x_in, int Nx, double **val_vec, int *len) {
  long chk_packet;
  
  MLPutFunction(stdlink, "EvaluatePacket", 1);
  MLPutFunction(stdlink, typ, 1);
  MLPutReal64List(stdlink, x_in, Nx);
  MLEndPacket(stdlink);
  MLCheckFunction(stdlink, "ReturnPacket", &chk_packet);
  MLGetReal64List(stdlink, val_vec, len);
  //  printf("input len = %d output = %d\n", Nx, *len);
  //  for (int i=0; i< *len; i++) {printf("%f , ", (*val_vec)[i]);}printf("\n");

}*/



void LDFFCorr( 
	      int *orbNpt, int Nunique, // number of points in distinct clusters or orbits; # of orbits e.g. {2,2,2} for 3 pair interactions
	      int *lmaxlist, int Nlmax, // for npt = 2, 3, 4; Nlmax==3
	      int *ncorrlist, int Nncorrlist, // number of correlation functions for each cluster/orbit; Nunique = Nncorrlist
	      int *multilist, int NCorr, // NCorr = Total[ncorrlist]
	      int *ffidxlist, int Nffidxlist,
	      double *flatOrbPos, int NflatOrbPos, // [# of flattened orbits, MAXNPT, 3] in cartesian coordinates
	      double *dx, int dxN, // dxN = 3* nAtom
	      int *allTyp, int NallTyp, // which orbit
	      int *allFlatOrbIdx, int NC, // NC = NallTyp = # of actual, translated clusters in supercell
	      int *allClus, int NallClus, // [NC, MAXNPT]
              double * A)
{
  int nAtom, npt, maxnpt, npair, pairs[MAXPAIR*2];
  double r0[MAXPAIR], r[MAXPAIR], dr[MAXPAIR], rvec[MAXPAIR*NDIM];
  double br[MAXPAIR*MAXL], dbr[MAXPAIR*MAXL];
  int nbr, ndbr;
  double *posThis, val, val1;
  int dim[] = {0, 0}, *clusThis;
  int *startidx, *startff, *ffThis, ifree, icorr, ipos, iC, lmax, ipair;
 

#ifdef DEBUG
  printf("started\n");
#endif


  if (Nunique!= Nncorrlist) {
    printf("Nunique %d !=Nncorrlist %d\n", Nunique, Nncorrlist);
    exit(-1);
  }

  if (Nlmax!= 3) {
    printf("Nlmax %d !=3\n", Nlmax);
    exit(-1);
  }

  if (NC!= NallTyp) {
    printf("NC %d !=NallTyp %d\n", NC, NallTyp);
    exit(-1);
  }

  if (NallClus!= NC*MAXNPT) {
    printf("NallClus %d != NC %d * MAXNPT %d\n", NallClus, NC, MAXNPT);
    exit(-1);
  }


  nAtom= dxN/NDIM;
  maxnpt= NallClus/NC;
  dim[0]= dxN;
  startidx= new int[Nunique];
  startidx[0]=0;
  for(int i=0; i<Nunique-1; i++) startidx[i+1]=startidx[i] + ncorrlist[i];
  dim[1]= startidx[Nunique-1] + ncorrlist[Nunique-1];
  if (dim[1]!=NCorr) {
    printf("dim[1] %d != NCorr %d\n", dim[1], NCorr);
    exit(-1);
  }
  startff= new int[Nunique];
  startff[0]=0;
  icorr=0;
  for(int i=0; i<Nunique-1; i++) {
    startff[i+1]= startidx[i];
    for (int j=0; j<ncorrlist[i]; j++, icorr++) {
      startff[i+1] += multilist[icorr]*(orbNpt[i])*(orbNpt[i]-1)/2;
    }
  }

  // append total energy correlation (one row) at the end
  dim[0]++;

  // A=new double[dim[0]*dim[1]];
  for (int i=0; i<dim[0]*dim[1]; i++) A[i]=0.0;

#ifdef DEBUG
  printf("A initialized\n");
#endif


  for(int iclus=0;  iclus< NC; iclus++) {
    iC=allTyp[iclus];
    npt = orbNpt[iC];
    if (npt> MAXNPT) {
      printf("npt %d  >MaxNpt %d", npt, MAXNPT);
      exit(-1);
    }
    npair = npt*(npt-1)/2;
    ipos=allFlatOrbIdx[iclus];
    lmax= lmaxlist[npt-2];
    posThis= flatOrbPos + ipos*MAXNPT*3;
    clusThis= allClus+ iclus*MAXNPT;
    ffThis= ffidxlist+ startff[iC];

#ifdef DEBUG
    printf("iclus %d, iC %d, npt %d, npair %d, lmax %d\n", iclus, iC, npt, npair, lmax);
#endif


    ipair=0;
    for (int ia=0; ia<npt; ia++) {
      for (int ib=ia+1; ib<npt; ib++) {
	pairs[ipair*2+0] = clusThis[ia];
	pairs[ipair*2+1] = clusThis[ib];
	
	for (int ix=0; ix<NDIM; ix++) rvec[ipair*NDIM+ix] = posThis[ib*NDIM+ix]-posThis[ia*NDIM+ix];
	r0[ipair] = vlen(rvec+ipair*NDIM);
	for (int ix=0; ix<NDIM; ix++) rvec[ipair*NDIM+ix]+= dx[ pairs[ipair*2+1]*NDIM+ix]-dx[pairs[ipair*2+0]*NDIM+ix];
	r[ipair]  = vlen(rvec+ipair*NDIM);
	dr[ipair]= r[ipair]-r0[ipair];
	
	ipair++;
      }
    }

    /* if (npair==1) {
      calc_bas_val("b2", dr, npair, &br,  &nbr);
      calc_bas_val("d2", dr, npair, &dbr, &ndbr);
    } else {
      calc_bas_val("b3", dr, npair, &br,  &nbr);
      calc_bas_val("d3", dr, npair, &dbr, &ndbr);
    } */


    tk::spline *b;
    if (npair==1) {
      b=b2;
    } else {
      b=b3;
    }
    nbr = lmax*npair;
    ndbr = lmax*npair;
    for (int inb=0; inb<npair; inb++) {
      for (int il=0; il<lmax; il++) {
	br[inb*lmax + il] = b2[il](dr[inb]);
	dbr[inb*lmax + il] = b2[il].derivative(dr[inb]);
      }
    }
#ifdef DEBUG
    printf("obtained from Kernel nbr=%d ndbr=%d\n", nbr, ndbr);
#endif
 
    for(int ii=0; ii<ncorrlist[iC]; ii++){
#ifdef DEBUG
      printf("  ii= %d\n", ii);
#endif
      icorr= startidx[iC]+ii;
      for (int imulti=0; imulti<multilist[icorr]; imulti++, ffThis+=npair) {
	// b_l1(r_1) ... b'_iq() .. b_ln(r_n)
	// calculate force correlation
	for (int p=0; p<npair; p++) {
	  val= 1.0/r[p];
	  int iA= pairs[p*2+0], iB= pairs[p*2+1];
			
	  for (int q=0; q<npair; q++) {

	    if (p==q) {
	      val*= dbr[q*lmax+ ffThis[q]];
	    } else {
	      val*=  br[q*lmax+ ffThis[q]];
	    } 	      
	  }
	    
	  for(int ix=0; ix<NDIM; ix++) {
	    val1= val*rvec[p*NDIM+ix];
	    ifree= NDIM*iA+ix;
	    A[ifree*dim[1]+ icorr]+= val1;
	    ifree= NDIM*iB+ix;
	    A[ifree*dim[1]+ icorr]-= val1;
	  }
	} // end for p

	// calculate total energy correlation
	val= 1.0;
	for (int q=0; q<npair; q++) {
	  val*=  br[q*lmax+ ffThis[q]];
	}
	A[((int) dxN)*dim[1]+ icorr] += val;
#ifdef DEBUG
	printf("iclus= %d val= %g\n", iclus, val);
#endif
      }
    }

  }


  delete [] startidx;
  delete [] startff;

  /* finally return CorrMat */
//  MLPutReal64Array(stdlink, A, dim, NULL, 2);
//  delete [] A;


}


#ifdef THIS_IS_MAIN
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "spline.h"

int main(int argc, char** argv) {

   std::vector<double> X(5), Y(5);
   X[0]=0.1; X[1]=0.4; X[2]=1.2; X[3]=1.8; X[4]=2.0;
   Y[0]=0.1; Y[1]=0.7; Y[2]=0.6; Y[3]=1.1; Y[4]=0.9;

   tk::spline s;
   s.set_points(X,Y);    // currently it is required that X is already sorted

   double x=1.5;
   printf("spline at %f is %f\n", x, s(x));

   int l=1;
   init_basis(2, l, 5, X.data(), Y.data());
   printf("b2(%d, %f)= %f\n", 0, x, b2[0](x));

   double YY[] = {-100, -99, -96, -91, -80, 1, 2, 3, 4, 6};
   init_basis(2, 2, 5, X.data(), YY);
   for (x=-1.; x<3.; x+=0.1)
     printf("b2(%d, %f )= %f\n", 1, x, b2[0](x));
   

   return EXIT_SUCCESS;
}
#endif
