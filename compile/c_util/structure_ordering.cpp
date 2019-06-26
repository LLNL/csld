#include <math.h>
//#define DEBUG
#ifdef DEBUG
#include <stdio.h>
#include <stdlib.h>
#endif


void structure_ordering(int dim, int n1, double *p1, int n2, double *p2, int p1noduplicate, double tol, int *pi) {
  int * found = new int[n2]();
  double dr, *fcrd, min_dr, tmp;
  int is_small;

  fcrd = p1;
  for (int i =0; i<n1; i++) {
    // if tol > 0, then two atoms match only when they are close within tol
    // if tol < 0, then select the closest j


    // NOTE: if any element of pi is negative, the search has failed!
    pi[i] = -1;
    for (int j=0; j<n2; j++) {
      min_dr = 1E99;
      if (p1noduplicate && found[j]) continue;
      if (tol>0) {
        is_small = 1;
        for (int k=0; k<dim; k++) {
          dr = fcrd[k] - p2[j*dim+k];
          dr -= round(dr);
          if (fabs(dr)>tol) {
            is_small = 0;
            break;
          }
        }
        if (is_small) {
          pi[i] = j;
          found[j] = 1;
          break;
        }
      } else {
        dr = 0;
        for (int k=0; k<dim; k++) {
            tmp = fcrd[k] - p2[j*dim+k];
            tmp-=round(tmp);
            dr+= tmp*tmp;
        }
        dr=sqrt(dr);
        if (dr < min_dr) {
          pi[i] = j;
          min_dr = dr;
        }
      }
    }
    if (tol<=0) found[pi[i]] = 1;
#ifdef DEBUG
    printf("%d: %d ", i, pi[i]);
#endif
    fcrd+= dim;
  }
#ifdef DEBUG
  printf("\n");
#endif

  delete [] found;
}

