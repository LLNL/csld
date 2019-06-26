#include <stdio.h>
//#include <list>
//using namespace std;

#define EPS 1E-12
#define NOT_ZERO(x) (x >EPS)||(x <-EPS)

// typedef struct {
//   int j;
//   double* v;
// } row_ele;

// typedef struct {
//   slist<row_ele> * plil;
//   double * prow;
// } lil_row;
// typedef list<row_ele> lil_row;

typedef struct {
  int     nrow;
  int     ncol;
  double * val;
  int * nz;
  double ** prow;
  int ** pnz;
} sMatrix;

typedef sMatrix * Matrix;

Matrix NewMatrix(int x_dim, int y_dim)
{
  Matrix m;
  m = new sMatrix;
  
  m->nrow = x_dim;
  m->ncol = y_dim;
  m->val = new double[x_dim * y_dim];
  m->nz  = new int[x_dim * y_dim];
  m->prow = new double*[x_dim];
  m->pnz = new int * [x_dim];
  for(int i=0; i<x_dim; i++) {
    // m->rows[i].prow = m->val+ i*y_dim;
    // m->rows[i].plil = new slist<row_ele>;
    // m->rows[i] = new lil_row;
    m->prow[i] = m->val + i*y_dim;
    m->pnz[i] = m->nz + i*y_dim;
  }
  return m;
}

void FreeMatrix(Matrix m) {
  delete [] m->val;
  delete [] m->nz;
  // for (int i =0; i< m->nrow; i++) {
  //   delete m->rows[i]->plil;
  // }
  delete [] m->prow;
  delete [] m->pnz;
  delete m;
}
 
// void MtxSetRow(Matrix m, int irow, double *v)
// {
//  }
 
Matrix InitMatrix(int x_dim, int y_dim, double* v)
{
  Matrix m;
  m = NewMatrix(x_dim, y_dim);
  double val;
  for (int ix=0; ix<x_dim; ix++) {
    for (int iy=0; iy<y_dim; iy++) {
      int off = ix*y_dim + iy;
      val = v[off];
      if (NOT_ZERO(val)) {
	m->prow[ix][iy] = val;
	m->pnz[ix][iy] = 1;
      } else {
	m->prow[ix][iy] = 0.0;
	m->pnz[ix][iy] = 0;
      }
    }
  }

  return m;
}
 
void MtxDisplay(Matrix m)
{
  const char *sc;
  for (int ix=0; ix<m->nrow; ix++) {
    printf("   ");
    sc = " ";
    for (int iy=0; iy<m->ncol; iy++) {
      printf("%s %6f", sc, m->prow[ix][iy]);
      sc = ",";
    }
    printf("\n");
  }
  printf("\n");
}
 
void MtxMulAndAddRows(Matrix m, int ixrdest, int ixrsrc, double mplr)
{
  double *drow, *srow;
  int * dnz, *snz;
  drow = m->prow[ixrdest];
  srow = m->prow[ixrsrc];
  dnz = m->pnz[ixrdest];
  snz = m->pnz[ixrsrc];
  for (int iy=0; iy<m->ncol; iy++) {
    if (snz[iy]) {
      drow[iy] += mplr * srow[iy];
      if (NOT_ZERO(drow[iy])) {
	dnz[iy] = 1;
      } else {
	dnz[iy] = 0;
	drow[iy] = 0.0;
      }
    }
  }
}
 
void MtxSwapRows(Matrix m, int rix1, int rix2)
{
  if (rix1 == rix2) return;
  double * tmprow;
  int * tmpnz;
  tmprow = m->prow[rix1];
  tmpnz = m->pnz[rix1];
  m->prow[rix1] = m->prow[rix2];
  m->pnz[rix1] = m->pnz[rix2];
  m->prow[rix2] = tmprow;
  m->pnz[rix2] = tmpnz;

  // r1 = m->mtx[rix1];
  // r2 = m->mtx[rix2];
  // for (ix=0; ix<m->dim_x; ix++)
  // {
  //     temp = r1[ix]; 
  // 	r1[ix]=r2[ix]; 
  // 	r2[ix]=temp;
  // }
}
 
void MtxNormalizeRow(Matrix m, int rix, int lead)
{
  // int ix;
  // EL_Type *drow;
  // EL_Type lv;
  // drow = m->mtx[rix];
  // lv = drow[lead];
  // for (ix=0; ix<m->dim_x; ix++) {
  // 	drow[ix] /= lv;
  // }
  double * pr = m->prow[rix];
  int * pnz = m->pnz[rix];
  double invv = 1.0/pr[lead];
  for (int j=0; j<m->ncol; j++) {
    if (pnz[j]) 
      pr[j]*= invv;
  }
}
 
#define MtxGet( m, rix, cix ) m->prow[rix][cix] 

void MtxToReducedREForm(Matrix m)
{
  int lead;
  int rix, iix;
  double lv;
  int rowCount = m->nrow;
 
  lead = 0;
  for (rix=0; rix<rowCount; rix++) {
    if (lead >= m->ncol)
      return;
    iix = rix;
    while (0 == m->pnz[iix][lead]) {
      iix++;
      if (iix == rowCount) {
	iix = rix;
	lead++;
	if (lead == m->ncol)
	  return;
      }
    }
    // MtxSwapRows(m, iix, rix );
    if (iix != rix) {
      double * tmprow;
      int * tmpnz;
      tmprow = m->prow[iix];
      tmpnz = m->pnz[iix];
      m->prow[iix] = m->prow[rix];
      m->pnz[iix] = m->pnz[rix];
      m->prow[rix] = tmprow;
      m->pnz[rix] = tmpnz;
    }
    // MtxNormalizeRow(m, rix, lead );
    double * pr = m->prow[rix];
    int * pnz = m->pnz[rix];
    double invv = 1.0/pr[lead];
    for (int j=0; j<m->ncol; j++)
      if (pnz[j]) pr[j]*= invv;

    for (iix=0; iix<rowCount; iix++) {
      if ( iix != rix ) {
	lv = MtxGet(m, iix, lead );
	MtxMulAndAddRows(m,iix, rix, -lv) ;
      }
    }
    lead++;
  }
}

// Matrix MtxToReducedREForm2(int x_dim, int y_dim, double* v) {
//   Matrix m = InitMatrix(x_dim, y_dim, v);
//   MtxToReducedREForm(m);
//   return m;
// }


int calc_nullspacebasis(Matrix mat, double *b_out) {
  int nfree = 0, nnotf=0;
  int m= mat->nrow;
  int n= mat->ncol;

  int *row, *col;
  row = new int[m];
  for (int i = 0; i < m; i++ ) row[i] = 0;
  col = new int[n];
  for (int j = 0; j < n; j++ ) col[j] = - ( j + 1 );

  // count non-vanishing rows
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      if (mat->pnz[i][j]) {
        row[i] = ( j + 1 );
        col[j] = ( j + 1 );
	nnotf++;
	break;
      }
    }
  }

  nfree = n - nnotf;
  if (nfree <=0 ) return 0;

  //  nullspace = r8mat_zero_new ( n, nullspace_size );
  for (int i=0; i<nfree*n; i++) b_out[i]=0.0;

  int j2 = 0, i2;
//
//  If column J does not contain a leading 1, then it contains
//  information about a null vector.
//
  for (int j = 0; j < n; j++ )
  {
    if ( col[j] < 0 )
    {
      for (int i = 0; i < m; i++ )
      {
        if (mat->pnz[i][j])
        {
          i2 = row[i] - 1;
          b_out[i2+j2*n] = - mat->prow[i][j];
        }
      }
      b_out[j+j2*n] = 1.0;
      j2++;
    }
  }
  delete [] col;
  delete [] row;

/*  for (int i=0; i<nfree; i++) {
    for (int j=0; j<n; j++) printf("%6f ", b_out[i*n+ j]);
    printf("\n");
  }*/

  return nfree;
}

int nullspace(int x_dim, int y_dim, double* v, double *b_out) {
  Matrix m = InitMatrix(x_dim, y_dim, v);
  MtxToReducedREForm(m);
  int nfree = calc_nullspacebasis(m, b_out);
  FreeMatrix(m);
  return nfree;
}


int main()
{
  Matrix m1;
  static double im[][4] = {{1,2,-1,-4},{2,3,-1,-11},{-2,0,-3,22}};
 
  m1 = InitMatrix( 3, 4, (double *)im);
  printf("Initial\n");
  MtxDisplay(m1);
  MtxToReducedREForm(m1);
  printf("Reduced R-E form\n");
  MtxDisplay(m1);
  printf("Null space\n");
  double * bas = new double[4*4];
  calc_nullspacebasis(m1, bas);
  return 0;
}
