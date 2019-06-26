#include "c_util.h"

#include <cstdio>

template<typename number> int load_file(const char* fname,
                                        number*&     output,
                                        const char* fmt)
{
  unsigned long pos = 0;
  unsigned long len = 0;
  
  FILE* fh = fopen(fname, "r");
  
  if (!fh)
  {
    printf("Couldn't open file '%s'\n", fname);
    return 0;
  }
  
  int count = 0;
  
  printf("Loading '%s'\n", fname);
  
  fseek(fh, 0, SEEK_END);
  len = ftell(fh);
  fseek(fh, 0, SEEK_SET);
  
  while (!feof(fh))
  {
    number v;
    
    pos = ftell(fh);
    if (pos % 100 == 0)
    {
      printf("\r  Scanning....%3%%", (pos/(float)len)*100.0);
      fflush(stdout);
    }
    
    if (fscanf(fh, fmt, &v) == 1)
    {
      count++; }
  }
  
  printf("\r  Scanning....Done!\n", (pos/(float)len)*100.0);
  
  printf("  Found %d entries\n", count);
  
  output = new number[count];
  
  fseek(fh, 0, SEEK_SET);
  
  count = 0;
  pos   = 0;
  while (!feof(fh))
  {
    number v;
    
    pos = ftell(fh);
    if (pos % 100 == 0)
    {
      printf("\r  Loading....%3%%", (pos/(float)len)*100.0);
      fflush(stdout);
    }
    
    if (fscanf(fh, fmt, &v) == 1)
    {
      output[count] = (number)v;
      count++;
    }
  }
  
  printf("\r  Loading....Done!", (pos/(float)len)*100.0);
  
  fclose(fh);
  
  printf("\n\n");
  
  return count;
  
};

#ifdef __MAIN__
int main(int argc, char** argv)
{
  double* dx;         int dxN;
  int*    clusALLpts; int clusALLptsN;
  int*    clusALLunq; int nClus;
  int*    clusALLsym; int clusALLsymN;
  int*    uniqueNpt;  int nUnique;
  double* uniqueFac;  int uniqueFacN;
  double* sgopALLmat; int sgopALLmatN;
  
  double* A;
  int maxnpt;
  
  
  if (argc < 2)
  {
    printf("Usage: %s [maxnpt]\n", argv[0]);
    return -1;
  }
  
  sscanf(argv[1], "%d", &maxnpt);
  
  dxN         = load_file<double>("data/dx_data",         dx,         "%lf");
  clusALLptsN = load_file<int>("data/clusALLpts_data",    clusALLpts, "%d");
  clusALLsymN = load_file<int>("data/clusALLsym_data",    clusALLsym, "%d");
  nClus       = load_file<int>("data/clusALLunq_data",    clusALLunq, "%d");
  sgopALLmatN = load_file<double>("data/sgopALLmat_data", sgopALLmat, "%lf");
  nUnique     = load_file<int>("data/uniqueNpt_data",     uniqueNpt,  "%d");
  uniqueFacN  = load_file<double>("data/uniqueFac_data",  uniqueFac,  "%lf");
  
  A = new double[(dxN+1) * 1417];
  
  for (int i=0; i < 100; i++)
  {
    ForceCorrML(dx, dxN,
                clusALLpts, clusALLptsN,
                clusALLunq, nClus,
                clusALLsym, clusALLsymN,
                uniqueNpt,  nUnique,
                uniqueFac,  uniqueFacN,
                sgopALLmat, sgopALLmatN,
                maxnpt,
                A);
  }
  
  
  FILE* fh = fopen("A.dat", "w");
  for (int i=0; i < (dxN+1) * 1417; i++)
  {
    fprintf(fh, "%g\n", A[i]);
  }
  fclose(fh);
  
  return 0;
};
#endif
