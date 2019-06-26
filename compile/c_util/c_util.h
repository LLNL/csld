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
                  );
void init_basis(int tgt, int lmax, int npt, const double *x, const double *y);
void calc_val(double x, int lmax, double *val);
void calc_der(double x, int lmax, double *val);

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
              double * A);

void fct_trans(int npt, int dim, double *gamma, int *pi, double *Gm);

int nullspace(int x_dim, int y_dim, double* v, double *b_out);

void structure_ordering(int dim, int n1, double *p1, int n2, double *p2, int p1noduplicate, double tol, int *pi);

void CECorr(
                 double *sigma, int nAtom,
                 int *clusALLpts, int clusALLptsN,
                 int *clusALLunq, int nClus,
                 int *uniqueNpt, int nUnique,
                 int maxnpt,
                 double *A
                  );

void magphocorr( 
                 double *dx, int dxN,
                 int *clusLDpts, int clusLDptsN,
                 int *clusLDsym, int nClus,
                 double *uniqueFac, int uniqueFacN,
                 double *sgopALLmat, int sgopALLmatN,
                 // above 6 lines for the LD cluster
                 int maxnpt, int maxnptMAG,
                 double *vec, int vecN,
                 int nDOF, // degrees of freedom; should be 2 
                 int *orbits, int nOrbit,
                 int *orders_of_orbit, int nOrbit2,
                 int *singleorb_orbit, int nOrbit3,
                 // overall orbit id of composite cluster
                 int *clusMAGpts, int clusMAGptsN,
                 double m0,
                 // above 2 lines for the MAG cluster
                 // WARNING: so far assuming magnetic interaction is v1.v2 * v3.v4 * ...
                 // i.e. ignoring other permutations of the sites, for simplicity
                 double *A
                  );

