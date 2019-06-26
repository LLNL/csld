"""
C code for csld, including
fct_trans
"""

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np
import scipy.sparse

# if you want to use the Numpy-C-API from Cython
# (not strictly necessary for this example)
np.import_array()

# ctypedef np.int_t DTYPE_t
# DTYPE = np.int


# cdefine the signature of our c function
cdef extern from "c_util.h":
    void fct_trans(int npt, int dim, double *gamma, int *pi, double *Gm)

# create the wrapper code, with numpy type annotations
def fct_trans_c(npt, dim, np.ndarray[np.double_t, ndim=2, mode="c"] gamma not None,
                     pi_in not None):
    cdef np.ndarray[np.double_t, ndim=2] GM = np.zeros((dim**npt, dim**npt))
    cdef np.ndarray[np.int32_t, ndim=1] pi = np.array(pi_in, dtype=np.int32)
    fct_trans(npt, dim, <double*> gamma.data, <int*> pi.data, <double*> GM.data)

    return scipy.sparse.dok_matrix(GM)





# force+energy correlation
cdef extern from "c_util.h":
    void ForceCorrML(
                 double *dx, int dxN,
                 int *clusALLpts, int clusALLptsN,
                 int *clusALLunq, int nClus,
                 int *clusALLsym, int clusALLsymN,
                 int *uniqueNpt, int nUnique,
                 double *uniqueFac, int uniqueFacN,
                 double *sgopALLmat, int sgopALLmatN,
                 int maxnpt,
                 double *A)

# create the wrapper code, with numpy type annotations
def ld_get_correlation(natom, norbit, maxnpt, nfct_tot,
        np.ndarray[np.double_t, ndim=2, mode="c"] dx not None,
        np.ndarray[np.int32_t, ndim=2, mode="c"] clusALLpts not None,
        np.ndarray[np.int32_t, ndim=1, mode="c"] clusALLunq not None,
        np.ndarray[np.int32_t, ndim=1, mode="c"] clusALLsym not None,
        np.ndarray[np.int32_t, ndim=1, mode="c"] uniqueNpt not None,
        np.ndarray[np.double_t, ndim=1, mode="c"] uniqueFac not None,
        np.ndarray[np.double_t, ndim=3, mode="c"] sgopALLmat not None):

    cdef np.ndarray[np.double_t, ndim=2] A = np.zeros((dx.size+1, nfct_tot))
    ForceCorrML(<double*> dx.data, dx.size,
                <int *> clusALLpts.data, clusALLpts.size,
                <int *> clusALLunq.data, clusALLunq.size,
                <int *> clusALLsym.data, clusALLsym.size,
                <int *> uniqueNpt.data, uniqueNpt.size,
                <double *> uniqueFac.data, uniqueFac.size,
                <double*> sgopALLmat.data, sgopALLmat.size,
                maxnpt,
                <double*> A.data)

    return A



cdef extern from "c_util.h":
    int nullspace(int x_dim, int y_dim, double* v, double *b_out)

def get_nullspace(sp_mat):
    (m, n) = sp_mat.shape
    cdef np.ndarray[np.double_t, ndim=1] b_out = np.zeros(n* n)
    cdef np.ndarray[np.double_t, ndim=2] mat = sp_mat.todense()
    nfree = nullspace(m, n, <double*> mat.data, <double*> b_out.data)
    if nfree<=0:
        return scipy.sparse.csr_matrix((0, n))
    else:
        # cdef np.ndarray[np.double_t, ndim=1] c_out = b_out[:nfree*n].copy()
        # print(b_out.shape)
        return scipy.sparse.csr_matrix(b_out[:nfree*n].copy().reshape((nfree, n)))



cdef extern from "c_util.h":
    void init_basis(int tgt, int lmax, int npt, const double *x, const double *y)
    void calc_val(double x, int lmax, double *val)
    void calc_der(double x, int lmax, double *val)
    void LDFFCorr(
      int *orbNpt, int Nunique,
      int *lmaxlist, int Nlmax,
      int *ncorrlist, int Nncorrlist,
      int *multilist, int NCorr,
      int *ffidxlist, int Nffidxlist,
      double *flatOrbPos, int NflatOrbPos,
      double *dx, int dxN,
      int *allTyp, int NallTyp,
      int *allFlatOrbIdx, int NC,
      int *allClus, int NallClus,
      double * A)

# create the wrapper code, with numpy type annotations
def init_ldff_basis(tgt, lmax, np.ndarray[np.double_t, ndim=1, mode="c"] x not None,
    np.ndarray[np.double_t, ndim=2, mode="c"] y not None):
    init_basis(tgt, lmax, x.size, <double *> x.data, <double *> y.data)

def get_val(x, lmax):
    cdef np.ndarray[np.double_t, ndim=1] val = np.zeros(lmax)
    calc_val(<double> x, <int> lmax, <double*> val.data)
    return val

def get_der(x, lmax):
    cdef np.ndarray[np.double_t, ndim=1] val = np.zeros(lmax)
    calc_der(<double> x, <int> lmax, <double*> val.data)
    return val

def ldff_get_corr(
        np.ndarray[np.int32_t, ndim=1, mode="c"] orbNpt not None,
        np.ndarray[np.int32_t, ndim=1, mode="c"] lmaxlist not None,
        np.ndarray[np.int32_t, ndim=1, mode="c"] ncorrlist not None,
        np.ndarray[np.int32_t, ndim=1, mode="c"] multilist not None,
        np.ndarray[np.int32_t, ndim=1, mode="c"] ffidxlist not None,
        np.ndarray[np.double_t, ndim=3, mode="c"] flatOrbPos not None,
        np.ndarray[np.double_t, ndim=2, mode="c"] dx not None,
        np.ndarray[np.int32_t, ndim=1, mode="c"] allTyp not None,
        np.ndarray[np.int32_t, ndim=1, mode="c"] allFlatOrbIdx not None,
        np.ndarray[np.int32_t, ndim=2, mode="c"] allClus not None):
    cdef np.ndarray[np.double_t, ndim=2] A = np.zeros((dx.size+1, multilist.size))
    LDFFCorr(<int *> orbNpt.data, orbNpt.size,
             <int *> lmaxlist.data, lmaxlist.size,
             <int *> ncorrlist.data, ncorrlist.size,
             <int *> multilist.data, multilist.size,
             <int *> ffidxlist.data, ffidxlist.size,
             <double *> flatOrbPos.data, flatOrbPos.size,
             <double *> dx.data, dx.size,
             <int *> allTyp.data, allTyp.size,
             <int *> allFlatOrbIdx.data, allFlatOrbIdx.size,
             <int *> allClus.data, allClus.size,
             <double*> A.data)
    return A


cdef extern from "c_util.h":
    void structure_ordering(int dim, int n1, double *p1, int n2, double *p2, int p1noduplicate, double tol, int *pi);

def get_structure_ordering(np.ndarray[np.double_t, ndim=2, mode="c"] p1 not None,
                     np.ndarray[np.double_t, ndim=2, mode="c"] p2 not None, p1noduplicate=1, tol=1E-4):
#    print("p1,2 shape", p1.__class__, p2.__class__, p1.shape[0],p1.shape[1],p2.shape[0],p2.shape[1])
    cdef np.ndarray[np.int32_t, ndim=1] pi = np.zeros(p1.shape[0], dtype=np.int32)
    structure_ordering(p1.shape[1], p1.shape[0], <double*> p1.data, p2.shape[0], <double*> p2.data, p1noduplicate, <double> tol, <int*> pi.data)
    return pi


