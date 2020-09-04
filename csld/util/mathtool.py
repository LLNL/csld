import math
import numpy as np
import itertools
import scipy, scipy.sparse, scipy.sparse.linalg
from scipy.sparse import lil_matrix as spmat
# from .rref import rref
#from .cy_rref import cy_rref
from _c_util import fct_trans_c, get_nullspace


def cofactor(M):
    return np.linalg.det(M)* np.linalg.inv(M).T

def vec_linspace(ls, n, remove_duplicate_middle_point=False):
    if remove_duplicate_middle_point:
        return np.vstack([np.array([np.linspace(i,j,n) for i,j in zip(ls[i], ls[i+1])]).T[0 if i<=0 else 1:] for i in range(len(ls)-1)])
    else:
        return np.vstack([np.array([np.linspace(i,j,n) for i,j in zip(ls[i], ls[i+1])]).T for i in range(len(ls)-1)])

#  Normalize[RandomVariate[NormalDistribution[0,1],DIM ]]RandomReal[{rmin,rmax}]
def random_displacements(natom, rmin, rmax=None, dim=3):
    dx = np.random.normal(size=(natom, dim))
    veclen = np.full(natom, rmin) if rmax is None else np.random.uniform(rmin, rmax, natom)
    dx*= (veclen/np.linalg.norm(dx, axis=1))[:,None]
    return dx

def entropy_prob(p_in, scaled=True):
    """
    :param p: probability
    scaled: whether to return entropy/ln(N), i.e. between 0 & 1
    return entropy
    """
    p = p_in/np.sum(p_in)
    s = -np.dot(p, np.log(p))
    return s/np.log(len(p)) if scaled else s

def multi_index_fac(list1):
    """
    multi-index factorial 
    """
    list2=[math.factorial(i) for i in list1]
    return np.prod(list2)

def remove_duplicates(values, compare):
    """
    remove duplicate, keeping order of original elements, using CUSTOM comparison method
    :param values: input list
    :param compare: functions to compare the values
    :return: values without duplicate
    """
    output = []
    for value in values:
        # If value has not been encountered yet,
        # ... add it to output
        if value in output:
            continue
        match = False
        for v2 in output:
            if compare(value, v2):
                match= True
                break

        if not match:
            output.append(value)
    return output


def allindex(value, inlist):
    """
    Mathematica Positions -equivalent
    :param value:
    :param inlist: list from which to find value
    :return: all indices
    """
    # indices = []
    # idx = -1
    # while True:
    #     try:
    #         idx = qlist.index(value, idx+1)
    #         indices.append(idx)
    #     except ValueError:
    #         break
    # return indices
    return [i for i, x in enumerate(inlist) if x == value]


def relativePosition1d(origMappedOne2One, final):
    """
    indexing function such that
    map of pi(i) th element in original list = i-th list in transformed list
    :param origMappedOne2One:
    :param final:
    :return:
    """
    resl=[]
    for i in final:
        resl.append(allindex(i, origMappedOne2One))

    for i in range(len(final)):
        flag=len(resl[i])>1
        resl[i]=resl[i][0]
        if flag:
            for j in range(i+1,len(final)):
                resl[j]= [x for  x in resl[j] if x != resl[i]]
    return resl


def relativePosition(origMappedOne2One, final, dim=1):
    """
    If input is 2d (list of sublists), then return results for each sublist (listable)
    :param origMappedOne2One:
    :param final:
    :return:
    """
    if dim>1:
        return [relativePosition1d(origMappedOne2One[i], final[i]) for i in range(len(final))]
    else:
        return relativePosition1d(origMappedOne2One, final)


#define general tally
def Tally(listin):
    listdif=[]
    for i in listin:
        if not i in listdif:
            listdif.append(i)
    counter=[]
    for i in listdif:
        cot=0
        for j in listin:
            if i==j:
                cot+=1
        counter.append(cot)
    listout=[]
    for i in range(len(counter)):
        listout.append([listdif[i],counter[i]])
    return listout


def Union(listin, return_index=False):
    listout=[]
    idxout=[]
    for i, val in enumerate(listin.tolist() if isinstance(listin, np.ndarray) else listin):
        if not val in listout:
            listout.append(val)
            idxout.append(i)
    return idxout if return_index else listout

def Position(listin,el):
    tmp=[]
    for i in range(len(listin)):
        if listin[i]==el:
            tmp.append(i+1)
    return tmp

def ListFlat(listin):
    # tmp=[]
    # for i in listin:
    #     tmp=tmp+i
    # return tmp
    return list(itertools.chain.from_iterable(listin))

#List
#input: list=[1,2,3]
#output: listlist=[[1],[2],[3]]
def List(listin):
    listout=[]
    for i in listin:
        listout.append([i])
    return listout

def IntegerDigits(n, b, len_digits):
    """
    define IntegerDigits same as in Mathematica
    :param n: integer
    :param b: base
    :param len_digits:
    :return: list of integer base b with paddled length equals len_digits
    """
    # tmp=base10ton(n, b)
    # tmp=[int(i) for i in str(tmp)]
    # assert len(tmp)<=lg, "Wrong with base transformation"
    # if len(tmp)==lg:
    #     return np.array(tmp)
    # else:
    #     return np.array([0 for _ in range(lg-len(tmp))] + tmp)
    if len_digits==0:
        return np.array([])
    digits = np.base_repr(n, 3, padding=10)[-len_digits:]
    return np.array([int(i) for i in digits])

# def base10ton(num, base):
#     """Change ``num'' to given base
#     Upto base 36 is supported."""
#
#     converted_string, modstring = "", ""
#     currentnum = num
#     if not 1 < base < 37:
#         raise ValueError("base must be between 2 and 36")
#     if not num:
#         return '0'
#     while currentnum:
#         mod = currentnum % base
#         currentnum = currentnum // base
#         converted_string = chr(48 + mod + 7*(mod > 10)) + converted_string
#     return converted_string


#define Permutations same as in Mathematica
#input: list elements
#output: unioned permutations of the element
#def Permutations(listin):
def perm(n, i):
    if i == len(n) - 1:
        return n
    else:
        for j in range(i, len(n)):
            n[i], n[j] = n[j], n[i]
            perm(n, i + 1)
            n[i], n[j] = n[j], n[i] # swap back, for the next loop
#http://www.daniweb.com/software-development/python/code/216696/a-simple-recursive-permutation-function-python
def permutate(seq):
    """permutate a sequence and return a list of the permutations"""
    if not seq:
        return [seq]  # is an empty sequence
    else:
        temp = []
        for k in range(len(seq)):
            part = seq[:k] + seq[k+1:]
            #print k, part  # test
            for m in permutate(part):
                temp.append(seq[k:k+1] + m)
                #print m, seq[k:k+1], temp  # test
        return temp
#Union the permutated lsit
#same as the Mathematica Version(Permutations)
def Permutations(l):
    return Union(permutate(l))




def checkdim(inlist, dim):
    dataout=[]
    for i in range(len(inlist)):
        tmp=[]
        for j in range(dim):
            tmp.append(inlist[i][j])
        dataout.append(tmp)
    return dataout



#define FromDigits same as in Mathematica
#input: list of integer
#output: integer for specific base
def FromDigits(listin, b):
    listin= list(map(int,listin))
    lg=len(listin)
    tmp=0
    for i in range(len(listin)):
        tmp=tmp+listin[i]*math.pow(b,lg-(i+1))
    return int(tmp)


def MixedIndexForImproper(c, DIM):
    """

    :param c: cluster or list of items. Example [0 1 0] or [pos1, pos1]
    :param DIM:
    :return: [B matrix to symmetrize the tensor, orbits of the FLATTENED tensor elements]
    """
    npt=len(c)
    t=Tally(c)
    dim= DIM**npt
    # print("len", npt, len(set(c)))
    if npt == len(set(c)):
        return [scipy.sparse.csr_matrix((0, dim)), np.arange(dim)[:,None]]
    allimproper1 = [i[0] for i in t if i[1] >1]
    # print("allimproper1",allimproper1)
    allimproper1= list(set(allimproper1))
    # print("allimproper1",allimproper1)
    allimproper=[[j for j, x in enumerate(c) if x == i] for i in allimproper1]
    allimproperFlat=ListFlat(allimproper)
    # print("all duplicates", allimproper1, allimproper, allimproperFlat)

    includedInOrbit = [0 for _ in range(dim)]
    orbits=[]
    # print("Tuples1([[1,2],[10,20],[300,500]])", Tuples1([[1,2],[10,20],[300,500]]))
    for i in range(dim):
        if includedInOrbit[i] == 1:
            continue

        idx = IntegerDigits(i, DIM, npt)
        #do permutation
        permu=[list(map(list, set(itertools.permutations(idx[j])))) for j in allimproper]
        # print("permu", permu)
        permu0 = [list(itertools.chain.from_iterable(x)) for x in itertools.product(*permu)]
        # print("permu0", permu0)

        # permu0=[]
        # #do Tuple and join
        # for j in Tuples1(permu):
        #     permu0=permu0+j
        # print("idx", idx, permu0, permu, ListFromIndex(idx, allimproper[0]), allimproper[0])

       #permu0 equals permu in the last argument
#        print "After permutation"
#        print permu
#        print "After tuple"
#        print Tuples1(permu)
#        print "After tuple and join"
#        print permu0
        row = idx
        tabletmp=[]
        for p in permu0:
            for k in range(len(p)):
                row[allimproperFlat[k]] = p[k]
#            ListFromIndex(row, allimproperFlat)=p
            irow=FromDigits(row,DIM)
            includedInOrbit[irow]=1
            tabletmp.append(irow)
#        print tabletmp
        orbits.append(tabletmp)

    #construct eqIdx
    eqIdx=[]
    for i in orbits:
        if len(i)>=2:
            tmp=[]
            for j in range(1,len(i)):
                tmp.append([i[0],i[j]])
            eqIdx.append(tmp)
    eqIdx=ListFlat(eqIdx)
#    print "eIdx:"
#    print eqIdx

    #construct B
    B=scipy.sparse.lil_matrix((max([len(eqIdx),1]), dim))
    for i in range(len(eqIdx)):
        B[i, eqIdx[i][0]]=1.0
        B[i, eqIdx[i][1]]=-1.0
    # print(B.todense())
    return [B, orbits]

# #ListArray
# #input: number, dimx, dimy
# #2-dimensional array
# def List2D(n,x,y):
#     # x=int(x)
#     # y=int(y)
#     # tmp=[]
#     # for i in range(x):
#     #     tmp.append([n]*y)
#     # return tmp
#     return np.zeros((x,y))+ n


# #define FCTrans
# #npt,DIM,sgopALL[[iso[[iC,isym,3]],1]], iso[[iC,isym,1]]
# #input: npt; DIM; ith Symmetry Operator(transformed);
# #output: chop needed
# def FCTrans(npt, DIM, gamma, pi):
#     # print("debug pi=", pi, "npt=", npt, "gamma=", gamma)
#     dimTensor = DIM**npt
#     Gm = np.zeros((dimTensor, dimTensor))
#     for i in range(dimTensor):
#         idx1=list(array(IntegerDigits(i,DIM, npt))+1)
#         for j in range(dimTensor):
#             idx2=list(array(IntegerDigits(j,DIM, npt))+1)
#             Gm[i,j]=1.0
#             for k in range(npt):
#                 tmp1=idx1[k]
#                 tmp2=idx2[pi[k]-1]
#                 Gm[i, j]=Gm[i, j]*gamma[tmp1-1][tmp2-1]
#     return scipy.sparse.dok_matrix(Gm)
#



# def getrref(listin):
#     stdr=1E-10
#     use_c_rref = True
#     red= rref(listin) if use_c_rref else cy_rref(listin).todense().tolist()
#     tmp=[]
#     for i in red:
#         # print("debug norm=%g" %( np.linalg.norm(i)))
#         if np.linalg.norm(i)>=stdr: # find nonzero row
#             for ii in range(len(i)): # find first nonozero element
#                 if abs(i[ii])>=stdr:
#                     tmp.append(int(ii))
#                     break
#     return [array(red), tmp]
#
# def nullspace_rref(listin):
#     # print("processing", listin)
#     (xDim, yDim) = listin.shape
#     a=time.time()
# #    reduced=sympy.Matrix(listin).rref()
#     reduced=getrref(listin)
#     # print("reduced", reduced)
#     b=time.time()
#     Bmat=reduced[0].tolist()
#     Bfix=reduced[1]
#     if Bfix==[]:
#         return scipy.sparse.identity(yDim)
#     else:
# #        print "Bfix"
# #        print Bfix
#         Bfre= list(range(yDim))
#
#         for item in Bfix:
#             Bfre.remove(item)
# #        print "Bfre"
# #        print Bfre
#     #construct inmat
#         inmat=[]
#         for i in range(len(Bfix)):
#             tmp=[]
#             for j in Bfix:
#                 tmp.append(Bmat[i][j])
#             inmat.append(tmp)
# #        print "inmat"
# #        print "\n".join(map(str,inmat))
#
#     #construct sol
#         sol=[]
#         for i in range(yDim-len(Bfix)):
#             tmp=[0]*yDim
#             tmp[Bfre[i]]=1.0
#             sol.append(tmp)
#     #print "try-sol"
#     #print "\n".join(map(str,sol))
#
#     #construct y
#         y=[]
#         for i in range(yDim-len(Bfix)):
#             tmp=[]
#             for j in range(len(Bfix)):
#                 tmp.append(-1*Bmat[j][Bfre[i]])
#             y.append(tmp)
#     #print "y"
#     #print "\n".join(map(str,y))
#
#     #solve linear equation
#         for i in range(len(y)):
#             resul=(np.linalg.solve(inmat,y[i])).tolist()
#         #print "resul"
#         #print resul
#             for j in range(len(Bfix)):
#                     sol[i][Bfix[j]]=resul[j]
#     #print "\n"
#         c=time.time()
#         # print("******Time for rref: ", b-a)
#         # print("******Time for lieq: ", c-b)
#         if len(sol)<=0:
#             return scipy.sparse.csr_matrix((0, yDim))
#         else:
#             return scipy.sparse.csr_matrix(sol)


def mychop(arr, tol=1E-10):
    """

    :param arr: numpy array
    :param tol:
    :return:
    """
    if scipy.sparse.issparse(arr):
        nz=arr.nonzero()
        nonzero_mask = np.array(np.abs(arr[nz]) < tol)[0]
        rows = nz[0][nonzero_mask]
        cols = nz[1][nonzero_mask]
        arr[rows, cols] = 0
        arr.eliminate_zeros()
        return arr
    if np.iscomplex(arr).any():
        arr.real[abs(arr.real) < tol] = 0
        arr.imag[abs(arr.real) < tol] = 0
    else:
        arr[abs(arr)<tol]=0
    return arr


def myfrac(arr, tol=1E-10):
    """
     Returns the fractional part (positive or VERY small negative) of array. If fractional
     part > 1-tolerance, return arr-1 (small negative)
    :param arr: array
    :param tol:
    :return:
    """
    b= np.mod(arr, 1)
    b[np.where(b>=1.0-tol)]-=1
    return b


def mkgrid_surfaceintegral_spherical(grid, inv_symm=False):
    thetagrid=np.arccos(np.linspace(0 if inv_symm else -1, 1, grid[0], endpoint=False))
    phigrid=np.linspace(0, 2*np.pi, grid[1], endpoint=False)
    return np.outer(thetagrid, phigrid)


def tensor_constraint(dim, rank, rots, mappings=None, other_constraits=None):
    """
    rank: of tensor
    rots: list of rotation matrix (Cartesian)
    mappings: mapping function (0..N-1 to 0..N-1 one-to-one)
    """
    if rank <= 0:
        null = scipy.sparse.identity(1)
    else:
        dimtensor = dim ** rank
        pi = [list(range(rank))]*len(rots) if mappings is None else mappings
        assert len(pi) == len(rots), ValueError('Number of matrices and mappings %d != %d'%(len(rots), len(pi)))
        Bmats = [fct_trans_c(rank, dim, rots[i], pi[i]) - scipy.sparse.identity(dimtensor) for i in
                 range(len(rots))]
        if other_constraits is not None:
            Bmats.extend([m for m in other_constraits if m.shape[0] > 0])
        if len(Bmats) > 0:
            Bmat = scipy.sparse.vstack(Bmats)
            null = mychop(get_nullspace(Bmat), 1e-12)
        else:
            null = scipy.sparse.identity(dimtensor)
    return null


def get_symmetrized_lsq(mat, bval, alarm_tol=1e-5):
    """
    Least square fittinb of mat.x = bval
    """
    # print('debug using np', np.linalg.lstsq(mat.todense(), bval))
    fit_s = scipy.sparse.linalg.lsqr(mat, bval)[0]
    fit = mat.dot(fit_s)
    err = np.max(fit - bval)
    if err > alarm_tol:
        print("WARNING: large errors in symmetrized least square fitting (max=%8g)" % (err))
    return fit


# root mean square error (RMSE): relative and absolute values
def RMS(l):
    return np.linalg.norm(l)/np.sqrt(len(l)) if len(l)>0 else 1E99


def ndigits(n):
    """
    Number of digits: 1-9:1, 10-99:2, etc
    :param n: integer
    :return:
    """
    return int(np.ceil(np.log10(n+1)))


def groups_of_two(x, removeduplicate=False):
    """
    group a list into groups of 2, no order
    :param x: list of even length
    :return:
    """
    def _group2(l):
        if len(l)<=2:
            return [[sorted(l)]]
        else:
            return [[sorted([l[0], l[i]])] + x for i in range(1,len(l)) for x in _group2(l[1:i]+l[i+1:])]
    x2 = sorted(map(sorted, _group2(x)))
    return list(k for k,_ in itertools.groupby(x2)) if removeduplicate else x2


