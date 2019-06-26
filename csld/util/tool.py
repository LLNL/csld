# to include all module here in order to cite
from numpy import *
from numpy.linalg import *
import string
import os
import scipy
import scipy.sparse
#import rwposcar
#import anaxdat
import math

import numpy as np
import re
import itertools
import os


def list2index_unique(l, uni=None):
    u = list(set(l)) if uni is None else uni
    index_of_uniq = [l.index(i) for i in u]
    return [index_of_uniq[u.index(i)] for i in l]

def arr2str1d(m):
    return ' '.join(map(str, m))

def matrix2text(m_in):
    try:
        m=np.array(m_in)
    except:
        # m is not a matrix, e.g. [[1,2], [3,4,5]]
        return '\n'.join([arr2str1d(x) for x in m_in])
    if len(m.shape) <= 1:
        return arr2str1d(m)
    else:
        return '\n'.join([arr2str1d(x) for x in m])
#    return re.sub(r'\n *', r'\n', re.sub(r'[\[\]]', r'', str(m)).strip(), flags=re.M)

def convert_to_matrix(m):
    if len(m.shape) <= 1:
        return m.reshape((1,-1))
    else:
        return m

def allsubsets(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1))

def my_flatten(list2d):
    return list(itertools.chain.from_iterable(list2d))

def pad_right(a, shape_pad, filling=0):
    """

    :param a: the list to pad
    :param shape_pad: desired shape. Can be integer for 1d array
    :param filling: NOT used yet
    :return:
    """
    # return np.append(listin, [filling]*int(tol-len(listin)))
    shape_p = [shape_pad] if isinstance(shape_pad, int) else shape_pad
    return np.lib.pad(a, tuple((0, max(0,shape_p[i]-a.shape[i])) for i in range(a.ndim)), 'constant')


def non_1to1(pi):
    """

    :param pi: a mapping/indexing array
    :return: non unique index in range(len(pi)), non unique pi elements
    """
    counts = np.bincount(pi, minlength=len(pi))
    idxj = np.where(counts!=1)[0].tolist()
    return np.where([pi[i] in idxj for i in range(len(pi))])[0].tolist(), idxj


# these two functions are taken from the monty module
from functools import wraps


def singleton(cls):
    """
    This decorator can be used to create a singleton out of a class.

    Usage::

        @singleton
        class MySingleton():

            def __init__():
                pass
    """

    instances = {}

    def getinstance():
        if cls not in instances:
            instances[cls] = cls()
        return instances[cls]
    return getinstance


def cached_class(klass):
    """
    Decorator to cache class instances by constructor arguments.
    This results in a class that behaves like a singleton for each
    set of constructor arguments, ensuring efficiency.

    Note that this should be used for *immutable classes only*.  Having
    a cached mutable class makes very little sense.  For efficiency,
    avoid using this decorator for situations where there are many
    constructor arguments permutations.

    The keywords argument dictionary is converted to a tuple because
    dicts are mutable; keywords themselves are strings and
    so are always hashable, but if any arguments (keyword
    or positional) are non-hashable, that set of arguments
    is not cached.

    Example::

        >>> @cached_class
            class A(object):
                def __init__(self, val):
                    self.val = val
            ...
        >>> a1a = A(1)
        >>> a1b = A(1)
        >>> a2 = A(2)
        >>> id(a1a) == id(a1b)
        True
        >>> id(a1a) == id(2)
        False
    """
    cache = {}

    @wraps(klass, assigned=("__name__", "__module__"), updated=())
    class _decorated(klass):
        # The wraps decorator can't do this because __doc__
        # isn't writable once the class is created
        __doc__ = klass.__doc__

        def __new__(cls, *args, **kwargs):
            key = (cls,) + args + tuple(kwargs.items())
            try:
                inst = cache.get(key, None)
            except TypeError:
                # Can't cache this set of arguments
                inst = key = None
            if inst is None:
                # Technically this is cheating, but it works,
                # and takes care of initializing the instance
                # (so we can override __init__ below safely);
                # calling up to klass.__new__ would be the
                # "official" way to create the instance, but
                # that raises DeprecationWarning if there are
                # args or kwargs and klass does not override
                # __new__ (which most classes don't), because
                # object.__new__ takes no parameters (and in
                # Python 3 the warning will become an error)
                inst = klass(*args, **kwargs)
                # This makes isinstance and issubclass work
                # properly
                inst.__class__ = cls
                if key is not None:
                    cache[key] = inst
            return inst

        def __init__(self, *args, **kwargs):
            # This will be called every time __new__ is
            # called, so we skip initializing here and do
            # it only when the instance is created above
            pass

    return _decorated


#define touch file
def touch(file):#input string
    if os.path.isfile(file):
        os.system(str("rm"+" "+file))
        os.system(str("touch"+" "+file))
    else:
        os.system(str("touch"+" "+file))

def mkdir(dir):
    if os.path.isdir(dir):
        os.system(str("rm"+" -r "+dir))
        os.system(str("mkdir"+" "+dir))
    else:
        os.system(str("mkdir"+" "+dir))
if False:
    mkdir("xixi/")
#define rm file
def rm(file):
    if os.path.isfile(file):
        os.system(str("rm"+" "+file))
    else:
        print("No file found, dont need to rm")
    
#define check file(1 exist; else0)
def check(file):
    if os.path.isfile(file):
        return int(1)
    else:
        return int(0)

#define check the file status (print the status)
def checkfile(file):
    if os.path.isfile(file):
        print(str(file)+" exists :)")
    else:
        print(str(file)+" not found :(")

#define readallline function
def readinline(file):
    dataout=[]
    if check(file):
        fin=open(file,"r")
        for line in fin:
            dataout.append(line.split())#map(float,line.split()))
        fin.close()
    else:
        print(str(file)+" not found :(")
    return array(dataout)

#define write1dmat
def write1dmat(datain, file):
    if check(file):
        rm(file)
        touch(file)
    else:
        touch(file)
    fout=open(file, "w")
    fout.writelines("\n".join(map(str,datain)))
    fout.close()

#define write2dmat
def write2dmat(datain, file):
    if check(file):
        rm(file)
        touch(file)
    else:
        touch(file)
    fout=open(file, "w")
    #cout line number
    fout.writelines(str(len(datain))+"\n")
    for i in datain:
        fout.writelines(" ".join(map(str,i))+"\n")
    fout.close()

#define write2dMTX
def write2dMTX(datain, file):
    if check(file):
        rm(file)
        touch(file)

    else:
        touch(file)
    fout=open(file, "w")
    fout.writelines("%%MatrixMarket matrix coordinate real general\n")
    fout.writelines("%Created by Wolfram Mathematica 9.0 : www.wolfram.com\n")
    print("Transfering to sparse matrix----")
    BB=scipy.sparse.coo_matrix(datain)
    print("Spare matrix obtained!")
#    print BB.row
#    print BB.col
#    print BB.data
    fout.writelines(str(len(datain))+" "+str(len(datain[0]))+" "+str(len(BB.data))+"\n")
    for i in range(len(BB.data)):
        fout.writelines(str(BB.row[i]+1)+" "+str(BB.col[i]+1)+" "+str(BB.data[i])+"\n")
    #for i in range(len(datain)):
        #for j in range(len(datain[0])):
            #fout.writelines(str(i+1)+" "+str(j+1)+" "+str(datain[i][j])+"\n")
    fout.close()
    
def read2dMTX(file):
    if check(file):
        counter=0
        for line in open(file):
            counter=counter+1
            if counter <=2:
                continue
            if counter ==3:
                inlist=list(map(int,line.split()))
                nrow=inlist[0]
                ncol=inlist[1]
                dataout=array([[0.0]*ncol]*nrow)
                continue
            if counter >=4:
                tmp=line.split()
                #print str(tmp)+", "+str(tmp[2])
                dataout[int(tmp[0])-1][int(tmp[1])-1]=float(tmp[2])
                #print "\n"
        return dataout.tolist()
    else:
        print(str(file)+" not found :(")

#test
if False:
    Amat=[[0,1],[2,0],[0,0],[0,16]]
    print(Amat)
    write2dMTX(Amat, "test.mtx")
    print(read2dMTX("test.mtx"))

#define read1dmat
#read float
def read1dmat(file):
    mat=[]
    if check(file):
        for line in open(file):
            mat.append(float(line))
        return mat
    else:
        print(str(file)+" not found :(")

if False:
    haha=[1,2,3,4,5]
    write1dmat(haha, "haha")
    xixi=read1dmat("haha")
    print(xixi)
#define read2dmat (this is a relatively fast way: iter or chunck read)
def read2dmat(file,icomplex=False):
    mat=[]
    if check(file):
        print("Read matrix start")
        for line in open(file):
            if not icomplex:
                mat.append(list(map(float,line.split())))
            else:
                mat.append(list(map(complex,line.split())))
        print("Read matrix end")
        #delete line counter
        del mat[0]
        return mat    
    else:
        print(str(file)+" not found :(")


#test
#mat=read2dmat("C-isoo.mat")
#print len(mat)
#print len(mat[0])

def clusstr(clus):
    dataout=""
    for item in clus:
        dataout=dataout+str(item[0])+" "+str(item[1])+" "+str(item[2])+"\n"
    return dataout

def lptstr(lpt):
    dataout=""
    for item in lpt:
        dataout=dataout+str(item[0][0])+" "+str(item[0][1])+" "+str(item[0][2])+" "+str(item[1])+"\n"
    return dataout


#define writeorb(orb)
def writeorb(orbset, file):
    if check(file):
        rm(file)
        touch(file)
    else:
        touch(file)
    fout=open(file, "w")
    fout.write(str(len(orbset))+"\n\n")
    for orb in orbset:
        fout.write(str(len(orb))+"\n\n")
        for item in orb:
            npt=len(item[0])
            fout.write(str(npt)+"\n")
            fout.write(clusstr(item[0]))
            fout.write(str(item[1])+"\n")
            fout.write(str(item[2])+"\n")
            fout.write(lptstr(item[3]))
        fout.write("\n")
    fout.close()

def writeclus(clus, file):
    if check(file):
        rm(file)
        touch(file)
    else:
        touch(file)
    fout=open(file,"w")
    fout.write(str(len(clus))+"\n\n")
    for item in clus:
        fout.write(str(len(item))+"\n")
        fout.write(clusstr(item))
        fout.write("\n")
    fout.close()

def writeSCinfo(SCinfo, file):
    if check(file):
        rm(file)
        touch(file)
    else:
        touch(file)
    fout=open(file, "w")
    tmp=[SCinfo['SC'], SCinfo['invSC'], SCinfo['SCref'], SCinfo['SCpos'], SCinfo['SCmat'], SCinfo['invSCmat'], SCinfo['order']]
    lentmp=[len(i) for i in tmp]
    fout.write(" ".join(map(str,lentmp))+"\n")
    for i in tmp:
        if i==SCinfo['order']:
            fout.write("\n".join(map(str,i))+"\n")
        else:
            for j in i:
                fout.write(" ".join(map(str,j))+"\n")
    fout.close()

def readSCinfo(file):
    SCinfo={}
    if check(file):
        fin=open(file, "r")
        lenlist=list(map(int,(fin.readline()).split()))
#        tmp=[SCinfo['SC'], SCinfo['invSC'], SCinfo['SCref'], SCinfo['SCpos'], SCinfo['SCmat'], SCinfo['invSCmat'], SCinfo['order']]
        tmp=[]
        for i in range(7):
            tmp1=[]
            for j in range(lenlist[i]):
                if i in [0,1,3,4,5]:
                    tmp1.append(list(map(float,(fin.readline()).split())))
                elif i in [2]:
                    tmp1.append(list(map(int,(fin.readline()).split())))
                else:
                    tmp1.append(map(int,(fin.readline()).split())[0])
            tmp.append(tmp1)
        SCinfo['SC']=tmp[0]
        SCinfo['invSC']=tmp[1]
        SCinfo['SCref']=tmp[2]
        SCinfo['SCpos']=tmp[3]
        SCinfo['SCmat']=tmp[4]
        SCinfo['invSCmat']=tmp[5]
        SCinfo['order']=tmp[6]
    else:
        print(str(file)+" not found :(")
    return SCinfo

    
def readclus(file):
    if check(file):
        fin=open(file, "r")
        nclus=int(fin.readline())
        clus=[]
        for i in range(nclus):
            item=[]
            fin.readline()
            npt=int(fin.readline())
            for j in range(npt):
                item.append(list(map(float, fin.readline().split())))
            clus.append(item)
        return clus
    else:
        print(str(file)+" not found :(")
#writeclus(clus,"uniqueC")
#print "\n".join(map(str, readclus("uniqueC")))

def readorb(file):
    if check(file):
        orbset=[]
        fin=open(file, "r")
        Norb=int(fin.readline())
        for i in range(Norb):
            orb=[]
            fin.readline()
            nitem=int(fin.readline())
            fin.readline()
            for j in range(nitem):
                item=[]
                npt=int(fin.readline())
                clus=[]
                lpt=[]
                for k in range(npt):
                    line=fin.readline()
                    clus.append(list(map(float,line.split())))
                item.append(clus)
                item.append(int(fin.readline()))
                item.append(int(fin.readline()))
                for k in range(npt):
                    line=fin.readline()
                    tmp=list(map(float,line.split()))
                    tmp=list(map(int, tmp))
                    lpt.append([[tmp[0],tmp[1],tmp[2]],tmp[3]])
                item.append(lpt)
                orb.append(item)
            orbset.append(orb)
        fin.close()
        return orbset
    else:
        print(str(file)+" not found :(")


#def read fit.ou
def readfit(file):
    if check(file):
        counter=0
        readflag=False
        for line in open(file):
            counter=counter+1
            if counter==1:
                nstruc=map(int, line.split())[1]
                fitlist=[0.0]*nstruc
            if len(line.split())>=1 and (line.split())[0]=="found":
                readflag=True
                continue
            if readflag:
                index=int((line.split())[0])
                resl=float((line.split())[1])
                fitlist[index-1]=resl
        print("Fit.our read successfully and length: "+str(len(fitlist)))
        return fitlist
    else:
        print(str(file)+" not found :(")
#test:
if False:
    print(readfit("fit.out-mu1"))
