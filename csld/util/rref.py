"""
  Reduced Row Echelon - Python interface to C implementation.
  
  To compile the C code do:
  
    gcc -fPIC -shared rref.c -o _rref.so
  
  Usage:
  
    from rref import rref
    rref([[1,2,-1,-4], [2,3,-1,-11], [-2,0,-3,22]])
    
  Output:
    [[1, 0, 0, -8], [0, 1, 0, 1], [0, 0, 1, -2]]
"""
import ctypes
#import scipy.sparse
import numpy as np

try:
    lib = ctypes.cdll.LoadLibrary("_rref.so")
except OSError:
    lib = ctypes.cdll.LoadLibrary("./_rref.so")

class Matrix(ctypes.Structure):
    _fields_ = [("dim_x",   ctypes.c_int32),
                ("dim_y",   ctypes.c_int32),
                ("m_stor",  ctypes.POINTER(ctypes.c_double)),
                ("mtx",     ctypes.POINTER(ctypes.POINTER(ctypes.c_double)))]


lib.MtxToReducedREForm2.restype  = ctypes.POINTER(Matrix)
lib.MtxToReducedREForm2.argtypes = [ctypes.c_int32, ctypes.c_int32,
                                    ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]

lib.NewMatrix.argtypes = [ctypes.c_int32, ctypes.c_int32]
lib.NewMatrix.restype  = ctypes.POINTER(Matrix)

def rref(matrix):
    # Convert from python list to c arrays
    im = (ctypes.POINTER(ctypes.c_double)*len(matrix))()
    for i, row in enumerate(matrix):
        im[i] = (ctypes.c_double * len(row))(*row)

    # Do RREF
    m = lib.MtxToReducedREForm2(len(matrix[0]), len(matrix), im)

    # Convert back to python list
    om = []
    for i in range(m.contents.dim_y):
        row = []
        for j in range(m.contents.dim_x):
            row.append(round(m.contents.mtx[i][j],5))
        om.append(row)

    # Clean up
    lib.FreeMatrix(m)
    del im
    del m

    # And we're done
    return np.array(om)

def test():
    import random

    im = (ctypes.POINTER(ctypes.c_int)*81)()
    for i in range(81):
        r = (ctypes.c_int*81)()
        for j in range(81):
            r[j] = random.randint(1,10)
        im[i] = r

    lib.MtxToReducedREForm2(81,81,im)

if __name__ == "__main__":
    test()
