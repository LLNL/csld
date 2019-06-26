#!/usr/bin/env python3
"""
test_bregman.py - Mimic the behavior of test_bregman.f90 but in Python.

To compile the Breagman library, do the following::

    ifort -nus -fpic -cpp -fast -shared -o libbregman.so bregman.f90

"""

import scipy.io
import numpy as np
from scipy.sparse.linalg.eigen.arpack import eigsh
try:
    from bregman import bregman
except ImportError:
    print("Failed to import Bregman subroutines")
    pass


def bregman_func(A_in, En_in, method=5, mu=0.001, lbd=3., maxIter=100, tol=1E-4):
    # Open the A matrix file
    if isinstance(A_in, str):
        A = scipy.io.mmread(A_in).astype(np.double)
    else:
        A = A_in
    Nstruc, Ncorr = A.shape
    if isinstance(En_in, str):
        En = np.loadtxt(En_in).astype(np.double)
    else:
        En = En_in
    # print(A, En, En.shape)
    assert Nstruc == En.shape[0]
    ecis = np.zeros(Ncorr)

    if method <=2:
        # FPC requires that max eigenvector of A.A^T <=1
        rescale_factor = eigsh(A.dot(A.T), 3, which='LM', return_eigenvectors=False)[-1]
        rescale_factor = np.sqrt(rescale_factor)
        A = A/rescale_factor
        En = En/rescale_factor
    tau = min(1.999, -1.665*float(Nstruc)/float(Ncorr) + 2.665)
        # print("FPC step size tau =", tau)

    # print("method=", method)
    # print("mu=", mu)
    # print("lambda=", lbd)

    # Values in Fortan are passed by reference, so we need to convert the types
    # and then pass by reference
    if (method in (1,3,5)) and scipy.sparse.issparse(A):
            A = A.todense()

    # Now scale the mu parameter by Nstruc for convenience
    mu2 = mu*Nstruc

    if method == 1:
        ecis= bregman.bregmanfpc(maxIter, tol, tau, mu, A, En)
    elif method == 2:
        ecis= bregman.sparsebregmanfpc(maxIter, tol, tau, mu, A, En)
    elif method ==3:
        ecis= bregman.splitbregman(maxIter, tol, tau, mu, A, En)
    elif method == 4:
        ecis= bregman.sparsesplitbregman(maxIter, tol, tau, 1/mu, A, En)
    elif method == 5:
        # right preconditioning
        ecis= bregman.bregmanrprecond(maxIter, A, En, 1/mu, lbd, tol)
    elif method == 6:
        # sparse right preconditioning
        print(" WARNING: DO NOT USE; NOT IMPLEMENTED YET")
        ecis= bregman.sparsebregmanrprecond(A, En, 1/mu, lbd)
    else:
        raise ValueError("Unknown bregman method %d", method)

    return ecis

if __name__ == "__main__":
    # TODO: Support command line arguments
    print([bregman_func('A.mtx', 'En.dat', mu=mu, lbd=10) for mu in (1E-7, 1E-6, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1., 10., 300.)])


