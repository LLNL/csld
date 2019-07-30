#!/usr/bin/env python3

"""
Fit a linear model
    A x = b
"""


import numpy as np
from numpy.linalg import norm
import scipy as sp
import scipy.sparse
try:
    from cssolve.bregman_func import bregman_func
except ImportError:
    print("Failed to import bregman_func")
    pass
try:
    from bcs_driver import bcs_driver
except ImportError:
    #print("Failed to import bcs solver")
    pass
try:
    from sklearn import linear_model
except ImportError:
    print("Failed to import sklearn.linear_model")
    pass
try:
    import matplotlib
    matplotlib.use('AGG')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
except ImportError:
    print("WARNING: cannot import pyplot; plotting disabled")
    plt = None
    pass

import logging
logger = logging.getLogger(__name__)

debug_level = 10

# root mean square error (RMSE): relative and absolute values
def RMS(l):
    return norm(l)/np.sqrt(len(l)) if len(l)>0 else 1E99

# leave-one-out cross validation
# A.x=b
def LOOCV(A, x, b):
    hmat= A.dot(np.linalg.inv(A.T.dot(A))).dot(A.T)
    return RMS((b-A.dot(x))/(1-np.diag(hmat)))

def get_errors(l0, l1):
    if len(l0)<=0:
        return [100, 1E99]
    return [100*norm(l0-l1)/norm(l0) if norm(l0)>1E-8 else float('inf'), RMS(l0-l1)]

def L01norm(ls, tol=1E-16):
    sum0=0
    sum1=0.0
    for x in ls:
        sum1 += abs(x)
        if abs(x)> tol:
            sum0 += 1
    return [sum0, sum1]

def L1norm(ls):
    return L01norm(ls)[1]

def predict_holdout(A, fval, sol, errbar=[]):
    f_pred = A.dot(sol)
    err = get_errors(fval, f_pred)
    return err + [fval, f_pred, sol, errbar]


def analyze_solutions(sols, threshold=0.9):
    """

    :param sols: list of solution vectors
    :param threshold: set to zero if abs(<x>)< threshold* std(x)
    :return: averaged solutions
    """
    # remove possible nan
    dim = len(sols[0])
    sols= sols[np.where(np.any(np.isnan(sols), axis=1)==False)]
    # cannot analyze if there is only one solution!
    if len(sols) ==1:
        return (sols[0], np.zeros_like(sols[0]))
    elif len(sols)==0:
        return (np.zeros(dim), np.zeros(dim))

    sol= np.average(sols, axis=0)
    err= np.std(sols, axis=0)
    idx = np.where(abs(sol) < threshold*err)
    sol[idx] = 0.0
    err[idx] = 0.0
    return (sol, err)


def tikhonov_reg(A, b, mu):
    """
    Tikhonov Regularization or ridge regression, |Ax-b|^2 + mu x^2
    :param A:
    :param b:
    :param mu:
    :return:
    """
    return sp.linalg.inv(A.T.dot(A) + mu* np.eye(A.shape[1])).dot(A.T.dot(b))



def fit_in_steps(ld, steps, sol0, *args, **kwargs):
    """
    fit in steps, each steps =
       [description, input_orders, submodels]
    """
    for i,step in enumerate(steps):
        print("Step %d %s"%(i, step[0]))
        knownsol= ld.select(sol0 if i==0 else solutions[ibest], step[1])
        ibest, solutions = csfit(args, kwargs, submodels=step[2],knownsol=knownsol)



def csfit(Amat, flist_in, wt, mulist, method=5, submodels=None, nSubset=1, subsetsize=0.8,
          trunc_threshold=0.9,
          reweight=True, penalty='arctan', jcutoff=1E-3,sigma2=-1.0,eta=1E-8,
          holdsize=0.06, lbd=3.0, pdfout=None, maxIter=400, tol=1E-4, fitf=None):
    """
    all data= holdout +training
    training = cross validation + the rest or subset, i.e. actual training
    holdout set will be excluded from any fitting, good if there are plenty of data points
    cross-validation error will be performed, for each subset fitting, on all-holdout-subset data
    Choice of best solution is based on holdout error if present, or cv error otherwise
    :param Amat:
    :param flist_in:
    :param wt:
    :param mulist:
    :param method:
    :param submodels:
    :param nSubset:
    :param subsetsize:
    :param trunc_threshold:
    :param reweight:
    :param penalty:
    :param jcutoff:
    :param holdsize:
    :param lbd:
    :param pdfout:
    :param maxIter:
    :param tol:
    :param fitf:
    :return:
    """
    # print("Amat", Amat.__class__)
    nVal, nCorr = Amat.shape
    wtlist = np.ones(nVal) if wt==1 else wt
    holdsize= int(round(nVal*holdsize)) if isinstance(holdsize, float) else holdsize
    holdsize= max(0,min(nVal, holdsize))
    holdout= sorted(np.random.choice(nVal, holdsize, replace=False))
    subsetsize= int(round(nVal*subsetsize)) if isinstance(subsetsize, float) else min(nVal-holdsize, subsetsize)
    alltrain = list(set(range(nVal))-set(holdout))
    cvsize= nVal-holdsize-subsetsize
    print("  No. of data points: all=%d  holdout=%d  cross-v=%d  train=%d"%(nVal, holdsize, cvsize, subsetsize))

    #print(submodels)
    if submodels is None:
        submodels = [['All', scipy.sparse.identity(nCorr), np.zeros(nCorr)]]
    if nSubset<=1:
        print("WARNING: at least two subsets required to estimate fitting parameter variance")
    #if cvsize<=0:
    #    print("cross validation set empty, skipping subsets")
    #    nSubset=1

    sol_all = []
    rms_all = []
    inRMS= RMS(flist_in)
    # use theC to select only some parameters for fitting, give weight to parameters, etc
    for model_name, theC, sol0_in in submodels:
        print("  considering", model_name, " nvar=", theC.shape[1])
        sol0 = np.ones(nCorr) * sol0_in  # in case sol0_in == 0
        flist0 = Amat.dot(sol0)
        flist = flist_in - flist0
        modelErr= get_errors(flist_in, flist0)
        modelRMS = modelErr[1]
        err_min_model = 1.E99
        if np.linalg.norm(sol0)>1E-16:
            print("  Error after subtracting known solution:", modelErr)
        AC = Amat.dot(theC)
        # each "fit" is: [rms_err, fval, f_pred, sol, errbar]
        fit_model= []
        if method==101:
            # BCS, so mu is not needed
            mulist=[0]
        for mu in mulist:
            # Bayesian CS
            if method == 101:
               # sols,error_bars= bcs_driver.f90wrap_do_wrapped(AC.todense()[alltrain], flist[alltrain], 1E-12, 1E-18)
               # sols=sols.reshape((1,-1))
               # cv_rms= [get_errors(flist[holdout], AC[holdout].dot(sols[0]))]
               # print('debug sols, errr', sols.shape,error_bars.shape)
               # print('debug reweight,penalty,jcutoff,sigma2,eta', reweight,penalty,jcutoff,sigma2,eta)
                fit_rms, fit_abserr, cv_rms, cv_abserr, error_bars, sols = bcs_driver.bcs_fit(cvsize, reweight,
                               penalty, jcutoff, np.ones(nSubset)*sigma2, eta, AC.todense()[alltrain], flist[alltrain])
  #              rmse = np.mean(hold_rms)
  #              holdout_rel_abs_err = [100*rmse/(norm(flist)/np.sqrt(len(flist))), rmse]
  #               print("fit_rms, fit_err, cv_rms, cv_err", fit_rms, fit_abserr, cv_rms, cv_abserr)
  #               print("fit_rms, fit_err, cv_rms, cv_err", fit_rms.shape, fit_abserr.shape, cv_rms.shape, cv_abserr.shape)
  #               print('error_bars, sols=', error_bars, sols)
  #               print('error_bars, sols=', error_bars.shape, sols.shape)
                cv_errs = cv_rms
                # we don't know the random sample in bcs, so just make one up
                cv1=sorted(np.random.choice(alltrain, cvsize, replace=False))
            else:
                sols = []
                cv_errs = []
                for itrain in range(nSubset):
                    ts1= sorted(np.random.choice(alltrain, subsetsize, replace=False))
                    cv1=[x for x in alltrain if x not in ts1]
                    if method == 201:
                        # Tikhonov Regularization or ridge regression, |Ax-b|^2 + mu x^2
                        sol= tikhonov_reg(AC[ts1], flist[ts1], mu*subsetsize)
                    else:
                        sol = bregman_func(AC[ts1], flist[ts1], method, mu, lbd=lbd,
                                           maxIter=maxIter, tol=tol)
                    sols.append(sol)
                    cv_errs.append(get_errors(flist[cv1], AC[cv1].dot(sol))[1])
                #errs.extend([get_errors(flist[holdout], AC[holdout].dot(sol)) for sol in sols])
                #holdout_rel_abs_err = np.mean(errs, axis=0)
                sols = np.array(sols)

            sol, errbar = analyze_solutions(sols, trunc_threshold)
            l01 = L01norm(sol)
            fit_hold = predict_holdout(AC[holdout], flist[holdout], sol, errbar)[1:]
            fit_cv = predict_holdout(AC[cv1], flist[cv1], sol, errbar)[1:]
            fit_cv[0] = RMS(cv_errs)
            #err[:2] = holdout_rel_abs_err
            fit_mu = fit_hold if holdsize>0 else fit_cv
            # if fit_mu[0] < err_min:
            #     err_min = fit_mu[0]
            #     ibest = isol

            # sol = sol[1]
            # tmp= [mu, sol[0], l01[0], l01[1], np.mean([x[0] for x in fitmu]), np.min([x[0] for x in fitmu]), sol]
            fit_mu=[mu, fit_mu[0]/inRMS*100, fit_mu[0], fit_cv[0], l01[0], l01[1]]+fit_mu
            fit_model.append(fit_mu)
            print("    mu= %.3g rel_err= %5.2f %%  err= %9.5f err_cv= %9.5f L0= %d L1= %8.3f"%tuple(fit_mu[:6]))
            np.savetxt(open('solution_ALL','ab'), (theC.dot(sol)+sol0)[None,:])

        bestsol = min(fit_model, key=lambda x: x[2])
        # print('debug bestsol=', bestsol)
        i_mu= fit_model.index(bestsol)
        sol_all.append(theC.dot(bestsol[-2]) + sol0)
        rms_all.append(bestsol[2])
#        np.savetxt('solution_all_'+model_name, [theC.dot(x[-2])+sol0 for x in fit_model])
        print("     mu    err%%      err    err_cv     L0   L1(scaled)     L1")
        for i, s in enumerate(fit_model):
            print("%8.2g %7.2f  %7.4g  %7.4g  %6d   %7g   %7g" % tuple(s[:6]+[L1norm(theC.dot(s[-2]))]), '  \033[92m<==BEST\033[0m' if i==i_mu else '')

        if plt is not None:
            def writepdf():
                if pdfout is not None:
                    plt.savefig(pdfout, format='pdf')
                    plt.close()

            if len(mulist)>1:
                plt.figure()
                plt.plot([x[0] for x in fit_model], [x[1] for x in fit_model], 'bo-')
                plt.xscale('log')
                plt.yscale('log')
                plt.xlabel("mu")
                plt.ylabel("rel. Err")
                writepdf()

            plt.figure().add_subplot(111, aspect=1)
            plt.plot(bestsol[7], bestsol[8], 'bo')
            plt.xlabel('Input')
            plt.ylabel('Prediction')
            plt.title('RMSE= %6f  (%.2f%%)'% (bestsol[2], bestsol[1]))
            writepdf()
            if fitf is not None:
                np.savetxt(fitf, np.array([flist_in, Amat.dot(sol_all[-1])]).T)
                # with open(fitf, 'ab') as f:
                #     np.savetxt(f, np.array(bestsol[7:9]).T)

            plt.figure()
            plt.vlines(np.arange(1,bestsol[9].shape[0]+1), 0, bestsol[9], color='k', linestyles='solid')
            plt.axes().errorbar(np.arange(1,bestsol[9].shape[0]+1), bestsol[9], yerr=bestsol[-1], fmt='o')
            plt.ylabel("fitted values")
            plt.axhline(0, color='black', lw=1)
            writepdf()

            # plt.show()

        # fittings.append(bestsol[-1])

    print(rms_all)
    rel_err = fit_mu[1]
    return np.argmin(rms_all), np.array(sol_all), rel_err
