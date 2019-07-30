#!/usr/bin/env python3

"""
Misc subroutines for linear fitting:
 * init correlation matrix A, value vector b
 * fit to get unknown coefficients x
 * predict using fitted x
"""

import re
import glob
import logging
import numpy as np
import scipy
import scipy.sparse

from cssolve.csfit import csfit, predict_holdout


def print_version(version):
    print(" " * 3, version, "\n")


def print_end():
    print("Done")


def add_common_parameter(parser):
    """
    Add a few command-line parameters common to different models
    :return:
    """
    parser.add_argument("--log_level", type=int, help="Logging level. Default 1", default=1)
    parser.add_argument("--symm_step", type=int, help="Space group. 1 = file, *2 = spglib, 3 = spglib & save", default=2)
    parser.add_argument("--symm_prim", action='store_false',
                        help="Symmetrize primitive cell with space group and save to POSCAR_sym. Default: True", default=True)
    parser.add_argument("--clus_step", type=int, help="Clusters. 0 = ignore & exit, 1 = file, 2 = generate, *3 = generate & save", default=3)
    parser.add_argument("--symC_step", type=int, help="Independent parameters. 0 = ignore & exit, 1 = file, 2 = generate, *3 = generate & save", default=3)
    parser.add_argument("--train_step", type=int, help="Correlation matrix. 0 = ignore & exit, 1 = file, 2 = generate, *3 = generate & save, 4 = skip", default=3)
    parser.add_argument("--fit_step", type=int, help="Fitting. 0 = ignore & exit, 1 = file, 2 = generate, *3 = generate & save", default=3)
    parser.add_argument("--pred_step", type=int, help="Prediction. *0 = skip, 1 = file, 2 = generate, 3 = generate & save", default=0)
    parser.add_argument("--refit", action='store_true',
                        help="Perform another fitting, equivalent to \"--clus_step 1 --symC_step 1 --train_step 1\". Default: False", default=False)
    parser.add_argument("--cont", "-c", action='store_true',
                        help="Continue from previous run, equivalent to \"--clus_step 1 --symC_step 1 --train_step 4 --fit_step 1\". Default: False", default=False)
    parser.add_argument("--predict", action='store_true',
                        help="predict, equivalent to \"--cont --pred_step 2 and skipping magnon/phonon steps\". Default: False", default=False)
    parser.add_argument("--override", "-o", action="append", help="Override setting file from command line, e.g. \"[structure] epsilon_inf=eps.txt\"")


def process_common_options(options):
    # Quick prediction mode
    if options.predict:
        options.cont= True
        options.pred_step = 2
    # Assuming fitting is available from previous run
    if options.cont:
        options.clus_step = 1
        options.symC_step = 1
        options.train_step = 4
        options.fit_step = 1
    # Another fit reusing existing info
    if options.refit:
        options.clus_step = 1
        options.symC_step = 1
        options.train_step = 1


def override_from_commandline(override_list, settings):
    if not override_list:
        return
    pattern = re.compile(r"\[\w+\] *\w+ *=.*")
    for x in override_list:
        if not pattern.match(x):
            print('ERROR: %s is not a valid override. Expecting e.g. [sec] name=val'%(x))
            return
        sec, other = re.split(r'\]', x[1:], maxsplit=1)
        sec = sec.lower()
        other = other.strip()
        tag, val = re.split(r'=', other, 1)
        tag = tag.rstrip()
        val = val.lstrip()
        if sec not in settings:
            settings.add_section(sec)
        settings.set(sec, tag, val)


def upon_exit(pdfout):
    if pdfout is not None:
        pdfout.close()


def init_training(model, setting, step, **kwargs):
    """
    training structures from which to obtain the Correlation matrix
    :return:
    """
    from csld.util.io_utils import co, load_matrix
    if step <= 0:
        exit(0)
    elif step == 4:
        Amat = None
        fval = None
    elif step == 1:
        Amat = scipy.sparse.csr_matrix(load_matrix(setting['corr_in']))
        fval = np.zeros((Amat.shape[0], 3))
        fval[:, 0] = np.loadtxt(setting['fval_in'])
    elif step in [2, 3]:
        traindat= [y.split() for x, y in setting.items() if re.match('traindat.*', x) is not None]
        Amat, fval = model.get_correlation([[sc[0],
                                        [f for subs in sc[1:]
                                        for f in sorted(glob.glob(subs))]] for sc in traindat],
                                           corrtype=setting['corr_type'], setting=setting, **kwargs)
        if step == 3:
            scipy.io.mmwrite(setting['corr_out'], Amat)
            np.savetxt(setting['fval_out'], fval[:, 0])
    else:
        print("ERROR: Unknown corr_step: ", step)
        exit(-1)
    print("+ Correlation matrix for training done")
    return Amat, fval


def fit_data(model, Amat, fval, setting, step, pdfout):
    """
    Fitting
    :param model
    :param Amat:
    :param fval:
    :param setting:
    :param step:
    :param pdfout:
    :return: optimal solution
    """
    if step <= 0:
        exit(0)
    elif step == 1:
        solutions = model.load_solution(setting['solution_in'])
        if Amat is not None:
            err = [np.std(Amat.dot(solutions[i])-fval[:,0]) for i in range(solutions.shape[0])]
            ibest = np.argmin(err)
        else:
            ibest = 0
            if solutions.size <= 0:
                logging.error("ERROR: empty solution")
                exit(-1)
            if solutions.shape[0] > 1:
                logging.warning("More than 1 solutions found. Returning the first.")
    elif step in [2, 3]:
        mulist = list(map(float, setting['mulist'].split()))
        submodels = [y.split() for x, y in setting.items() if re.match('submodel.*', x) is not None]
        submodels = [[x[0], ' '.join(x[1:])] for x in submodels]
        knownsol = setting.get('solution_known', '')
        submodels = model.get_submodels(submodels, setting=setting, knownsol=knownsol)

        ibest, solutions, rel_err = csfit(Amat, fval[:,0], 1, mulist,
                method=int(setting['method']),
                maxIter=int(setting['maxiter']),
                tol=float(setting['tolerance']),
                nSubset=int(setting['nsubset']),
                subsetsize=float(setting['subsetsize']),
                holdsize=float(setting['holdsize']),
                lbd=float(setting['lambda']),
# bcs options
                reweight=setting.getboolean('bcs_reweight', False),
                penalty=setting.get('bcs_penalty', 'arctan'),
                jcutoff=setting.getfloat('bcs_jcutoff',1E-7),
                sigma2=setting.getfloat('bcs_sigma2',-1.0),
                eta=setting.getfloat('bcs_eta',1E-3),

                fitf=setting.get('true_v_fit'),
                submodels=submodels, pdfout=pdfout)
        if step == 3:
            np.savetxt(setting['solution_out'], solutions)
            np.savetxt(setting['solution_out']+'_full', model.Cmat.T.dot(np.array(solutions)[:,:model.Cmat.shape[0]].T).T)

    else:
        print("ERROR: Unknown fit_step: ", step)
        exit(-1)
    print("+ Fitting done. Best solution", ibest)
    return ibest, solutions, rel_err


def predict(model, sols, setting, step):
    """

    :param model:
    :param sol:
    :param setting:
    :param step:
    :return:
    """
    if step <= 0:
        return
    elif step in [1, 2, 3]:
        Amat, fval = init_training(model, setting, step, delForce=0)
    else:
        print("ERROR: Unknown pred_step: ", step)
        exit(-1)

    errs = []
    for i in range(len(sols)):
        err = predict_holdout(Amat, fval[:, 0], sols[i])
        errs.append(err[0])
        print("  sol# %d: err= (%.2f%%) %f" % (i, err[0], err[1]))
        np.savetxt("%s_%d"%(setting['fval_out'],i), np.transpose(err[2:4]))
    print("+ Prediction done")
    return np.argmin(errs)
