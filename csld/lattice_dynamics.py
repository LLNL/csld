#!/usr/bin/env python3

"""
  Lattice dynamics model
"""

import os
from itertools import product
import numpy as np
from numpy.linalg import norm
import scipy.linalg
import scipy.io
import scipy.sparse
from scipy.sparse import lil_matrix as spmat

from .basic_lattice_model import BasicLatticeModel
from .interface_vasp import Poscar
from .cluster import Cluster
from .structure import SupercellStructure
from .util.mathtool import MixedIndexForImproper, IntegerDigits, tensor_constraint, RMS
from .util.tool import pad_right, matrix2text
from _c_util import fct_trans_c, ld_get_correlation, get_nullspace, init_ldff_basis, ldff_get_corr
from .coord_utils import ReadPBC2Cart
from .util.string_utils import str2arr, str2bool
from f_phonon import f_phonon

import logging
logger = logging.getLogger(__name__)

debug_level = 10

class LDModel(BasicLatticeModel):
    """
    Lattice dynamics
        atomic variables: vector displacements on each atom

        property:   total potential energy (electronic + nucleic potential, but NOT kinetic of nucleic)
                    force
    """
    def __init__(self, prim, raw_clusters, **kwargs):
        """

        :param prim:
        :param raw_clusters:
        :return:
        """
        BasicLatticeModel.__init__(self, prim, raw_clusters, **kwargs)
        self.all_range={o: max(self.proper_range[o],self.imp_range[o]) for o in self.proper_range}
        self.dipole_force = False if (not kwargs.get('dipole_force', True)) or \
                                          ('born_charge' not in self.prim[0].properties.keys()) or \
                                          (self.prim.intensive_properties['epsilon_inf'] is None) else True
        self.symm_residual_force = kwargs.get('symm_residual_force', True)
        self.Cmat = None
        self.ldff = None
        self._dpcor = None


    def translational_invariance(self):
        """
        Translational invariance Matrix
        :return: C matrix
        """
        print("Applying translational invariance.")
        crange = self.all_range
        prim = self.prim
        orbits = self.orbits
        rot_mats = [s.rot for s in prim.spacegroup]
        maxnpt = max([orb.cluster.order for orb in orbits])
        Cmat = self.Cmat1
        BmatCollect = []
        #self.prim._init_all_nb_list([crange[ord] for ord in range(2, 1+max(crange.keys()))])

        for iorb in range(len(orbits)):
            orb =  orbits[iorb]
            clus = orb.cluster
            if clus.order >= maxnpt:
                continue
            npt_ex = clus.order + 1

            cut = crange[npt_ex]
            if debug_level > 5555:
                print("debug> order %d cutoff %f" %(npt_ex, cut))
            if clus.diameter > cut:
                if debug_level > 9990:
                    print("  %d diameter %.4f out of range %.4f" % (iorb, clus.diameter, cut))
                continue

# find all points within range of clus
            if clus.order <=0:
                # empty cluster
                sumPts = [[0,0,0, l] for l in range(prim.num_sites)]
            else:
                sumPts = prim.find_nb_cluster(np.array(clus.ijkls), crange[npt_ex])
            if debug_level > 55:
                print("debug> translation of", repr(clus))

            dimTensor = 3** npt_ex
            Bmat = spmat((dimTensor, self.nfct_tot))
            foundClus = False
            for sumpt in sumPts:
                clusSum = clus.append_site(sumpt)
                if debug_level > 5555:
                    print("  ", iorb, " searching ", sumpt)
                #idSum= Cluster2MI[clusSum]//Sort;

# find the orbit that each summed cluster corresponds to *)
                [found, ioF, icF, igF, pi] = self.identify_cluster(clusSum)
                # print("identified",[found, ioF, igF, pi])
                if found:
                    foundClus = True
                    Bmat[:, self.orb_idx_full[ioF]:self.orb_idx_full[ioF+1]]+=\
                        fct_trans_c(npt_ex, 3, rot_mats[igF], pi)
            if not foundClus:
                print('  ',iorb, " nothing found for ", clus)
                continue

            BmatCollect.append(Bmat)
            Bmat = Bmat.dot(Cmat.T)
            # Bmat[abs(Bmat)< 1E-9] = 0

            if npt_ex > 999999:
                print("bypassing ", clus)
     #           BmatCollect.extend( Select[RowReduce[Chop[Bmat,10.^-10]],(Norm[#]> 10.^-6)&]])
            else:
                if debug_level > 99999:
                    print(" calc null bmat")
                    print(Bmat)
                # if not (scipy.sparse.isspmatrix(Bmat) and Bmat.getnnz()<=0):
                # Cmat = nullspace_rref(Bmat.toarray()).dot(Cmat)
                Cmat = get_nullspace(Bmat).dot(Cmat)
                print("  %4d + sum(a) %d remaining" % (iorb, Cmat.shape[0]))

        self.Cmat = Cmat
        return BmatCollect


    def symmetrize(self):
        self.isotropy_derivative_constraint()
        self.translational_invariance()
        #self.CheckNumericTranslationalInvariance()
    #    self.process_fct_order()


    def prepare_index_full(self):
        for orb in self.orbits:
            orb.ncorr_full = 3**orb.cluster.order
        super().prepare_index_full()
        self.nfct_tot = self.ncorr_full

    def prepare_index(self):
        # self.process_all_fct()
        self.nfct = self.Cmat.shape[0]
        allfct_ord = np.hstack([np.full(3**o.cluster.order, o.cluster.order, dtype=np.int) for o in self.orbits])
        self.fct_ord = [allfct_ord[row[0]] for row in self.Cmat.tolil().rows]
        self.ord_range = {o: len(self.fct_ord) - self.fct_ord[::-1].index(o) - 1 if o in self.fct_ord else 0
                          for o in range(self.maxorder+1)}
        self.ord_range[-1] = -1
        self.fct_ord = np.array(self.fct_ord)
        self.allfct_orbidx = np.hstack([np.full(3**o.cluster.order, i,dtype=int) for i, o in enumerate(self.orbits)])
        self.fct_orbidx = [self.allfct_orbidx[row[0]] for row in self.Cmat.tolil().rows]
        np.savetxt('num_fct_ord.txt', [self.ord_range[o]-self.ord_range[o-1] for o in range(self.maxorder+1)], fmt='%d')


    def isotropy_derivative_constraint(self):
        """
        Apply isotropy and derivative commutativity constraints for all clusters in this model
        :return:
        """
        clusters, iso_list, pi_list = zip(*[[orb.cluster, [x[0] for x in orb.isotropy[1:]],
                                             [x[1] for x in orb.isotropy[1:]]] for orb in self.orbits])
        Cmats = self.calc_isotropy_derivative_constraint(self.prim.spacegroup, clusters, iso_list, pi_list, self.nfct_tot, self.symm_residual_force)
        self.Cmat1 = scipy.sparse.vstack(Cmats)
        print("Isotropy/deriv. constraints done. After dim/before=", self.Cmat1.shape)
        return Cmats

    @staticmethod
    def calc_isotropy_derivative_constraint(ops, clusters, iso_list, pi_list, nfct_tot, symm_residual_force=True):
        """
        Apply isotropy and derivative commutativity constraints for all clusters in this model
        symm_residual_force: Whether to symmetrize the point cluster, i.e. residual forces. Default to True. If fitting phonons with small displacement and supercells that do not preserve symmetry (i.e. non-cubic low symmetry supercell for cubic system)
        :return:
        """
        print("Applying point group symmetry")
        Cmats = []
        ltothis=0

        for iorb in range(len(clusters)):
            clus = clusters[iorb]
            npt = clus.order
            dimThis = 3**npt
            ltothis+=dimThis
            if npt <= 0:
                null = scipy.sparse.identity(dimThis)
            else:
                idx_constr = MixedIndexForImproper(clus.vertices, 3)[0]
                igs = iso_list[iorb]
                pis = pi_list[iorb]
                if (not symm_residual_force) and (npt==1):
                    print("WARNING: symmetrization of point cluster (i.e. residual force) turned OFF")
                    igs = igs[:0]
                    pis = pis[:0]
                # print('debug isotropy2 igs, pis ', igs, list(igs), pis, 'zipped', *zip(*orb.isotropy[1:]))
                null = tensor_constraint(3, npt, [ops[ig].rot for ig in igs], pis, other_constraits=idx_constr)
            nfree = null.shape[0]
            print("  %4d null= %d/%d" %(iorb, nfree, dimThis), repr(clus), end='')
            if nfree>0:
                if npt <=2 and debug_level > 99999:
                    print([clus, null])
                if ltothis-dimThis>0:
                    null = scipy.sparse.bmat([[spmat((nfree, ltothis-dimThis)), null]])
                if nfct_tot-ltothis>0:
                    null = scipy.sparse.bmat([[null, spmat((nfree,  nfct_tot-ltothis))]])
                Cmats.append(null)
                print()
            else:
                print(' vanishing cluster!')
        return Cmats


    @staticmethod
    def write_mat(Cmat, outf):
        """
        Export Cmat
        :param outf:
        :return:
        """
        if Cmat is None:
            raise ValueError("Cmat not set for this model")
        print("writing matrix", Cmat.shape, "to", outf)
        scipy.io.mmwrite(outf, Cmat)

    @property
    def ncorr(self): return self.nfct + (0 if self.ldff is None else self.ldff.ncorr)

    def get_submodels(self, name_ord, u_list, lr_pair_penalty=0, ldffscale_list=[1], knownsol=None):
        """
        :param name_ord: list of [name, fct order]
        :param u_list: list of uscale
        :param lr_pair_penalty: penalty exp(-penalty*radius) for pair clusters
        :return: list of [name, matrix] defining the different fittings
        """
        all_ord=list(set([o.cluster.order for o in self.orbits]))
        ldff = self.ldff
        if ldff is None:
            ld_diag = np.arange(0)
            ld_scalelist = [1]
            maxfitord = self.maxorder
        else:
            ld_diag = np.ones(ldff.ncorr)
            ld_scalelist = ldffscale_list
            maxfitord = self.maxorder+1
            all_ord+=[maxfitord]
            self.ord_range[maxfitord] = self.ord_range[self.maxorder] + ldff.ncorr

        sol0 = np.zeros(self.ncorr)
        if knownsol:
             print("  Reading previous solution from %s"%(knownsol))
             input_sol = self.load_solution(knownsol).reshape(-1)
             sol0[:min(sol0.size, input_sol.size)] = input_sol[:min(sol0.size, input_sol.size)]

        print("No. of parameters", {o: self.ord_range[o]-self.ord_range[o-1]
                                              for o in range(maxfitord+1)})
        name_ord = self.process_name_ord(name_ord, all_ord)

        pair_r0 = np.min([orb.cluster.diameter for orb in self.orbits if orb.cluster.order==2 and orb.cluster.order_uniq==2])
        pair_diameter=np.array([self.orbits[idx].cluster.diameter-pair_r0 if self.orbits[idx].cluster.order==2 and self.orbits[idx].cluster.order_uniq==2 else 0 for idx in self.fct_orbidx])
        return [[nm+ ' uscale= %g'%(uscale) + str('' if ldff is None else " ldffscale=%g"%(ldffscale)),
                scipy.sparse.diags(np.hstack(((1/uscale)**(self.fct_ord-1)* np.exp(-pair_diameter*lr_pair_penalty), ldffscale*ld_diag)), 0).tocsr()[:,
                  sum([list(range(self.ord_range[i-1]+1, self.ord_range[i]+1)) for i in o if i<=maxfitord], []) if o[0]>=0 else list(range(-o[0],-o[1]))], sol0]
                for (nm, o) in name_ord for uscale in u_list for ldffscale in ld_scalelist]

    def CheckNumericTranslationalInvariance(self, trans=np.array([1.0, 2.0, 3.0])):
        """
        Apply uniform translation, calculate the force
        :return:
        """
        print("  To be implemented: checking translation", trans)


    def get_full_fct(self, sol_sym):
        """
        Return all FCT elements from symmetry reduced parameters
        :param sol_sym: symmetrized solution vector
        :return: expanded, full FCT's without symmetry
        """
        return self.Cmat.T.dot(sol_sym[:self.nfct])

    def get_correlation(self, sclist, wtFunc= lambda x:1, corrtype='f', delForce=1, shift=True, **kwargs):
        """
        :param sclist: [ [sc0, [sub0, sub1, ...]], [sc1, [sub11, ...]]]
        :param wtFunc:
        :param corrTyp: 'e' for energy, 'f' for force
        :param delForce: which components to delete
        :param shift: whether to subtract the shift (average force)
        :param residual_force: Whether to subtract residual forces of equilibrium structure, if found
        :return: correlation matrix A
        """
        import os.path
        ldff = self.ldff
        theC = self.Cmat.T
        ncorr = theC.shape[1]
        residual_force= str2bool(kwargs['setting'].get('residual_force', 'F'))
        if ldff is not None:
            ncorr += ldff.ncorr

        if corrtype == 'e':
            totNF = sum([len(sc[1]) for sc in sclist])
        elif corrtype == 'f':
            # ignore the last atom because of translational invariance
            totNF = sum([3*(Poscar.from_file(rd+"/POSCAR").structure.num_sites - delForce) for sc in sclist for rd in sc[1]])
        else:
            raise ValueError("ERROR: expecting to fit f(orce) or e(energy) but found %s"%(corrtype))
        print("  Total number of linear equations", totNF)
        assert totNF>0, ValueError("ERROR got no input data")

        Alist = np.zeros((totNF, ncorr))
        Flist = np.zeros((totNF, 3))
        Astart=0
        for sc in sclist:
            print("  reading supercell", sc[0])
            SCinfo= SupercellStructure.from_file(self.prim, sc[0])
            SCinfo.to_unperturbed()
            x0frac = SCinfo.frac_coords
            #SCinfo = SupercellStructure(self.prim, SCmat, None, x0frac)
            ncell = SCinfo.n_cell

            clusALL = self.translate_to_supercell(SCinfo)
            if corrtype=='f' and residual_force and os.path.exists(os.path.dirname(sc[0])+"/residual_force.txt"):
                print("  found 'residual_force.txt'")
                f0 = np.loadtxt(os.path.dirname(sc[0])+"/residual_force.txt")
                if shift:
                    f0-= np.mean(f0, axis=0)
            else:
                f0 = 0

            if self.dipole_force:
                fcmfile= os.path.dirname(os.path.abspath(sc[0]))+"/fcm_dp"
                if False and os.path.isfile(fcmfile):
                    print("    reading dipole FC for supercell "+fcmfile)
                    fcm_dp = np.loadtxt(fcmfile)
                else:
                    print("   computing long-range forces")
                    fcm_dp = self.get_hessian_dipole_corrected(SCinfo)
                    np.savetxt(fcmfile, fcm_dp)
            else:
                fcm_dp = None

            if ldff is not None:
                radialC= self.translate_to_supercell(SCinfo, ldff.orb_idx)
            if debug_level > 10:
                print("supercell clusters generated")

            for rd in sc[1]:
                rundir = rd
                dx= ReadPBC2Cart(rundir + "/POSCAR", x0frac)
                # weight= wtFunc(dx, uScale)
                # weight = np.ones(dx.shape[0])
                weight = 1
                if debug_level > 2:
                    print("    config",rundir, " weight=", weight, " max |dx|=", np.amax(norm(dx,axis=1)))
                dx_sort = dx
                Amat= self.calc_correlation(dx_sort, clusALL).dot(theC)

                if ldff is not None:
                    Amat = scipy.sparse.hstack((Amat, ldff.calc_correlation(dx_sort, radialC, ncell))).tocsr()

                if corrtype == "e":
                    thisNF = 1
                    if os.path.isfile(rundir+"/energy.txt"):
                        values = np.loadtxt(rundir+"/energy.txt", ndmin=1)/ncell
                    else:
                        print("WARNING: no energy.txt file found. Proceeding with 0...")
                        values = np.zeros(1)
                    valuesFit = values.copy()
                    if fcm_dp is not None:
                        en_dp = 0.5*np.dot(dx_sort.reshape(-1),fcm_dp.dot(dx_sort.reshape(-1)))
                        np.savetxt(rundir+"/energy.txt_dp", en_dp)
                        valuesFit -= en_dp/ncell
                    Amat = Amat[-1:]/ncell
                elif corrtype == 'f':
                    thisNF = 3*(len(dx)- delForce)
                    if os.path.exists(rundir+"/force.txt"):
                        values = np.loadtxt(rundir+"/force.txt")
                        shift_size= np.linalg.norm(np.sum(values, axis=0))
                        if shift_size >1e-3:
                            print("WARNING: large shift in force %.4f in %"(shift_size, rundir+"/force.txt"))
                        if shift:
                            values-= np.mean(values, axis=0)
                        valuesFit = values.copy() - f0
                    else:
                        print("WARNING: force.txt not found in %s ... setting to zero"%(rundir))
                        values = np.zeros((len(dx), 3))
                        valuesFit = values.copy()
                    assert values.shape == dx.shape, 'force [%d %d] coords [%d %d]' % (values.shape[0], values.shape[1], dx.shape[0], dx.shape[1])
                    if fcm_dp is not None:
                        f_dp = -fcm_dp.dot(dx_sort.reshape(-1)).reshape((-1,3))
                        #np.savetxt(rundir+"/force.txt_dp", f_dp)
                        valuesFit -= f_dp
                    values = values.flatten()[:thisNF]
                    valuesFit = valuesFit.flatten()[:thisNF]
                    if debug_level >30:
                        print("forces read in")
                    if thisNF != len(values):
                        raise ValueError("ERROR: expecting ", thisNF, " but found ", len(values), " force components")
                    Amat = Amat[:thisNF]

                if debug_level >9999:
                    print("      A size", Amat.shape, Amat.__class__)
                Alist[Astart:Astart+thisNF, :] = (Amat * weight).todense()
                Flist[Astart:Astart+thisNF, 0] = valuesFit * weight
                Flist[Astart:Astart+thisNF, 1] = np.full((thisNF), weight, dtype=np.double)
                Flist[Astart:Astart+thisNF, 2] = values-valuesFit
                Astart += thisNF

        return [spmat(Alist), Flist]



    def calc_correlation(self, dx, clusALL):
        """

        :param dx:
        :param clusALL: all clusters in the supercell
        :return:
        """
        maxnpt = self.maxorder
        return spmat(ld_get_correlation(dx.shape[0], len(self.orbits), maxnpt, self.nfct_tot,
                                  np.array(dx),
                                  np.array([pad_right(np.array(clus[0]), maxnpt) for clus in clusALL], dtype=np.int32),
                                  np.array([clus[1] for clus in clusALL], dtype=np.int32),
                                  np.array([clus[2] for clus in clusALL], dtype=np.int32),
                                  np.array([orb.cluster.order for orb in self.orbits], dtype=np.int32),
                                  np.array([orb.cluster.factorial for orb in self.orbits]),
                                  np.array([op.rot_inv for op in self.prim.spacegroup])))

    def save_fct(self, sol, outf, scmat, combine_improper=True):
        """
        :param sol: solution vector
        :param outf: output filename. Two files .lat and .pot, will be generated
        :param scmat: supercell integer 3x3 matrix
        :return:
        """
        print("  saving lattice and potential to", outf)
        scinfo = SupercellStructure.from_scmat(self.prim, scmat)
        self.save_fct_lat(outf+'.lat', scinfo)
        assert sol.shape[0] >= self.nfct
        self.save_fct_pot(outf+'.pot', self.get_full_fct(sol), sol[self.nfct:], scinfo,
                          combine_improper=combine_improper)

    def save_fct_lat(self, outf, scinfo):
        """

        :param outf:
        :param scinfo:
        :return:
        """
        # write lattice points
        fp = open(outf, 'w')
        natom = self.prim.num_sites
        ncell = scinfo.n_cell
        SCposFrac=  scinfo.frac_coords
        # print("debug mass", self.prim.atomic_masses)
        outs =[matrix2text(self.prim.lattice._matrix), matrix2text(scinfo.sc_mat), str(natom)]
        outs += ["%f %f %f %d" % tuple(self.prim.frac_coords[i].tolist()+[self.prim.atomic_numbers[i]]) for i in range(natom)]
        outs += [str(SCposFrac.shape[0]), '']
        fp.write('\n'.join(outs))

        for iA in range(natom):
            for jLat in range(scinfo.n_cell):
                fp.write("%d %s %d %s\n" %(iA*ncell+ jLat+1, matrix2text([scinfo.sc_ref[jLat]]),
                        iA+1, matrix2text([SCposFrac[iA*ncell+ jLat]])))
        fp.close()

    def save_fct_pot(self, outf, sol_fct, sol_ff, scinfo, tol=1.E-12, combine_improper=False, output_ijkl=True):
        """
        write potential file.
        :param outf:
        :param scinfo:
        :return:
        """
        ldff = self.ldff
        fp = open(outf, 'w')
        natom = self.prim.num_sites
        dim = 3
        ops = self.prim.spacegroup
        fc_norm = []
        if self.dipole_force:
            fcm_dp = self.get_hessian_dipole_corrected(scinfo)
        else:
            fcm_dp = np.zeros((scinfo.num_sites*3, scinfo.num_sites*3))
        flag_dp = np.zeros((scinfo.num_sites, scinfo.num_sites),dtype=np.int)

        for iO, orb in enumerate(self.orbits):
            clus0 = orb.cluster
            npt = clus0.order
            fac = clus0.factorial
            val = sol_fct[self.orb_idx_full[iO]:self.orb_idx_full[iO+1]]
            fc_norm.append([clus0.diameter, np.linalg.norm(val),clus0.order, clus0.order_uniq])

            if ldff is not None:
                ppout = ldff.tostr(sol_ff, iO)
            else:
                ppout = "0\n0"

            if fc_norm[-1][1]>tol or (ppout != "0\n0"):
                for ic, clus in enumerate(orb.clusters):
                    trans_cluster = np.array(BasicLatticeModel.translate_cluster_to_supercell(scinfo, clus))
                    valTrans = fct_trans_c(npt, 3, ops[orb.clusters_ig[ic]].rot, np.arange(npt, dtype=int)).dot(val)
                    if npt==2:
                        valTrans+= fcm_dp[trans_cluster[0,0]*3:trans_cluster[0,0]*3+3,trans_cluster[0,1]*3:trans_cluster[0,1]*3+3].flatten()
                        flag_dp[trans_cluster.T.tolist()]=1
                    if combine_improper:
                        valTrans= clus.reduce_improper_fct_output(valTrans)  # fewer terms to save!
                    # fctTrans = valTrans.reshape([dim for _ in range(npt)])
                    # If[npt==2,AppendTo[pairc, orbitUniq[[iC,icOrb,1]]]; AppendTo[pairFCM, fctTrans]];
                    # if ic <= 3:
                    #     print(iO, ic, clus)
                    #     print(valTrans)
                    #     print(trans_cluster)
                    clus_out = matrix2text(clus._ijkls if output_ijkl else clus.frac_coords)
                    outs = [str(npt), clus_out, str(len(trans_cluster)), matrix2text(trans_cluster+1),
                            LDModel.fct2str(npt, valTrans/fac, tol), ppout]
                    fp.write("\n".join(outs) + "\n\n")
        np.savetxt("fct_norm_vs_diameter.txt", fc_norm, header='col1=diameter col2=norm col3=npt col4=npt_uniq')
        return
        for i1 in range(scinfo.num_sites):
            for i2 in range(i1, scinfo.num_sites):
                if flag_dp[i1,i2]:
                    continue
                npt = 2
                # WARNING: TODO: convert pair of coords to minimal distance within the supercell periodic boundary condition
                clus = Cluster.from_coords(scinfo.cart_coords[[i1,i2]], self.prim, frac_coords=False)
                fac = clus.factorial
                clus_out = matrix2text(clus._ijkls if output_ijkl else clus.frac_coords)
                ppout = "0\n0"
                valTrans = fcm_dp[i1*3:i1*3+3, i2*3:i2*3+3].flatten()
                outs = [str(npt), clus_out, str(1), matrix2text(np.array([[i1,i2]])+1),
                        LDModel.fct2str(npt, valTrans/fac, tol), ppout]
                fp.write("\n".join(outs) + "\n\n")


    def save_fcshengbte(self, sol, ord, tol=1e-20):
        assert ord in (3,4), "Only order 3 or 4 FCTs are accepted by shengbte, got %d"%(ord)
        import io, re
        sol_fct = self.get_full_fct(sol)
        ops = self.prim.spacegroup
        fc_name="FORCE_CONSTANTS_%s"%({3:"3RD",4:"4TH"}[ord])
        fp=io.StringIO()
        icount=0
        for iO, orb in enumerate(self.orbits):
            clus0 = orb.cluster
            npt = orb.order
            if npt != ord:
                continue
            val = sol_fct[self.orb_idx_full[iO]:self.orb_idx_full[iO+1]]
            if np.amax(np.abs(val)) <= tol:
                continue
            for ic, clus in enumerate(orb.clusters):
                ijkls = clus._ijkls_np
                valTrans = fct_trans_c(npt, 3, ops[orb.clusters_ig[ic]].rot, np.arange(npt, dtype=int)).dot(val).reshape((3,)*ord)
#                print('debug iorb, ic iper', iO, ic, len(perms))
                for iper in clus0.permutations():
                    icount+=1
                    #print('debug', icount, clus0.ijkls, iO, ic)
                    ijk_other= matrix2text(self.prim.lattice.get_cartesian_coords(ijkls[iper[1:],:3]- ijkls[iper[0:1],:3]))
                    valPerm = np.transpose(valTrans, iper).reshape((-1))
                    fp.write("\n%d\n%s\n%s\n"%(icount, ijk_other, matrix2text(ijkls[iper,3]+1)))
                    fp.write(re.sub(r".*\n", r"",LDModel.fct2str(npt, valPerm, -1),count=1)+'\n')
        with open(fc_name, 'w') as modified: modified.write("%d\n"%(icount) + fp.getvalue())
        fp.close()

    def load_solution(self, sol_f, potential_coords_ijkl=True):
        """
        sol_f: file_name_of_solution [order_to_keep]
                  File format is either solution vector or potential file
                  order_to_keep is like 0,1,2 (no space)
        """
        solinf = sol_f.split()
        if solinf[0][-4:].lower()!='.pot':
            print("  Loading symmetrized FCT from %s"% (solinf[0]))
            sol= np.loadtxt(solinf[0], ndmin=2)
            if len(solinf)>1:
                print("WARNING!!!!  only order %s will be kept"%(solinf[1]))
                tmp = np.zeros_like(sol)
                for ord in eval("[%s]"%(solinf[1])):
                    print("ord= %d  corresponding idx="%(ord), (self.ord_range[ord-1]+1,self.ord_range[ord]+1))
                    tmp[:,self.ord_range[ord-1]+1:self.ord_range[ord]+1] = sol[:,self.ord_range[ord-1]+1:self.ord_range[ord]+1]
                sol=tmp
            return sol
        else:
            from .util.io_utils import read_nrecord_array
            print("  Loading symmetrized FCT from potential %s"% (solinf[0]))
            full_fct= np.zeros(self.nfct_tot)
            lines = open(solinf[0], 'r').readlines()
            line=0
            while line<len(lines):
                line, xyz=read_nrecord_array(lines, line)
                if potential_coords_ijkl:
                    clus= Cluster.from_ijkl(xyz.astype(int), self.prim)
                else:
                    clus= Cluster.from_coords(xyz, self.prim)
                line, clus_instances=read_nrecord_array(lines, line)
                line, ijval=read_nrecord_array(lines, line)
                line, rad1=read_nrecord_array(lines, line)
                line, rad2=read_nrecord_array(lines, line)
                line += 1 # empty line
                [found, ioF, icF, igF, pi] = self.identify_cluster(clus)
                if (not found) or (icF != 0) or (igF != 0):
                    continue
                ord = clus.order
          #      print("found cluster order=%d id=%d line=%d"%(ord, ioF, line))
                fct= np.zeros((3,)*ord)
                for x in ijval:
                    #print("  debug x=", x, tuple(x[:ord].astype(int)-1), x[-1])
                    fct[tuple(x[:ord].astype(int)-1)] = x[-1]
                fct = fct*clus.factorial
                fct = fct_trans_c(ord, 3, self.prim.spacegroup[igF].rot, pi).T.dot(fct.reshape((-1)))
                full_fct[self.orb_idx_full[ioF]:self.orb_idx_full[ioF+1]] = fct.reshape((-1))

            #print("debug full_fct=", full_fct)
            sol = scipy.sparse.linalg.lsqr(self.Cmat.T, full_fct,atol=1e-20,btol=1e-20)
            #print("debug sol=", sol)
            np.savetxt(solinf[0]+'_loaded_sol', sol[0])
            if sol[3] > 1E-4:
                print("WARNING large error %f loading potential to symmetrized FCT"%(sol[3]))
            return np.array([sol[0]])





    def get_pair_info(self, sol_fct, ord=2, tol=1.E-20):
        """
        extract pair interactions for phonon calculations.
        :param ord: usually 2
        :param sol_fct: solution vector
        :return:
        """
        natom = self.prim.num_sites
        ops = self.prim.spacegroup
        pairijk = []
        pairTyp = []
        pairFCM = []
        dim=3

        for iO, orb in enumerate(self.orbits):
            npt = orb.cluster.order
            if npt != ord:
                continue
            val = sol_fct[self.orb_idx_full[iO]:self.orb_idx_full[iO+1]]
            if (abs(val)<=tol).all():
                continue

            for ic, clus in enumerate(orb.clusters):
                valTrans = fct_trans_c(npt, 3, ops[orb.clusters_ig[ic]].rot, np.arange(npt, dtype=int)).dot(val)
                fctTrans = valTrans.reshape([dim]*npt)
                pairijk.append(clus._ijkls_np[0,:3] - clus._ijkls_np[1,:3])
                pairTyp.append(clus._ijkls_np[:,3])
                pairFCM.append(fctTrans)
        if len(pairijk)>0:
            return (np.array(pairijk), pairTyp, pairFCM)
        else:
            return (np.zeros((1,3),dtype=int), np.zeros((1,2),dtype=int), np.zeros((1,3,3)))


    @staticmethod
    def fct2str(npt, fct, tol=1.E-12):
        """
        Note fct is a 1-D array, NOT tensor
        """
        outs = ["%s %.15f"%(matrix2text([IntegerDigits(i, 3, npt)+1]), fct[i]) for i in range(3**npt) if abs(fct[i]) > tol]
        return "\n".join([str(len(outs))] + outs)


    @staticmethod
    def get_hessian_dipole(s):
        """
        :param s:  structure with epsilon_inf and born_charge
        """
        if s.intensive_properties['epsilon_inf'] is None:
            return np.zeros((3*s.num_sites,3*s.num_sites))
        return f_phonon.get_fcm_dipole(s.lattice.matrix.T, s.lattice.inv_matrix, 1E-18, s.cart_coords.T,
                                       np.array(s.site_properties['born_charge']).transpose([1,2,0]), s.intensive_properties['epsilon_inf'].T, [0,0,0]).real


    def get_hessian_dipole_corrected(self, s):
        """
        s: supercell
        """
        fcm_dp = self.get_hessian_dipole(s)
        if self._dpcor is None:
            return fcm_dp
        f_phonon.init(s.lattice._matrix, s.atomic_masses, s.frac_coords, *self.translate_pairinfo_to_supercell(s, *self._dpcor))
        fcm_cor= f_phonon.get_dm([[0.,0.,0.]],3*s.num_sites)[:,:,0].real
        #print(fcm_cor.shape, fcm_dp.shape, fcm_cor[0,0], fcm_dp[0,0], s.atomic_masses.__class__, s.atomic_masses)
        #for iA in range(s.num_sites):
        #    for ix in range(3):
        #        for jA in range(s.num_sites):
        #            for jx in range(3):
        #                fcm_cor[iA*3+ix,jA*3+jx]*=np.sqrt(s.atomic_masses[iA]*s.atomic_masses[jA])

        #np.savetxt('fcm_dp', fcm_dp)
        #np.savetxt('fcm_cor', fcm_cor)
        mass= np.sqrt(np.array(s.atomic_masses).repeat(3))
        fcm_cor*= np.outer(mass, mass)
        #np.savetxt('fcm_corScaled', fcm_cor)
        #np.savetxt('fcm_all', fcm_dp+fcm_cor)
        return fcm_dp+ fcm_cor


    def translate_pairinfo_to_supercell(self, sc, ijk_prim, typ_prim, fcm_prim):
        return LDModel.pairinfo_to_supercell(self.prim, sc, ijk_prim, typ_prim, fcm_prim)

    @staticmethod
    def pairinfo_to_supercell(prim, sc, ijk_prim, typ_prim, fcm_prim):
        nfcm_prim = len(ijk_prim)
        nfcm = nfcm_prim*sc.n_cell
        ijk=np.zeros((nfcm, 3),dtype=int)
        typ=np.zeros((nfcm, 2),dtype=int)
        fcm=np.zeros((nfcm, 3,3))
        for i in range(nfcm_prim):
            clus= Cluster.from_ijkl([np.append(ijk_prim[i],typ_prim[i][0]), [0,0,0,typ_prim[i][1]]], prim)
            tc= BasicLatticeModel.translate_cluster_to_supercell(sc, clus)
            newijk= np.round(np.dot(clus.frac_coords[0]-clus.frac_coords[1], sc.inv_sc_mat)-(sc.frac_coords[tc[0][0]]-sc.frac_coords[tc[0][1]])).astype(int)
            # print(i, "cluster", clus, tc, newijk)
            for j in range(sc.n_cell):
                fcm[i*sc.n_cell+j] = fcm_prim[i]
                typ[i*sc.n_cell+j] = tc[j]
                ijk[i*sc.n_cell+j] = newijk
        return ijk, typ, fcm


    def get_dpcor(self, bondlen, errtol=1e-7):
        offd=[[0,0,1],[1,2,2]]
        offdflatUp = [1,2,5]
        offdflatDn = [3,6,7]
        diagflat= [0,4,8]
        fcm_dp = self.get_hessian_dipole(self.prim)
        np.savetxt("prim_dp", fcm_dp)
        non_symm = [fcm_dp[3*i:3*i+3, 3*i:3*i+3]-fcm_dp[3*i:3*i+3, 3*i:3*i+3].T for i in range(self.prim.num_sites)]
        pts=self.l_point_cls()
        npt=len(pts)
        bvec = np.array([non_symm[pt][offd[0],offd[1]] for pt in pts]).reshape((-1))
        if np.linalg.norm(bvec)/np.sqrt(len(bvec)) <1E-13:
            print("++ no ASR violation in long-range force constant matrix")
            return None
        print('************** corrections to long-range force constant matrix *****************')
        # create an LD model using nearest neighbor only
        ldNN = init_ld_model(self.prim, {'model_type':'LD', 'max_order':2, 'cluster_diameter':str(bondlen),
          'proper_diameter':str(bondlen),'cluster_filter':'lambda cls: True'}, {}, 2, 2, 0, False)
        print(ldNN)
        C1mats = ldNN.isotropy_derivative_constraint()[1+2*npt:]
        C1 = spmat(scipy.sparse.vstack(C1mats))[:,1+12*npt:]
        nvar= C1.shape[0]

        Bmats= ldNN.translational_invariance()[1:]
        Bmats = [i[:,1+12*npt:] for i in Bmats]
        B1 = spmat(scipy.sparse.vstack(Bmats))
        Acorrection = spmat(np.zeros((len(bvec), nvar)))
        for i, pt in enumerate(pts):
            Acorrection[3*i:3*i+3] = (Bmats[i][offdflatUp]-Bmats[i][offdflatDn]).dot(C1.T)
        Acorrection=spmat(Acorrection)

#        from cssolve.bregman_func import bregman_func
#        solution = bregman_func(Acorrection, bvec, method=1, mu=1E-5, lbd=3,maxIter=2000, tol=1E-6)
#        print(get_errors(bvec, Acorrection.dot(solution)))
#        print(solution)
        dpcor_sol_f='dpcor_sol.dat'
        if False and os.path.isfile(dpcor_sol_f):
            print('++ Loading dpcor from %s'%(dpcor_sol_f))
            solution = np.loadtxt(dpcor_sol_f)
        else:
            solution = scipy.sparse.linalg.lsqr(Acorrection, bvec)[0]
            np.savetxt(dpcor_sol_f, solution)
#        solution = np.linalg.lstsq(Acorrection[:-3].todense(), bvec[:-3])[0]
        rmse = RMS(bvec - Acorrection.dot(solution))
        if rmse > errtol:
            raise ValueError('dpcor correction FAILED rmse= %5g. Check symmetry or increase dpcor_bond'%(rmse))
        # np.savetxt('Adpcor.out', Acorrection.todense())
        # np.savetxt('bdpcor.out', bvec)
        #print('correction FCM=',solution)
        print('************** corrections done (rmse= %5g) *****************'%(rmse))
        # send correction SR FCM to f_phonon
        full_sol_proper= np.array(C1.T.dot(solution))
        #print(C1.T.dot(solution))
        #print(np.zeros(1+3*len(pts)), "\nonsite", -B1.dot(full_sol_proper), "\npair", full_sol_proper)
        # print('debug onsite dpcor', B1.dot(full_sol_proper).reshape((-1,3,3)))
        full_sol = np.hstack((np.zeros(1+3*len(pts)), -B1.dot(full_sol_proper), full_sol_proper))
        # print('debug trans inv', ldNN.translational_invariance())
        # print('debug checking trans inv', scipy.sparse.vstack(ldNN.translational_invariance()).dot(full_sol))
        #print('DEBUG CORR_pair_info',ldNN.get_pair_info(full_sol))
        dpcor_pair = ldNN.get_pair_info(full_sol,2,1E-30)
      #  f_phonon.init_dpcor(*dpcor_pair)
        self._dpcor = dpcor_pair
        return dpcor_pair



def init_ld_model(prim, setting, setting_ldff, clus_step, symC_step, ldff_step, dpcor=True, pdfout=None):
    """
     model initialization
    :param prim:
    :param setting:
    :param ldff_setting:
    :param clus_step:
    :param symC_step:
    :param ldff_step:
    :return: LD model and its associated LDFF model
    """
    from scipy.io import mmread, mmwrite

    if clus_step <= 0:
        exit(0)
    assert setting['model_type'] == 'LD', ValueError("This script is intended for lattice dynamics only")
    maxorder = int(setting['max_order'])
    scale = prim.lattice._scale if str2bool(setting.get('fractional_distance','False')) else 1
    irange = (np.hstack(([0.1, 0.1], str2arr(setting['cluster_diameter'])[:maxorder-1]))*scale).tolist()
    prange_str = setting.get('proper_diameter', '')
    # use cluster_diameter if proper_diameter not specified
    if not prange_str:
        prange_str = setting['cluster_diameter']
    prange = (np.hstack(([0.1, 0.1], str2arr(prange_str)[:maxorder-1]))*scale).tolist()
    irange = dict(zip(range(maxorder+1), irange))
    prange = dict(zip(range(maxorder+1), prange))
    clus_sel = eval(setting['cluster_filter'])
    dipole_force = str2bool(setting.get('dipole_force', 'True'))
    symm_residual_force = str2bool(setting.get('symm_residual_force', 'True'))
    spec = {'maxorder':maxorder, 'prange':prange, 'filter':clus_sel, 'dipole_force':dipole_force, 'symm_residual_force':symm_residual_force}
    spec.update({'irange':irange})
    if clus_step == 1:
        model = LDModel.from_file(prim, setting['cluster_in'], **spec)
    elif clus_step in [2, 3]:
        model = LDModel.generate_clusters(prim, **spec)
        model.cleanup()
        if clus_step == 3:
            model.save_clusters(setting['cluster_out'])
    else:
        print("ERROR: Unknown clus_step: ", clus_step)
        exit(-1)
    print("+ Obtained %d proper clusters" %(len(model.clusters)), model.tally())
    model.generate_improper()
    model.cleanup()
    model.get_orbit_isotropy()
    model.prepare_index_full()

    # if we only need a NN model
    if not dpcor:
        return model
    if model.dipole_force and model._dpcor is None:
        model.get_dpcor(setting.getfloat('dpcor_bond', 2.8),setting.getfloat('dpcor_errtol', 1e-7))

    print(model)
    #model.save_clusters('cluster_all')

######## independent parameters
    if symC_step <= 0:
        exit(0)
    elif symC_step == 1:
        model.Cmat = mmread(setting['symC_in'])
    #    ld.process_fct_order(pdfout)
    elif symC_step in [2, 3]:
        model.symmetrize()
        if symC_step == 3:
            mmwrite(setting['symC_out'], model.Cmat)
    else:
        print("ERROR: Unknown symC_step: ", symC_step)
        exit(-1)
    model.prepare_index()
    print("+ LD symmetrization done. After dim/before=", model.Cmat.shape)

######## Force field on lattice
    if len(setting_ldff) <= 0:
        model.ldff = None
        return model
    entries = [k for k, v in setting_ldff.items()]
    if (ldff_step <= 0) or ('orbit_indices' not in entries):
        model.ldff = None
    elif ldff_step == 2:
        l234 = list(map(int, setting_ldff['num_basis'].split()))
        assert len(l234) >= 1
        xpts = list(map(float, setting_ldff['interpolation_pts'].split()))
        assert len(xpts) >= 3
        if 'polaron_force' in entries:
            ldfftype= PolaronFF
            model.ldff = PolaronFF(model, str2arr(setting_ldff['orbit_indices'], int).tolist(),
                         xpts=np.arange(*xpts[:3]),
                         lmax2=l234[0],
                         bas2 = eval(setting_ldff['basis_2']),
                         nradial=int(setting_ldff['nradial']),
                         dimer_indices=str2arr(setting_ldff['dimer_indices'],int).tolist(),
                         chgFunc=eval(setting_ldff['chgfunc']),
                         dchgFunc=eval(setting_ldff['dchgfunc']))
        else:
            ldfftype= LDFFmodel
            model.ldff = LDFFmodel(model, str2arr(setting_ldff['orbit_indices'], int).tolist(),
                         xpts=np.arange(*xpts[:3]),
                         lmax2=l234[0],
                         cut12=str2arr(setting_ldff.get('cut12','-0.7 0.7'),float,(-1,2)),
                         m12=str2arr(setting_ldff.get('m12','12 6'),int,(-1,2)),
                         bas2 = eval(setting_ldff['basis_2']))
        print("+ %s initialized %d parameters" %(ldfftype.__name__, model.ldff.ncorr))
    else:
        print("ERROR: Unknown ldff_step: ", ldff_step)
        exit(-1)
    return model


class LDFFmodel():
    """
    LD force field (only depends on interatomic distances)
    cut12: cutoffs for extrapolation. if r-r0<cut1 or >cut2, extrapolate to c0+c1/r**m
    m12: m1, m2 for cut1, cut2 respectively
    """
    def __init__(self, ld, orb_idx, lmax2=-1, bas2=[], lmax3=-1, bas3=[], xpts=np.array([]), cut12=np.array([[-0.7,0.7]]), m12=np.array([[12,6]])):
        """
        :param ld  the LD model
        :param orb_idx  indices of the selected orbits used in LDFF
        :param bas2  Either 1) a list of basis functions, each takes 1 parameter, dr=r-r0
                            2) a function b[l, dr] where l=0..lmax-1
        :param lmax2  is ignored if bas2 is a list
        :param xpts  a list of sampling points for dr, e.g. -1, -0.9, ..., 1
        """
        self.ld = ld
        self.orb_idx = orb_idx
        self.lmax2 = len(bas2) if isinstance(bas2, list) else lmax2
        self.bas2 = bas2
        self.lmax3 = len(bas3) if isinstance(bas3, list) else lmax3
        self.bas3 = bas3
        self.xpts = xpts
        self.cut12=cut12
        self.m12=m12
        n_xpts = len(xpts)
        ncorr_list = []
        multi_list = []
        ffidx_list = []
        npt_list = []
        ncorr = 0
        for i in orb_idx:
            orb = ld.orbits[i]
            npt = orb.cluster.order
            assert npt == 2, TypeError("Pair interactions for LDFF only")
            assert orb.cluster.factorial <= 1, TypeError("LDFF cannot accept improper orbit %d"%(i))
            # nvar = npt*(npt-1)/2
            ffidx = self.symmetrize_idx(npt)
            multi = np.array([x.shape[0] for x in ffidx])
            nc = len(ffidx)
            ncorr += nc
            ncorr_list.append(nc)
            multi_list.append(multi)
            ffidx_list.append(ffidx)
            npt_list.append(npt)
        self.ncorr_list = np.array(ncorr_list)
        self.ffidx_list = ffidx_list
        self.multi_list = np.array(multi_list)
        self.npt_list = npt_list
        self.ncorr = ncorr
        y2 = np.zeros((self.lmax2, n_xpts))
        for l in range(self.lmax2):
            for ix, x in enumerate(xpts):
                y2[l, ix] = self.eval_bas(bas2, l, x)
        np.savetxt('ldff_bas.txt', np.vstack((xpts,y2)).T)
        self.y2 = y2
        init_ldff_basis(2, self.lmax2, xpts, y2)

    @staticmethod
    def eval_bas(bas, l, x):
        return bas[l](x) if isinstance(bas, list) else bas(l, x)

    @staticmethod
    def eval_val(bas, ppval, x):
        return np.dot(ppval, [LDFFmodel.eval_bas(bas, l, x) for l in range(len(ppval))])
    
    def symmetrize_idx(self, npt):
        if npt == 2:
            return np.arange(self.lmax2)[:,None,None]
        elif npt == 3:
            print("3-body LDFF symmetrization TO BE IMPLEMENTED")
            return [[list(i)] for i in product(range(self.lmax3), range(self.lmax3), range(self.lmax3))]


    def tostr(self, sol, io, tol=1.E-12):
        """

        :param sol: the whole LDFF coefficients
        :param io:
        :return:
        """
        if io not in self.orb_idx:
            return "0\n0"
        iff = self.orb_idx.index(io)
        npt = self.ld.orbits[io].cluster.order
        lm = self.lmax2 if npt == 2 else self.lmax3
        ppval= sol[iff*lm:(iff+1)*lm]
        ppord= list(range(1,lm+1))
        outs = []
        # if (abs(ppval)> tol).any():
#           ppval=Transpose[{ppord, ppval}];
#           ppval=Select[ppval, (Abs[#[[2]]]>0)&];
        npp = 0
        for i in range(lm):
            if abs(ppval[i])>tol:
                npp += 1
                outs.append("%d %.12f" % (ppord[i], ppval[i]))
        outs.append("1\n"+ " ".join(map(str,self.extrapolations(iff, ppval))))
        return ("%d\n" % (npp)) + "\n".join(outs)


    def extrapolations(self, iff, ppval, dx=1E-4):
        exff = list(range(8))
        xfrac = self.ld.orbits[self.orb_idx[iff]].cluster.frac_coords
        r0 = np.linalg.norm(self.ld.prim.lattice.get_cartesian_coords(xfrac[1]-xfrac[0]))
        for i in range(2):
            exff[i*3] = self.cut12[0,i] if len(self.cut12)<=1  else self.cut12[iff,i]
            xa = exff[i*3]
            xb = xa+dx
            ya = self.eval_val(self.bas2, ppval, xa)
            yb = self.eval_val(self.bas2, ppval, xb)
            m= self.m12[0,i] if len(self.m12)<=1  else self.m12[iff,i]
            r=r0+xa
            exff[i*3+2] = -(yb-ya)/dx*(r**(m+1))/m
            exff[i*3+1] = ya - exff[i*3+2]/(r**m)
            exff[6+i] = m
        return exff


    def calc_correlation(self, dx, clusALL, ncell):
        """

        :param dx:
        :param clusALL: all clusters in the supercell
        :return:
        """
        len_orb = np.array([len(self.ld.orbits[i].clusters) for i in self.orb_idx])
        len_orb_sums = [len_orb[:i].sum() for i in range(len(self.orb_idx))]
        clus_id = [len_orb_sums[i] + j for i, ii in enumerate(self.orb_idx)
                    for j in range(len(self.ld.orbits[ii].clusters)) for _ in range(ncell)]
        return spmat(ldff_get_corr(
            np.array([self.ld.orbits[i].cluster.order for i in self.orb_idx], dtype=np.int32),
            np.array([self.lmax2, self.lmax3, 0], dtype=np.int32),
            np.array(self.ncorr_list, dtype=np.int32),
            np.array([ii for i in self.multi_list for ii in i], dtype=np.int32),
            np.array([iiii for i in self.ffidx_list for ii in i for iii in ii for iiii in iii], dtype=np.int32),
            np.array([pad_right(clus.coords,[4,3]) for i in self.orb_idx for clus in self.ld.orbits[i].clusters]),
            np.array(dx),
            np.array([self.orb_idx.index(clus[1]) for clus in clusALL], dtype=np.int32),
            np.array(clus_id, dtype=np.int32),
            np.array([pad_right(np.array(clus[0]), 4) for clus in clusALL], dtype=np.int32)))

    def plot_pairPES(self, sols):
         fname= 'ldff_PES.txt'
         header= 'col1=du'
         col=2
         mat=[self.xpts]
         for isol, sol0 in enumerate(sols):
             offset=0
             sol = sol0[-self.ncorr:]
             for i,npt in enumerate(self.npt_list):
                 if npt==2:
                     mat.append(np.dot(sol[offset:offset+self.lmax2], self.y2))
                     header+=" %d=sol_%d_clus_%d"%(col, isol+1, self.orb_idx[i]+1)
                     col+=1
                 offset+= self.ncorr_list[i]
         np.savetxt(fname, np.array(mat).T, header=header)
         print("  LDFF: pair PES exported to %s"%(fname))


