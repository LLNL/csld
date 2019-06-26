#!/usr/bin/env python3

"""
Phonon

Only harmonic Phonon subroutines implemented
All numerically heavy computations are in f90 (f_phonon.f90)
"""

#import os
#from itertools import islice, product
from ..util.tool import my_flatten
from ..util.string_utils import fn_available
from ..symm_kpts import HighSymmKpath
from ..interface_vasp import Poscar
#from ..structure import SupercellStructure

import numpy as np
from math import sqrt, pi
try:
    # import matplotlib
    # matplotlib.use('AGG')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
except ImportError:
    print("WARNING: cannot import pyplot; plotting disabled")
    plt = None
    pass

from f_phonon import f_phonon

#import logging
#logger = logging.getLogger(__name__)

#debug_level = 10

class NA_correction():
    """
    long range non-analytic dipole-dipole correction for force-constants based phonon calculation
    """
    def __init__(self, method, rho, Ncell=-1.0):
        """
        method: -1 (disabled), 1 (Parlinski) or 2 (mixed-space) or 0 (dipole, FC MUST exclude dipole force)
        rho: mixing parameter in Parlinski
        Ncell: number of cells in the mixed-space approach
        """
        self.id= method
        self.rho = rho
        self.Ncell = float(Ncell)
        self.uniq_ijk_idx = []

    @classmethod
    def from_dict(cls, d):
        rho = d.getfloat('NAC_rho', 0.05)
        Ncell = d.getfloat('NAC_Ncell', -1.0)
        return cls(d.getint('nac',-1), rho, Ncell)


class Phonon():
    """

    All phonon related stuff
    Note: all q-points (and hence real space coordinates) should be Cartesian unless explicitly specified

    """
    def __init__(self, prim, LD, sol, pdfout, NAC=None, etafac=8.0):
        """
        :param prim:
        :param LD:
        :param sol: solution vector
        :param NAC: long range dipole-dipole non-analytic correction
        :return:
        """
        self.units = {'THz': (1, 'Freq. (THz)'), 'meV': (4.1356668, 'En (meV)'), 'eV': (0.0041356668, 'En (eV)'),
                      'cm': (33.3564227583, ' Freq. (cm^{-1})')}
        self.prim = prim
        self.LD = LD
        #self.dim = prim.num_sites*3
        if prim.intensive_properties['epsilon_inf'] is None:
            NAC=None
        if NAC is not None:
            if (NAC.id != 0) and LD.dipole_force:
                raise ValueError('Must set nac=0 when dipole forces are considered explicitly')
        self.NAC = NAC
        self.pdfout = pdfout
# process pair interactions
        self.pairinfo = LD.get_pair_info(LD.get_full_fct(sol))
        self.setup_cell(self.prim)

    def setup_cell(self, cell):
        from ..util.mathtool import Union
        self.cell=cell
        self.dim=self.cell.num_sites*3
        #print("debug tralslate pair to supercell, size=", cell.num_sites)
        pairinfo = self.pairinfo if cell is self.prim else self.LD.translate_pairinfo_to_supercell(cell, *self.pairinfo)
        f_phonon.init(cell.lattice._matrix, cell.atomic_masses, cell.frac_coords, *pairinfo)
        if (self.NAC is not None) and (self.NAC.id in [0, 1, 2]):
            self.NAC.uniq_ijk_idx = Union(self.pairinfo[0], True)
            if self.NAC.Ncell < 0:
                self.NAC.Ncell = float(len(self.NAC.uniq_ijk_idx))
            else:
                if abs(self.NAC.Ncell - len(self.NAC.uniq_ijk_idx)) > 1E-6:
                    print("WARNING: input Ncell not equal to the given number of ijk", len(self.NAC.uniq_ijk_idx))
            f_phonon.init_nac(self.NAC.id, cell.intensive_properties['epsilon_inf'],
                cell.site_properties['born_charge'], self.NAC.rho,
                              cell.lattice.volume, self.NAC.uniq_ijk_idx, 1.0)
            if (self.NAC.id==0) and (self.LD._dpcor is not None):
                print("    Loading dipole correction FCM for phonon")
                dpcor= self.LD._dpcor if cell is self.prim else self.LD.translate_pairinfo_to_supercell(cell, *self.LD._dpcor)
                f_phonon.init_dpcor(*dpcor)


    # short hand to get cartesian wave vector
    def to_c(self, kpts, cart):
        return np.array(kpts) if cart else self.prim.reciprocal_lattice.get_cartesian_coords(kpts)


    def get_dm(self, k_in, cart=True):
        """

        :param k_in: list of cartesian K-points
        :return:
        """
# call fortran program
        return f_phonon.get_dm(self.to_c(k_in, cart), self.dim)


    def get_dmnomass(self, k_in, cart=True):
        """
        Get the fourier transformed force constant matrix, without sqrt(M1,M2)
        """
        mass= np.sqrt(self.cell.atomic_masses).repeat(3)
        return self.get_dm(k_in, cart)* (np.outer(mass, mass)[:,:,None])


    @staticmethod
    def _maketicks(dist, lbl, plt):
        """
        private utility method to add ticks to a band structure
        """
        def _get_greek(ls):
            return ["$" + l + "$" if l.startswith("\\") or l.find("_") != -1 else l for l in ls]


        uniq_d = []
        uniq_l = []
        temp_ticks = list(map(list,zip(dist, _get_greek(lbl))))
        for i in range(len(temp_ticks)):
#            if (i < len(temp_ticks)-1) and (temp_ticks[i][1] != temp_ticks[i+1][1]) and abs(temp_ticks[i][0]-temp_ticks[i+1][0])<1E-7:
            if (i < len(temp_ticks)-1) and abs(temp_ticks[i][0]-temp_ticks[i+1][0])<1E-7:
                temp_ticks[i][1] +="|" + temp_ticks[i+1][1] if (temp_ticks[i][1] != temp_ticks[i+1][1]) else ''
                temp_ticks[i+1][1] = "" 
            if i == 0:
                uniq_d.append(temp_ticks[i][0])
                uniq_l.append(temp_ticks[i][1])
#                logger.debug("Adding label {l} at {d}".format(
#                    l=temp_ticks[i][0], d=temp_ticks[i][1]))
            else:
                if temp_ticks[i][1] == temp_ticks[i - 1][1]:
#                    logger.debug("Skipping label {i}".format(
#                        i=temp_ticks[i][1]))
                    continue
                else:
#                    logger.debug("Adding label {l} at {d}".format(
#                        l=temp_ticks[i][0], d=temp_ticks[i][1]))
                    uniq_d.append(temp_ticks[i][0])
                    uniq_l.append(temp_ticks[i][1])

#        logger.debug("Unique labels are %s" % list(zip(uniq_d, uniq_l)))
        plt.gca().set_xticks(uniq_d)
        plt.gca().set_xticklabels(uniq_l)

        for i in range(len(lbl)):
            if lbl[i] is not None and lbl[i]:
                # don't print the same label twice
                if i != 0:
                    if lbl[i] == lbl[i - 1]:
#                        logger.debug("already print label... "
#                                     "skipping label {i}".format(
#                                     i=ticks['label'][i]))
                        continue
                    else:
#                        logger.debug("Adding a line at {d}"
#                                     " for label {l}".format(
#                            d=ticks['distance'][i], l=ticks['label'][i]))
                        plt.axvline(dist[i], color='g')
                else:
#                    logger.debug("Adding a line at {d} for label {l}".format(
#                            d=ticks['distance'][i], l=ticks['label'][i]))
                    plt.axvline(dist[i], color='g')
        plt.axhline(0, color='0.7')
        with open('wavevec_label.txt', 'w') as f: f.write(''.join([x[1]+'\n' for x in temp_ticks]))
        return plt

    #@staticmethod
    def str2kpt(self, s, cart):
        """

        :param s: example "[[10,  [0,0,0],'\\Gamma', [0.5,0.5,0.5], 'X', [0.5,0.5,0], 'K']]"
                  or simply "[[10, 'X', 'K']]"
                  or "Auto N_pt" for automatic K-path with N_pt points on each segment
                  or "Auto" for default density (20)
        :param cart: boolean for cartesian coordinates
        :return:
        """
        from csld.util.mathtool import vec_linspace
        s_list = s.split()
        symkp= HighSymmKpath(self.prim)
        if s_list[0] == 'Auto':
            kpts, lbls = symkp.get_kpoints(20 if len(s_list)<2 else int(s_list[1]))
            kpts = np.array(kpts)
            return kpts, lbls

        d = eval(s)
        lbls= []
        kpts= []
        kpt_dat= symkp.kpath['kpoints']
        lbl_only = isinstance(d[0][1], str)
        npt= d[0][0]
        for x in d:
            if lbl_only:
                print("  found only labels in kpts specs.")
                lb_line = x[1:]
                kp_line = [kpt_dat[k] for k in lb_line]
            else:
                lb_line = x[2::2]
                kp_line = x[1::2]
            for i in range(len(lb_line)-1):
                lbls.extend([lb_line[i]] + (['']*(npt-2)) + [lb_line[i+1]])
            kpts.extend(vec_linspace(kp_line, npt))
#            lbls.append(x[2])
#            for y in x[4::2]:
#                lbls.extend(['']*(x[0]-2))
#                lbls.append(y)
#        kpts = np.vstack([vec_linspace(x[1::2], x[0]) for x in d])
        kpts = np.array(kpts)
        if not lbl_only:
            kpts = self.to_c(kpts, cart)

        return kpts, lbls

    def get_dispersion(self, k_in_s, unit='THz', cart=True):
        """

        :param k_in_s: string "Auto" or K-points definition
        :return:
        """
        #rec = self.prim.lattice.reciprocal_lattice
        kpts, labels = self.str2kpt(k_in_s, cart)
        #np.savetxt('q_frac.txt', self.prim.reciprocal_lattice.get_fractional_coords(kpts))
        xval = np.zeros(len(kpts))
        for i in range(1, len(kpts)):
            xval[i] = xval[i-1] + (0 if labels[i-1] and labels[i] else np.linalg.norm(kpts[i:i+1] - kpts[i-1:i]))
        eigE = f_phonon.get_dispersion(kpts, 1, self.dim)* self.units[unit][0]
        np.savetxt(fn_available('phonon-dispersion.out'), np.hstack((np.array([xval]).T, eigE.T)), header=
          'col1: wavevector distance;  col 2 ...: band 1 ...: '+self.units[unit][1])
        if False:
            ref_eigE= np.loadtxt('ref-dispersion.txt')[:,1:].T
            print(eigE.shape, ref_eigE.shape)
        else:
            ref_eigE= None
        if plt is not None and self.pdfout:
            plt.figure()
            if ref_eigE is not None:
                plt.plot(*my_flatten([[xval, eigE[i,:], 'b-'] for i in range(self.dim)]+
                                     [[xval, ref_eigE[i,:], 'r-'] for i in range(self.dim)]))
            else:
                plt.plot(*my_flatten([[xval, eigE[i,:], 'b-'] for i in range(self.dim)]))
            plt.xlabel("wave vector")
            plt.ylabel(self.units[unit][1])
            #print(xval, labels)
            Phonon._maketicks(xval, labels, plt)
            plt.savefig(self.pdfout, format='pdf')
            plt.close()

        return eigE


    def get_eig_e_vec(self, k_in,unit='THz', cart=True):
        """

        :param k_in: list of K-points
        :return:
        """
        from scipy.io import mmwrite
        eigE, eigV = f_phonon.get_eig_e_vec(self.to_c(k_in, cart), self.dim)
        dm = f_phonon.get_dm(self.to_c(k_in, cart), self.dim)
        np.savetxt("eigE.txt", eigE)
        for i in range(len(k_in)):
            mmwrite("eigV_%d.mtx"%(i), eigV[:,:,i])
            mmwrite("dm_%d.mtx"%(i), dm[:,:,i])
        # np.savetxt("hessian.txt", self.get_FCM())
        return (eigE, eigV)


    def plot_dos(self, plt, x, y, title, unit):
        if plt is not None:
            plt.figure()
            plt.plot(x, y, 'b-')
            plt.xlabel(unit)
            plt.ylabel("Phonon DOS")
            plt.title(title)
            plt.savefig(self.pdfout, format='pdf')
            plt.close()


    def get_dos(self, mesh, nEdos, ismear, epsilon, unit='THz', pdos=False):
        """

        :param mesh: [nx ny nz]
        :param ngrid_en: number of energy points
        :param ismear: -1 for tetrahedron, 0 for Lorentzian smearing, 1 for Gaussian smearing
        :param epsilon: width of smearing
        :param unit:
        :return: numpy array [[t1, dos1], [t2, dos2], ...]
        """
        kgrid = np.mgrid[0:1:1./mesh[0], 0:1:1./mesh[1], 0:1:1./mesh[2]].transpose((1,2,3,0)).reshape(-1,3)
        kred, wt = zip(*self.prim.syminfo.get_ir_reciprocal_mesh(mesh))
        dos, pd = f_phonon.get_dos_new(pdos, mesh, self.to_c(kgrid, False), nEdos, ismear, epsilon, self.prim.num_sites)

        en = dos[:,0]* self.units[unit][0]
        self.plot_dos(plt, en, dos[:,1]/self.units[unit][0], "Total", self.units[unit][1])
        # returned dos always in eV per primitive cell (3 N_atom modes)
        dos[:,0] *= self.units['eV'][0]
        dos[:,1] /= self.units['eV'][0]
        np.savetxt(fn_available('phonon-total-dos.out'), dos, header='col1: energy in eV;  col 2: DOS')

        if pdos:
            for i in range(self.prim.num_sites):
                self.plot_dos(plt, en, pd[i]/self.units[unit][0], "Partial for atom %d"%(i+1), self.units[unit][1])
            pd /= self.units['eV'][0]
            np.savetxt(fn_available('phonon-partial-dos.out'), np.hstack((dos[:,0:1], pd.T))
                       , header='col1: energy in eV;  col 2: atom 1 DOE, etc')

        return dos

    @staticmethod
    def calc_thermal_QHA(dos, Tlist, outf):
        """
        Thermodynamic properties in the quasi-harmonic approximation
        :param dos:
        :param Tlist: list of temperature
        :param outf:
        :return:
        """
        dat = f_phonon.calc_thermal(dos, Tlist)
        np.savetxt(outf, dat, header='QHA Per primitive cell: T (K);  E (eV);  A=E-TS (eV);  S (kB);  Cv (kB)')
        return dat


    def debye_velocity(self, grid, x1rad, x0rad, average=True):
        """
        returns Debye velocity in m/s
        :param grid: [num_dcostheta, num_dphi]
        :param x1rad: fractional WRT the Weigner-Sietz cell of reciprocal space
        """
        from ..util.mathtool import mkgrid_surfaceintegral_spherical
        scale= self.prim.reciprocal_lattice.WS_radius
        samples= mkgrid_surfaceintegral_spherical(grid, True)
        qpt= np.array([[np.sin(p[0])*np.cos(p[1]), np.sin(p[0])*np.sin(p[1]), np.cos(p[0])] for p in samples])*scale
        q1= qpt*x1rad
        q0= qpt*x0rad
        eigE1 = f_phonon.get_dispersion(q1, 1, self.dim)[:3,:].T.reshape((-1))
        eigE0 = np.zeros_like(eigE1) if np.abs(x0rad)<1E-20 else f_phonon.get_dispersion(q0, 1, self.dim)[:3,:].T.reshape((-1))
        # note energy in THz unit (frequency, NOT angular)
        # return unit meter/second
        velocity = 100*2*np.pi*(eigE1-eigE0)*self.units['THz'][0]/((x1rad-x0rad)*scale)
        if np.any(velocity <0):
            print('WARNING: found negative velocity in Debye velocity calculation:\n', velocity)
        #print('db v=', 100/np.power(np.mean(np.power(velocity,-3)),1./3))
        debyeV=np.power(np.mean(np.power(velocity,-3)), -1./3)
        if average:
            # average over all acoustic branches
            return debyeV
        else:
            vmat=velocity.reshape((-1,3)).T
            return np.array([np.power(np.mean(np.power(vmat[i],-3)), -1./3) for i in range(3)])


    def debye_T(self, grid, x1rad, x0rad=0, average=True):
        # the constant here is hbar/k_B/1 angstrom
        const = 0.0763823511
        debye_v = self.debye_velocity(grid, x1rad, x0rad, average)
        debye_t = debye_v*np.power(6*np.pi**2*self.prim.num_sites/self.prim.lattice.volume,1./3)*const
        return debye_t, debye_v

    def get_modes_in_supercell(self, sc, nogamma=False):
        """

        :param sc: supercell
        :return: eigen values and vectors on q-point compatible with the supercell
        """
        qmesh = sc.compatible_kpoints()
        # remove Gamma point
        if nogamma:
            qmesh = np.delete(qmesh, np.where(np.linalg.norm(qmesh, axis=1)<= 1E-8), 0)
        qmesh = self.to_c(qmesh, False)
        eigE, eigV= self.get_eig_e_vec(qmesh)
        # to get hermitian displacement operators, eigV[q]=eigV^*[-q]
        from ..coord_utils import small_fractional
        is_opposite = np.array(np.where(np.linalg.norm(small_fractional(self.prim.reciprocal_lattice.get_fractional_coords(qmesh[:,None,:]+qmesh[None,:,:])),axis=2)<1E-8)).T
        print('debug is oppo', is_opposite)
        return eigE, eigV

    def get_dm_supercell(self, sc):
        Nsc = sc.n_cell
        dim =3*sc.num_sites
        qmesh = sc.compatible_kpoints()
        qmesh = self.to_c(qmesh, False)
        dms= self.get_dm(qmesh)
        lc = self.prim.lattice.get_cartesian_coords(sc.ijk_ref)
        phases= np.transpose(np.exp(-1j*np.dot(lc[None,:,:]-lc[:,None,:], qmesh.T)), (0,2,1))
        return np.transpose(dms.dot(phases)/Nsc,(0,2,1,3)).reshape((dim,dim)).real


    def get_hessian(self, sc, with_setup_cell=True):
        if with_setup_cell:
            self.setup_cell(sc)
            return self.get_dmnomass(np.zeros((1,3)))[:,:,0].real
        else:
        # using the method of getting dm from compatible k-points has problems with indices
        # do not use
            mass= np.sqrt(sc.atomic_masses).repeat(3)
            return self.get_dm_supercell(sc)*np.outer(mass, mass)

    def export_hessian_forshengbte(self, sc):
        from csld.util.tool import matrix2text
        na= self.prim.num_sites
        Nsc = sc.n_cell
        hmat = self.get_hessian(sc, True)
        hmat = hmat.reshape((na,Nsc,3,na,Nsc,3))
        #with open('FORCE_CONSTANTS_2short', 'w') as f:
        #    for ia1 in range(na):
        #        for ia2 in range(na):
        #            for l,ls in enumerate(sc.ijk_ref):
        #                f.write("%d %d %d %d %d %s\n"%(ia1+1,ls[0]+1,ls[1]+1,ls[2]+1,ia2+1,matrix2text(hmat[ia1,0,:,ia2,l,:].reshape(-1))))
        # new compact format
        with open('FORCE_CONSTANTS_2ND', 'w') as f:
            f.write("%d %d\n"%(na, na*Nsc))
            index=np.arange(na*Nsc).reshape((na,*(np.diag(sc.sc_mat)[::-1])))
            for ia1 in range(na):
                for ia2 in range(na):
                    for l,ls in enumerate(sc.ijk_ref):
                        f.write("%d %d\n%s\n"%(index[ia1,0,0,0]+1,index[ia2,ls[2],ls[1],ls[0]]+1,matrix2text(hmat[ia1,0,:,ia2,l,:])))


    def covariance_matrix_in_supercell(self, sc, T):
        from ..util.units import kBbyeV, eVAng_to_THz
        unit_to_ang2=0.505379009
        dim=3*sc.num_sites
#        np.set_printoptions(linewidth=150)
        ###dm_sc = self.get_dm_supercell(sc)
        ### now use get_hessian instead
        self.setup_cell(sc)
        dm_sc = self.get_dm(np.zeros((1,3)))[:,:,0].real
        #for i in range(Nsc):
         #   for j in range(Nsc):
          #      phases= np.exp(-1j*np.dot(qmesh, lc[j]-lc[i]))/Nsc
           #     dm_sc[i*dim0:(i+1)*dim0,j*dim0:(j+1)*dim0] = dms.dot(phases)
#        for i in range(Nsc):
    #    print(qmesh.shape,(lc[:,:,None]-lc.T[None,:,:]).shape, lc[:,:,None]-lc.T[None,:,:])
        
     #   print('phase', phases.shape, phases, 'lc')
 #           dm_sc[dim0*i:dim0*(i+1), :] = dms.dot(phases).reshape((dim0,-1))
      #  dm_sc = np.transpose(dms.dot(phases),(0,2,1,3)).reshape((dim,dim))
#        print(dm_sc)
        #print('test non-hermitian', np.linalg.norm(dm_sc - np.matrix.getH(dm_sc)), np.linalg.norm(dm_sc.real), np.linalg.norm(dm_sc.imag))
        #print('debug dm_sc diag', np.diag(dm_sc.real))
        eigE, eigV = np.linalg.eigh(dm_sc)
        eigE= eigE.reshape((-1))
        eigE= eVAng_to_THz*np.sqrt(np.abs(eigE))*np.sign(eigE)
 #       print('debug shape', eigE.shape, eigV.shape)
#        eigV= eigV.T.reshape((eigV.shape[0], -1))
        if np.any(eigE < -1E-4):
            ValueError('WARNING: imaginary phonon encountered and ignored')
        index_acoustic=np.where(np.abs(eigE)<=1E-4)[0]
        assert len(index_acoustic) == 3, '%d acoustic modes found in covariance_matrix'%(len(index_acoustic))
        eigE[np.where(eigE<1E-4)] = 100*np.max(eigE)
        kBT= kBbyeV*T
        expVal= np.exp(eigE*self.units['eV'][0]/kBT)
        nph=1.0/(expVal-1.0)
#        mass= (np.sqrt(np.array(self.prim.atomic_masses).repeat(3)))[None,:].repeat(Nsc, axis=0).reshape(-1)
        mass= np.sqrt(sc.atomic_masses).repeat(3)
        mass_fac= 1/np.outer(mass, mass)
        covmat=np.zeros((dim,dim))
#        print('debug sc eigE', eigE)
        for imode in range(dim):
            covmat+=(1+2*nph[imode])/eigE[imode]*np.outer(eigV[:,imode], eigV[:,imode])
#            print('sc', imode, eigE[imode], np.linalg.norm(eigV[:,imode]))
        covmat*= unit_to_ang2*mass_fac
        #print(np.sqrt(np.diag(covmat)), covmat.shape)
        return covmat


        eigE, eigV= self.get_modes_in_supercell(sc)
#        print('colume or row', (dms[:,:,0].dot(eigV[:,:,0]).real)/ (eigV[:,:,0].real), eigE[:,0])
        eigE= eigE.reshape((-1))
        # construct supercell eigenvectors from primitive cell ones at commensurate q-points
        eigv_phase = np.exp(1j* np.dot(lc, qmesh.T))
        print('debug phase=', eigv_phase)
        eigV_sc= np.array([ np.outer(eigV[:, iband, iqc], np.exp(1j* np.dot(lc, qmesh[iqc]))/sqrt(Nsc)).reshape(-1)   for iqc in range(Nsc) for iband in range(dim0)]).T
        print('with prim eigV', eigV_sc.shape, np.linalg.norm(eigV_sc.imag))
#        eigV=eigV.real
        if np.any(eigE < -1E-4):
            ValueError('WARNING: imaginary phonon encountered and ignored')
        index_acoustic=np.where(np.abs(eigE)<=1E-4)[0]
        assert len(index_acoustic) == 3, '%d acoustic modes found in covariance_matrix'%(len(index_acoustic))
        eigE[np.where(eigE<1E-4)] = np.Inf
        print('debug prim eigE', eigE)
        kBT= kBbyeV*T
        expVal= np.exp(eigE*self.units['eV'][0]/kBT)
        nph=1.0/(expVal-1.0)
        mass= (np.sqrt(np.array(self.prim.atomic_masses).repeat(3)))[None,:].repeat(Nsc, axis=0).reshape(-1)
        mass_fac= 1/np.outer(mass, mass)
#        covmat=np.zeros((dim,dim),dtype=np.complex_)
#        for imode in range(dim):
#            covmat+=(1+2*nph[imode])/eigE[imode]*np.outer(eigV_sc[:,imode], np.conj(eigV_sc[:,imode]))
        print('debug np.diag((1+2*nph)/eigE)', (1+2*nph)/eigE)
        eigV_sc=eigV_sc*np.sqrt((1+2*nph)/eigE)
        covmat = np.dot(eigV_sc, np.matrix.getH(eigV_sc))
        covmat*= unit_to_ang2*mass_fac
        print('debug imag of covma', np.linalg.norm(covmat.imag))
        print(np.sqrt(np.diag(covmat.real)), covmat.shape)

        print('now literally')
        covmat=np.zeros((dim,dim),dtype=np.complex_)
        for iqc in range(Nsc):
            for iband in range(dim0):
                imode= iqc*dim0+iband
                covmat+=(1+2*nph[imode])/eigE[imode]*np.outer(eigV[:,iband,iqc], np.conj(eigV[:,iband,iqc]))
        covmat*= unit_to_ang2*mass_fac
        print('debug imag of covma', np.linalg.norm(covmat.imag))
        print(np.sqrt(np.diag(covmat.real)), covmat.shape)
        return covmat

    def supercell_snapshot(self, sc, T, Nfile, classical=False, path='.'):
        """

        :param sc: supercell
        :param T: temperature
        :param classical: whether phonon MSD is classical
        :return: Supercell with displacement according to phonon spectrum thermalized at T
        """
        import scipy
        assert classical==False, 'classical not implementeded yet'
        covmat= self.covariance_matrix_in_supercell(sc, T)
        Lcovmat = scipy.linalg.cholesky(covmat,lower=True)
        ndim=len(covmat)
        na_sc=sc.num_sites
        Poscar(sc).write_file('POSCAR_snapshot_nodisplacement')
        for i in range(Nfile):
            disp=np.random.normal(size=ndim)
            disp=np.dot(Lcovmat, disp).reshape((na_sc, 3))
            Poscar(sc.from_displacement(disp)).write_file('POSCAR_snapshot_%05d'%(i))

    def export_phononmode(self, sc, amplist, modes, posRef=None, cart=True):
        """

        :param sc: supercell
        :param amplist: list of amplitudes (may be negative)
        :param modes: list of [kp1, kp2, kp3, band_index]
        :return:
        """
        from cmath import exp
        qmesh = sc.compatible_kpoints()
        # # remove Gamma point
        # qmesh = np.delete(qmesh, np.where(np.linalg.norm(qmesh, axis=1)<= 1E-8), 0)
        q_c = self.to_c(qmesh, False)
        print('    supercell qmesh', qmesh)
#        print('mass', self.prim.atomic_masses)
        na = self.prim.num_sites
        rec = self.prim.lattice.reciprocal_lattice
        scref_c = self.prim.lattice.get_cartesian_coords(sc.sc_ref)
        cor0 = np.array(sc.cart_coords)
        out_ord = sc.get_order_wrt(posRef, True)

        #print(sc, amplist, modes)
        for i in range(len(modes)):
            kpt = self.to_c([modes[i, :3]], cart)
            iband = int(modes[i, 3])
            qdist = [rec.get_all_distances(q_c[i:i+1], kpt)[0,0] for i in range(len(qmesh))]
            iq = np.argmin(qdist)
            if qdist[iq] > 1E-5:
                print('*** WARNING MATCHING Q-point NOT found WARNING ***')
            q = q_c[iq:iq+1]
            eigE, eigV = self.get_eig_e_vec(q)
        #    print('debug EIGEN EN=', eigE, "iband=", iband, eigE.shape,eigV.shape)
            eig_en= eigE[iband, 0]
            eig_vec=eigV[:, iband, 0]
            print(' selecting', qmesh[iq],  'for mode %d '%(i), modes[i, :3], 'eigE=', eig_en)
            tmp_dm = self.get_dm(q)
            np.savetxt('dm_%d.txt'%(i), tmp_dm)
            u = np.zeros((na, sc.n_cell, 3), dtype=complex)
            for ia in range(na):
                for icell in range(sc.n_cell):
                    u[ia, icell, :]= sqrt(1/(2*sc.n_cell*self.prim.atomic_masses[ia]*abs(eig_en)))* \
                        exp(complex(0., np.dot(q, scref_c[icell])))* eig_vec[ia*3:3*ia+3]

            maxu = u.flatten()[np.argmax(np.abs(u))]
            du = np.reshape(np.real(u/maxu), (-1, 3))
            du /= np.max(np.linalg.norm(du, axis=1))
#            print(du)
            poscar = Poscar(sc)
            for aop in amplist.tolist():
                #print(sc)
                sc.set_coords((cor0 + aop* du)[out_ord], cart=True)
                #print(sc)
#                print(sc.cart_coords)
                poscar.write_file("mode%d_disp%s_band%d.POSCAR" % (i, aop, iband), direct=True)

#        from ..structure import SupercellStructure
#        SCinfo = SupercellStructure(self.prim, np.eye(3,dtype=int), None, self.prim.frac_coords)
#        print("  calc prim FCM") 
#        fcm_dp = self.LD.get_fcm_dipole(SCinfo)
#        np.savetxt("prim_dp", fcm_dp)


    def mean_square_displacement(self, T, outf):
        """

        :param T:
        :param outf:
        :return:
        """
