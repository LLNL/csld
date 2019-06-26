#!/usr/bin/env python3

"""
This module provides classes used to
  * define a phonon mode
  * apply symmetry op to a phonon mode
"""


import numpy as np
from numpy.linalg import norm

from .interface_vasp import Poscar
from .cluster import Cluster
from .coord_utils import ReadPBC2Cart, match_p1p0
from .util.mathtool import relativePosition
import logging
logger = logging.getLogger(__name__)

debug_level = 10

class PhononMode():
    """

    """

    def __init__(self, prim, clus, u, clus_dominant):
        """
        a proper cluster decorated with displacement u
        """
        self.prim = prim
        self.cluster = clus
        self.u = u
        self.dominant_cluster = clus_dominant
        # print(self.cluster)
        # print(self.dominant_cluster)

    def operated_by(self, isym):
        return PhononMode(self.prim, Cluster.from_ijkl(self.cluster.operated_by(isym), self.prim),
                          np.dot(self.u, self.prim.spacegroup[isym].rot.T),
                          Cluster.from_ijkl(self.dominant_cluster.operated_by(isym), self.prim))


    def __repr__(self):
#        outs = [""]
#        outs.append(self.frac_coords.__str__())
#        outs.append("ijk= "+ self.ijk.__repr__()+ " l="+ self.l.__repr__()+"]")
#        outs.append(str(self.ijkl))
#        return " ".join(outs)
        return str(self.cluster)+ repr(self.u)

    @classmethod
    def from_file(cls, prim, p1f, disp_cut, natom_dominant=2, p0f=None, dominant_specie=None, dominant_dr=True):
        """

        :param p0: read an ideal structure,
         p1: disturbed structure
         disp_cut: retain only atoms with displacement larger than this cutoff
         dominant_specie: which atoms may be the representative
        :paramdominant_dr: Select dominant pairs by largest change in distance
        :return: a phonon mode
        """
        if p0f is None:
            x0 = Poscar.from_file(p1f).structure.map_to_prim(prim)
        else:
            assert isinstance(p0f, str), TypeError("require p0f filename")
            x0 = Poscar.from_file(p0f).structure
        scmat = prim.get_scmat(x0.lattice.matrix)

        if debug_level > 10:
            print("supercell clusters generated")

        u = ReadPBC2Cart(p1f, x0.frac_coords)
        # center of mass
        u -= np.average(u, axis=0)
        unorm = norm(u, axis=1)
        # print(unorm)
        # print(u)
        uidx = np.where(unorm > disp_cut)[0]
        # print(uidx)
        if len(uidx) < natom_dominant:
            print("ERROR: found only %d displaced atoms" % (len(uidx)))
            exit(-1)

        clus_xyz = x0.frac_coords[uidx].copy()
        # move cluster atoms close to one atom considering PBC
        imax = np.argmax(unorm)
        pc0 = x0.frac_coords[imax:imax+1].copy()
        for i in range(len(uidx)):
            # print("debug", i, ui, clus_xyz[i], clus_xyz[i:i+1], pc0)
            clus_xyz[i] = match_p1p0(clus_xyz[i:i+1], pc0)[0]
            # print(clus_xyz[i], x0.frac_coords[ui])
        clus_xyz = np.dot(clus_xyz, scmat)

        # now find out the dominant cluster
        if dominant_dr and natom_dominant == 2:
            max_dr = -1.
            sel = [True for _ in uidx]
            if isinstance(dominant_specie, list):
                ele = x0.elements
                sel = [ele[i] in dominant_specie for i in uidx]
#            print(sel, ele)
            for i1, i in enumerate(uidx):
                if not sel[i1]:
                    continue
                for i2, j in enumerate(uidx):
                    if (not sel[i2]) or (i1 >= i2):
                        continue
                    r0 = prim.lattice.norm(clus_xyz[i2]-clus_xyz[i1])
                    r1 = np.linalg.norm(prim.lattice.get_cartesian_coords(clus_xyz[i2]-clus_xyz[i1]) + u[j]-u[i])
                    dr = abs(r0-r1)
                    if dr > max_dr:
                        max_dr = dr
#                        print(i, j, max_dr)
                        atom_dominant = np.array([i1, i2])
        else:
            if dominant_specie is None:
                udominant = unorm[uidx]
            else:
                ele = x0.elements
                assert isinstance(dominant_specie, list)
                # print("debug species=", dominant_specie, ele, "for finding dimer |u|=", [unorm[i] if ele[i] in dominant_specie else -1. for i in uidx])
                udominant = [unorm[i] if ele[i] in dominant_specie else -1. for i in uidx]
            atom_dominant = np.argsort(udominant)[-natom_dominant:]
#        print('dominant=', atom_dominant)
        for i in uidx[atom_dominant]:
            print("atom %4d element %3s   disp %8.3f" % (i, x0[i].specie, unorm[i]), x0.frac_coords[i])


        return cls(prim, Cluster.from_coords(clus_xyz, prim), u[uidx],
                   Cluster.from_coords(clus_xyz[atom_dominant], prim))


    def apply_to(self, dest, dest_clus, inv_sc_mat):
        """
        Apply phonon mode to a supercell
        :param dest: destination supercell
        :param dest_clus: dominant cluster by which to identify or match the phonon
        :param inv_sc_mat: inverse sc_mat belonging to the supercell
        :return:
        """
        prim = self.prim
        syms = prim.spacegroup
        for isym, sym in enumerate(syms):
            s_pl = self.operated_by(isym)
            is_eq, tr = dest_clus.equivalent_by_lattranslation(s_pl.dominant_cluster.ijkls)
            if is_eq:
                # print(isym, tr, s_pl.dominant_cluster.ijkls, dest_dimer.ijkls)
                # found a sym op to move polaron to desired position
                break
        assert is_eq, ValueError("cannot apply phonon mode")
        pi = relativePosition(dest_clus.ijkls, tr)
        tr_ijk = np.array(dest_clus.ijkls[0]) - np.array(s_pl.dominant_cluster.ijkls[pi[0]])
        s_pl.cluster.move_to_prim(dijk=tr_ijk[:3])
        dest_polaron = [dest.frac2ijkl(np.dot(s_pl.cluster.frac_coords[i], inv_sc_mat))[3] for i in
                        range(s_pl.cluster.order)]

        out_sc = dest.copy()
        for i, isite in enumerate(dest_polaron):
            out_sc[isite].move_by(s_pl.u[i])
        return out_sc
