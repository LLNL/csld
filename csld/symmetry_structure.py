#!/usr/bin/env python3

"""
This module implements symmetry-related structure forms.
"""


import numpy as np
from scipy.linalg import block_diag
import scipy.sparse

from .structure import Structure
from .util.mathtool import tensor_constraint, get_symmetrized_lsq, mychop
from _c_util import fct_trans_c


class SymmetrizedStructure(Structure):
    """
    This class represents a symmetrized structure, i.e. a structure
    where the spacegroup and symmetry operations are defined. This class is
    typically not called but instead is typically obtained by calling
    pymatgen.symmetry.SymmetryFinder.get_symmetrized_structure.

    Args:
        structure (Structure): Original structure
        spacegroup (Spacegroup): An input spacegroup from SymmetryFinder.
        equivalent_positions: Equivalent positions from SymmetryFinder.

        NOTE: conceptually we define **** orbits = sublattices ****
    .. attribute: equivalent_indices

        indices of structure grouped by equivalency
    """

    def __init__(self, structure, spacegroup, equivalent_positions, syminfo, tol=1E-3, warn_prim=2):
        Structure.__init__(self, structure.lattice,
                           [site.species_and_occu
                            for site in structure],
                           structure.frac_coords,
                           site_properties=structure.site_properties)

        self.spacegroup = spacegroup
        self.syminfo=syminfo
        self._pointgroup = None
        self.spacegroup.check_primitive(warn_prim)
        self.nsymm = len(spacegroup)
        u, inv = np.unique(equivalent_positions, return_inverse = True)
        self.l_of_orbit = u
        self.orbit_of_l = inv
        self.reprL_of_l = equivalent_positions
        self._l_list_of_orb = None # wait till setup of orbits
        self.wyckoffs = syminfo.get_symmetry_dataset()['wyckoffs']
        # print("DEBUG equiv", equivalent_positions,'spg=',spacegroup,spacegroup.__class__)
        # print("debug u=", u)
        # print("debug inv=", inv)
        # self.equivalent_indices = [[] for i in range(len(u))]
        # for i, inv in enumerate(inv):
        #     self.equivalent_indices[inv].append(i)
        # self.orb_id = equivalent_positions
        self.n_orb = len(u)
        # print("debug equivalent_indices=", self.equivalent_indices)
        # print("debug> coords ", self.cart_coords, self.frac_coords)
        # print("debug> orb_id=", equivalent_positions)

        #print("DEBUG symmops=",[x for x in spacegroup])
#
#  define mapping of point l by group g
#   x = ijk, l
# g x = ijk', l' = g( ijk + t_l) = Ag ijk + Ag t_l + t_g = Ag ijk + _ijk_tbl(l, g) + t_(_l_tbl)(l, g)
# so
# ijk'= Ag ijk + _ijk_tbl(l,g)
# l' = _l_tbl(l, g)
# The look-up tables _ijk_tbl and _l_tbl can be constructed by looking at the primitive cell only
# Also need the integer matrices Ag, i.e. rotation matrix in fractional coordinates
#
        ijkl_tbl = np.zeros((structure.num_sites, self.nsymm, 4), dtype=np.int)
        # l_tbl = np.zeros((structure.num_sites, self.nsymm), dtype=np.int)

        for il in range(structure.num_sites):
#            print("debug> site ", il, self.sites[isite])
            for ig in range(self.nsymm):
#                mapped_site= spacegroup.symmops[jsym].operate(self.frac_coords[isite])
                #print("debug", isite, jsym, mapped_site, self.identify_site(mapped_site))
#        self.verify_symmetry()
                x = self.frac_coords[il]
                gx = self.spacegroup[ig].operate(x, True)
#                print(il, ig, self.spacegroup.symmops[ig].operate(x, True), self.spacegroup.symmops[ig].operate(x, False))
#                try:
                ijklp = self.frac2ijkl(gx, True, tol)
#                except ValueError:
#                    print("position",x, "mapped to", gx, " NOT in the cell. Incompatible structure & symmetry!!")
#                    exit(-1)
                #print("DEBUG ",il, ig, ijklp)
                ijkl_tbl[il, ig] = ijklp
                # l_tbl[il, ig] = lp
        self._ijkl_tbl = ijkl_tbl
        # self._l_tbl = l_tbl
        self._Arot_int = [block_diag(g.rot_f.astype(np.int), [[0]]) for g in self.spacegroup]
        #self.alloy_species=None
        self._orbits = None


    def _setup_orbit(self):
        from .cluster import Cluster
        self._orbits = Cluster.clusters_to_orbit_isotropy(self, [Cluster.from_ijkl([[0, 0, 0, l]], self) for l in self.l_of_orbit])
        self._l_list_of_orb = [[c._l_list[0] for c in o.clusters] for o in self._orbits]
        assert self.n_orb == len(self._orbits), \
            ValueError('ERROR: %d orbits from spglib but %d from orbits'%(self.n_orb, len(self._orbits)))
        elements = self.types_of_elements
        id = np.zeros(len(elements), dtype=int)
        name_orb=[]
        for io, ol in enumerate(self._l_list_of_orb):
            idx_ele = elements.index(self.elements[ol[0]])
            id[idx_ele] += 1
            name_orb.append("%s%d"%(elements[idx_ele], id[idx_ele]))
        self.name_orb = name_orb
        # print('debug lattice orbits=', self._orbits, '_l_list_of_orb', self._l_list_of_orb)
        # print('debug spg rots=', self.syminfo.get_symmetry_dataset()["rotations"])
        # print('debug spg all=', self.syminfo.get_symmetry_dataset())
        # print('debug repr=', self.__repr__())


    @property
    def equivalent_indices(self):
        return self.l_list_of_orb


    @property
    def l_list_of_orb(self):
        """
        Finds all symmetrically equivalent indices for a particular index

        l: integer between 0 and num_sites-1
        Returns: List of all symmetrically equivalent indices.
        """
        if self._l_list_of_orb is None:
            self._setup_orbit()
        return self._l_list_of_orb

    @property
    def orbits(self):
        if self._orbits is None:
            self._setup_orbit()
        return self._orbits

    @property
    def pointgroup(self):
        if self._pointgroup is None:
            self._pointgroup = 0
        return self._pointgroup


    def symmetrize_coordinates(self,tol=1E-3):
        """
        OBSOLETE. Use symmetrize() instead
        coordinates of each atom become average of its isotropy group mappings
        :return:
        """
        s_pos = np.zeros((self.num_sites, 3))
        for ia in range(self.num_sites):
            x = self.frac_coords[ia]
            for ig in range(self.nsymm):
                gx = self.spacegroup[ig].operate(x, True)
                ijklp = self.frac2ijkl(gx, True, tol)
                ib = ijklp[-1]
                s_pos[ib] += gx + np.around(self.frac_coords[ib] - gx)
        for ia in range(self.num_sites):
            dist=self.sites[ia].distance(s_pos[ia]/self.nsymm)
            if dist>1E-12:
                print("WARNING: symmetrization moved site %4d by %6g"%(ia+1, dist))
            self.sites[ia].set_coords(s_pos[ia]/self.nsymm, cart=False)

    def symmetrize(self, symprec=1e-12):
        """
        using Structure.standardize_cell()
        :return the symmetrized primitive cell
        """
        import scipy
        prim_norefine= self.standardize_cell(True, True, symprec)
        if prim_norefine.num_sites != self.num_sites:
            raise ValueError("ERROR this is NOT a primitive cell")
        prim_refine= self.standardize_cell(True, False, symprec)
        basis_rot = prim_norefine.lattice.inv_matrix.dot(prim_refine.lattice.matrix)
        if np.linalg.norm(basis_rot-np.eye(3))<1e-1:
            newlat= prim_refine.lattice.matrix
        else:
            q,r = scipy.linalg.qr(basis_rot)
            # print('debug q,r',q,r, 'rot', basis_rot)
            newlat = prim_refine.lattice.matrix.dot(q.T)
        newprim= Structure(newlat, prim_refine.species_and_occu,prim_refine.frac_coords)
        # print('original', self.lattice.matrix)
        # print('F, F',self.standardize_cell(False, False, symprec).lattice.matrix)
        # print('F, T',self.standardize_cell(False, True, symprec).lattice.matrix)
        # print('T, F',self.standardize_cell(True, False, symprec).lattice.matrix)
        # print('T, T',self.standardize_cell(True, True, symprec).lattice.matrix)
        scmat= newprim.lattice.matrix.dot(self.lattice.inv_matrix).round().astype(np.int)
        if not (scmat==np.eye(3,dtype=int)).all():
            print("WARNING Symmetrizing primitive cell: lattice shape changed by", scmat, newprim.lattice.matrix.dot(self.lattice.inv_matrix))
        else:
            dlat=newprim.lattice.matrix-self.lattice.matrix
            if np.linalg.norm(dlat)>1e-12:
                print("WARNING Symmetrizing primitive cell: change in lattice\n", dlat)
        # print('debug cart coords=', np.hstack((newprim.cart_coords,self.cart_coords,newprim.cart_coords-self.cart_coords)))
        # np.savetxt('tmp.txt', np.hstack((newprim.cart_coords,self.cart_coords,newprim.cart_coords-self.cart_coords)),fmt='%8g')
        # np.savetxt('tmpf.txt', np.hstack((self.lattice.get_fractional_coords(newprim.cart_coords),self.frac_coords)),fmt='%8g')
        # np.savetxt('tmpjustf.txt', np.hstack((newprim.frac_coords,self.frac_coords)),fmt='%8g')
        dfrac = self.lattice.get_fractional_coords(newprim.cart_coords-self.cart_coords)
        dfrac-=np.round(dfrac)
        dist= np.linalg.norm(dfrac.dot(self.lattice.matrix), axis=1)
        for ia, dis in enumerate(dist):
            if dis>1E-12:
                print("WARNING Symmetrizing primitive cell: moved site %4d by %6g"%(ia+1, dis))
        return newprim


    def symmetrize_born(self, bornlist):
        ## ASR for effective charges, i.e. charge neutrality
        asr_err = np.sum(bornlist, axis=0)
        if np.max(np.abs(asr_err)) > 1e-5:
            print('WARNING Born charge ASR sum=\n', asr_err)
        bsym = self.symmetrize_site_tensors(bornlist)
        # print('debug after sym asr=', np.sum(bsym, axis=0))
        bsym -= np.sum(bsym, axis=0)/self.num_sites
        return mychop(bsym)


    def symmetrize_site_tensors(self, tlist):
        """
        symmetrize tensors on atoms, one orbit at a time
        :param tlist:
        :return:
        """
        tensors=np.array(tlist)
        rank = np.array(tensors[0]).ndim
        dim = tensors[0].shape[0]
        tensors_s = np.zeros_like(tensors)
        for i,orb in enumerate(self.orbits):
            cmat= tensor_constraint(3, rank, [self.spacegroup[ig].rot for ig,pi in orb.isotropy if ig>0])
            # print('debug sublatice', i , orb.pointgroup_symbol[0], self.name_orb[i], 'shape', cmat.shape)
            # print('debug sublat', i, 'cmat', cmat, 'orbit', orb, 'isotropy', orb.isotropy)
            tr_mat = scipy.sparse.vstack([fct_trans_c(rank, dim, self.spacegroup[orb.clusters_ig[ic]].rot,
                      np.arange(rank, dtype=int)).dot(cmat.T) for ic in range(self.orbits[i].multiplicity)])
            t_s = get_symmetrized_lsq(tr_mat, tensors[self.l_list_of_orb[i]].reshape((-1)))
            tensors_s[self.l_list_of_orb[i]] = t_s.reshape(tensors[self.l_list_of_orb[i]].shape)
        return tensors_s


    def symmetrize_tensor(self, tensor):
        rank = np.array(tensor).ndim
        pgrots = [op.rot for op in self.syminfo.get_pointgroup()]
        cmat= tensor_constraint(3, rank, pgrots)
        t_s = get_symmetrized_lsq(cmat.T, tensor.reshape((-1)), 1e-5)
        return t_s.reshape(tensor.shape)


    def operate_ijkl(self, ig, ijkl):
        """

        :param ijkl:
        :return: New [ijk', l']
        """
        return np.dot(self._Arot_int[ig], ijkl) + self._ijkl_tbl[ijkl[3], ig]


    def verify_symmetry(self, tol=1E-5):
        """
        check to see if every atom is mapped back into the structure
        :param tol:
        :return:
        """
        print("WARNING: should be obsolete now. No longer needed. Superseded by the symmetry look-up table")
        for i in range(self.nsymm):
            g= self.spacegroup.symmops[i]
            for x in self.frac_coords:
                gx = g.operate(x, True)
                try:
                    mapped = self.frac2ijkl(gx, True, tol)
                except ValueError:
                    print("position",x, "mapped to", gx, " NOT in the cell. Incompatible structure & symmetry!!")
                    return False
        # if debug_level>10:
        #     print("  symmetry verified")
        return True

    @staticmethod
    def init_structure(setting, step, symm_prim, debug_level=0, in_cell=True, read_CE=True, check_prim=2):
        """

        :param setting:
        :param step:
        :param symm_prim: symmetrize the cell (to high precision)
        :param debug_leve:
        :param in_cell: move inside cell if True, otherwise report error
        :return:
        """

        # init primitive cell
        from csld.analyzer import SpacegroupAnalyzer
        from csld.interface_vasp import Poscar
        from csld.util.string_utils import setting2arr

        pos = Poscar.from_file(setting['prim'], read_CE=read_CE)
        struc = pos.structure
        symprec = float(setting.get('sym_tol', '1e-5'))
        for isite, site in enumerate(struc.sites):
            if any(np.mod(site._fcoords, 1) != site._fcoords):
                if in_cell:
                    site._fcoords = np.mod(site._fcoords, 1)
                else:
                    print("  ERROR site %d must be inside unit cell. Fractional coords= %.4f %.4f %.4f"%
                          (isite+1, site._fcoords[0], site._fcoords[1], site._fcoords[2]))
                    exit(-1)
        ## init space group of primitive cell
        if step == 0:
            exit(0)
        elif step == 1:
            ## read in sym.out
    #        sgop = SpaceGroup.from_file(settings['structure']['spacegroup_in'])
    #        prim = SymmetrizedStructure(pos, sgop)
            print("  reading sym.out is NOT supported. Please use --symm_step 2 or 3")
            exit(0)
        elif step in [2,3]:
            syminfo = SpacegroupAnalyzer(struc, symprec)
            prim = syminfo.get_symmetrized_structure(check_prim=check_prim)
            sgop = prim.spacegroup
            if step == 3:
                sgop.write(setting['spacegroup_out'])
                print("    Space group operators written to "+setting['spacegroup_out'])
        else:
            print("ERROR: Unknown symm_step: ", step)
            exit(-1)

        ## other properties
        prim.intensive_properties['epsilon_inf'] = setting2arr(setting, 'epsilon_inf', (3,3), 'epsilon_inf.txt')
        prim.add_site_property('born_charge', setting2arr(setting, 'born_charge', (-1,3,3), 'born_charge.txt'))
        print("+ Space group done. Found ", sgop, "%d symmops"%(len(sgop)))
        if symm_prim:
            #prim.symmetrize(symprec)
            prim.symmetrize_coordinates(symprec)
            pos.structure = prim
            print(" ++ primitive cell coordinates symmetrized")
            if setting.getboolean('writeprim',True):
                with open("POSCAR_symmetrized","w") as f:
                    pos.write_file(f)
        print(" ++ Found %d sub_lattices"%(prim.n_orb))
        if debug_level > 0:
            print('*'*10, setting['prim'], '*'*10)
            print(repr(prim))
            print('*'*10, 'END', '*'*10)
        return prim


    def __repr__(self):
        outs = ["Structure Summary (%s)"%(self.spacegroup), repr(self.lattice)]
        for i, o in enumerate(self.l_list_of_orb):
            outs.append("%4d (%2d%s %5s) %6s %s"%(i+1, len(o), self.wyckoffs[o[0]], self.orbits[i].pointgroup_symbol[0], self.name_orb[i], str(self.sites[o[0]])))
        return "\n".join(outs)

