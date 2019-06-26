#!/usr/bin/env python3

"""
This module provides the basic lattice model as used in e.g. cluster expansion, lattice dynamics
"""


import numpy as np
import logging
import scipy.sparse
from collections import Counter

from .cluster import Cluster
from .atomic_model import AtomicModel
from .util.mathtool import relativePosition
from .util.tool import pad_right, matrix2text
from _c_util import get_structure_ordering
logger = logging.getLogger(__name__)


class BasicLatticeModel(AtomicModel):
    """
        A generic model with variables on the lattice:
    """
    def __init__(self, prim, raw_clusters, irange=None, **kwargs):
        """

        :param prim:
        :param raw_clusters:
        :return:
        """
        self.clusters = raw_clusters
        self.prim = prim
        self.maxorder=kwargs['maxorder']
        self.proper_range=kwargs['prange']
        self.cluster_filter = kwargs['filter']
        self.orbits = None
        self.imp_range = irange


    @classmethod
    def from_file(cls, prim, clusf, **kwargs):
        raw = Cluster.read_cluster(prim, clusf)
        return cls(prim, raw, **kwargs)


    def save_clusters(self, fname):
        self.save_clusters_to_file(self.clusters, fname)

    def save_clusters_to_file(self, clusters, fname):
        scale = self.prim.lattice._scale
        elements= self.prim.elements
        sublat = self.prim.orbit_of_l
        out_s = []
        out_ijkl = []
        for cl in clusters:
            out_s.append("1\n%.8f\n%s\n%s\n" % (cl.diameter/scale, " ".join([str(i) for i in cl.orders] +
                                ([str(i[0]) for i in cl.id] if hasattr(cl, 'id') else [])), matrix2text(pad_right(
                np.array(cl.frac_coords, dtype=object), (cl.order, 5)))))
            out_ijkl.append(", ".join(["%d %d %d %d %s%d"%tuple(i.tolist()+[elements[i[-1]],sublat[i[-1]]]) for i in cl._ijkls_np]))
        open(fname, "w").write("\n".join(out_s))
        open(fname+'_ijkl', "w").write("\n".join(out_ijkl))


    @classmethod
    def generate_clusters(cls, prim, **kwargs):
        """
        Generating PROPER clusters.
        Iteratively generation starting from empty and point clusters
        :param prim: Primitive cell
        :return: LDModel with certain clusters
        """
        sites = kwargs['sites'] if 'sites' in kwargs.keys() else list(range(prim.num_sites))
        proper_range = kwargs['prange']
        maxord = max([ord for ord in proper_range.keys() if proper_range[ord]>0])
        raw= []
        clus_per_ord = []
        for ord in range(maxord+1):
            clus_per_ord = cls.generate_n(prim, ord, proper_range[ord], clus_per_ord, sites)
            raw.extend(clus_per_ord)
        return cls(prim, raw, **kwargs)


    def generate_improper(self):
        """

        :return: Improper clusters
        """
        improp = []
        assert self.imp_range is not None, ValueError('improper range not specified')
        for ord, cut in self.imp_range.items():
            if ord <= 1:
                continue
            for clus in self.clusters:
                if (clus.diameter <= cut) and (clus.order< ord) and clus.order >0:
                    # print("extending", ord, cut, clus)
                    improp.extend(Cluster.improper_from_proper(clus, ord, self.prim))
        print("  generated %d improper clusters" %(len(improp)), self.tally_clusters(improp, 2, self.maxorder))
        self.clusters.extend(improp)


    def tally(self, minord=0):
        return self.tally_clusters(self.clusters, minord, self.maxorder+1)

    @staticmethod
    def tally_clusters(clusters, minord=0, maxord=0):
        return dict(Counter([c.order for c in clusters]))


    @staticmethod
    def generate_n(prim, n, cutoff, clus0, sites):
        """
        Input:
        prim cell
        order n
        cutoff for diameter of clusters
        clus0 clusters at order n-1
        :param sites: allowed sites in the primitive cell. Default all sites
        """
        pts = [[0, 0, 0, l] for l in sites]
        if n == 0:
            uc = [Cluster([], prim)]
        elif n == 1:
            uc = [Cluster.from_ijkl([pt], prim) for pt in pts]
            uc = Cluster.remove_equivalent_clusters(uc)
            #assert len(uc) == prim.n_orb
        elif n >= 2:
            print("  generating order %d clusters..."% (n))
            uc = []
            for clus in clus0:
                if clus.diameter > cutoff:
                    continue
                ijkl0 = clus.ijkls
                # each cluster is now sorted by the orbit (sub-lattice) type
                max_id = prim.orbit_of_l[ijkl0[-1][-1]]
                sumPts = prim.find_nb_cluster(np.array(clus.ijkls), cutoff)
                clus_new = []
                for pt_n in sumPts:
                    pt = pt_n.tolist()
                    if not pt[-1] in sites:
                        continue
                    if prim.orbit_of_l[pt[-1]] < max_id:
                        continue
                    if pt in ijkl0:
                        continue
                    clusSum = clus.append_site(pt)
                    if clusSum.diameter > cutoff:
                        continue
                    clus_new.append(clusSum)
                clus_new = Cluster.remove_equivalent_clusters(clus_new)
                uc.extend(clus_new)
            uc = Cluster.remove_equivalent_clusters(uc)
        return uc





    def __str__(self):
        if self.orbits is None:
            outs= ["  No.  ord uniq  diam"]
            for i, c in enumerate(self.clusters):
                outs.append("%4d %4d %4d %9.5f" % (i, c.order, c.order_uniq, c.diameter))
            return "\n".join(outs)
        else:
            outs= ["  No.  ord  uniq mult  diam"]
            for i, orb in enumerate(self.orbits):
                c = orb.cluster
                outs.append("%5d %4d %4d %4d %9.5f" % (i, c.order, c.order_uniq, len(orb.clusters), c.diameter))
            return "\n".join(outs)


    def cleanup(self):
        """
        1. filter through clusters
        2. move first atoms of each cluster to primitive cell
        3. sort clusters
        :return:
        """
        self.clusters= list(filter(self.cluster_filter, self.clusters))
        for i in range(len(self.clusters)):
            self.clusters[i].move_to_prim()
        self.clusters.sort()


    def l_point_cls(self):
        return [orb.cluster.ijkls[0][-1] for orb in self.orbits if orb.cluster.order==1]



    def symmetrize(self):
        self.Cmat = scipy.sparse.identity(self.ncorr_full)


    @property
    def ncorr(self): return self.Cmat.shape[0]


    def prepare_index_full(self):
    # implementations should set the number of correlations for each cluster. 
    # by default set to 1 per cluster
        for orb in self.orbits:
            if orb.ncorr_full<0:
                print("WARNING: ncorr_full NOT set for cluster and set to 1:",orb.cluster)
                orb.ncorr_full = 1
        self.orb_idx_full=np.cumsum([0]+[orb.ncorr_full for orb in self.orbits])
        self.ncorr_full= self.orb_idx_full[-1]


    def prepare_index(self):
        npara_ord = [0]*(self.maxorder+1)
        for orb in self.orbits:
            npara_ord[orb.cluster.order]+= orb.ncorr_full
        self.ord_idx=np.cumsum([0]+npara_ord)
        self.orb_idx=self.orb_idx_full


    def get_orbit_isotropy(self):
        """
        1. Find orbit of each cluster, modulus lattice translation
        2. Find isotropy group of each representative (first) cluster of an orbit
        :return: list of [orbit, isotropy]
            orbit is a list of Clusters
            isotropy is list of [ig, relativePos]
                ig = symmop index
                relativePos indexing function under this symmop
        """
        if self.orbits is None:
            self.orbits = Cluster.clusters_to_orbit_isotropy(self.prim, self.clusters)
        return self.orbits


    def translate_to_supercell(self, sc, orb_idx=None, clus_idx=None):
        """

        :param sc: SupercellStructure
        :return: all clusters in the supercell. Each cluster is give as
                [list of the indices of the vertices, i.d. of orbit, i.d. of symop]
        orb_idx: default (None) to process all; explicit e.g. [0, 3, 5] to process selected few
        clus_idx: default (None) to process all; explicit e.g. [0] to process selected few
        """
        # natom = self.prim.num_sites

        # def match_ordering(l1, l2):
        #     """
        #     ordering such that l1[ordering] == l2
        #     """
        #     n= len(l1)
        #     ordering = list(range(n))


        orb_range = orb_idx if orb_idx is not None else range(len(self.orbits))
        allclus = []
        for iorb in orb_range:
            orbit = self.orbits[iorb]
            clus_range = clus_idx if clus_idx is not None else range(len(orbit.clusters))
            for ic in clus_range:
                clus = orbit.clusters[ic]
                ig = orbit.clusters_ig[ic]
                allclus.extend([cx, iorb, ig, ic] for cx in self.translate_cluster_to_supercell(sc, clus))
        return allclus

    @staticmethod
    def translate_cluster_to_supercell(sc, clus):
        """

        :param sc: SupercellStructure
        :return: all clusters in the supercell
        """
        if clus.order == 0:
            return [[]]*sc.n_cell
        coords= clus.frac_coords
        use_compiled_code = True
        if use_compiled_code:
            from f_util import f_util
            allclus= f_util.tr_cl_sc(sc._ijkl.T, sc.sc_mat.T, sc.inv_sc_mat.T, sc.sc_ref.T, sc.prim.frac_coords.T, clus._ijkls_np.T).T-1
        else:
            allclus = [get_structure_ordering((coords+ijk[None,:]).dot(sc.inv_sc_mat), sc.frac_coords,0) for ijk in sc.ijk_ref]
        return allclus


    def identify_cluster(self, clus_in):
        """
        Find if and which the input cluster matches the model clusters
        :param clus_try:
        :return: [matched cluster, iorb, ig, clus_sorted]
        """
        for iorb, orb in enumerate(self.orbits):
            # Quick return if not match obviously *)
            if orb.cluster._must_differ(clus_in):
                continue
            # logger.debug("  proceed with %d"%(iorb))
            for ic, clus_match in enumerate(orb.clusters):
                [flag, foundOrb] = clus_match.equivalent_by_lattranslation(clus_in.ijkls)
                if flag:
                    ig = orb.clusters_ig[ic]
                    # logger.debug("debug>  sum "+repr(clus_in)+" == orbit %d cluster %d sym %d matching="%(iorb, ic, ig) +str(foundOrb))
                    # logger.debug(" calc FC Trans mem=")#, MemoryInUse()/1E6)
                    return [True, iorb, ic, ig, relativePosition(clus_match.ijkls, foundOrb, dim=len(clus_in.orders))]
        return [False, -1, -1, None, None]


    def load_solution(self, sol_f):
        return np.loadtxt(sol_f, ndmin=2)


    def process_name_ord(self, name_ord, all_ord):
        if len(name_ord) == 0:
            name_ord = [['All', all_ord]]
        else:
            for i in name_ord:
                if i[0][-6:].lower()=='except':
                    i[1] = [j for j in all_ord if j not in i[1]]
                else:
                    i[1] = [j for j in i[1] if j in all_ord]
        return name_ord


    def get_submodels(self, name_ord, knownsol=None, **kwargs):
        """
        :param name_ord: list of [name, fct order], e.g. pair 0,2
        :return: list of [name, matrix] defining the different fittings
        """
        all_ord=list(set([o.cluster.order for o in self.orbits]))
        name_ord = self.process_name_ord(name_ord, all_ord)
        print("No. of parameters", {o: self.ord_idx[o+1]-self.ord_idx[o] for o in range(self.maxorder+1)})
        sol0 = np.zeros(self.ncorr)
        if knownsol:
             print("  Reading previous solution from %s"%(knownsol))
             input_sol = self.load_solution(knownsol).reshape(-1)
             sol0[:min(sol0.size, input_sol.size)] = input_sol[:min(sol0.size, input_sol.size)]
        return [[nm,
                scipy.sparse.identity(self.ncorr).tocsr()[:,
                  sum([list(range(self.ord_idx[i], self.ord_idx[i+1])) for i in o if i<=self.maxorder], [])], sol0]
                for (nm, o) in name_ord]


