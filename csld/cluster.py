#!/usr/bin/env python3

"""
This module provides classes used to define a 
  * site in a periodic structure for use in a cluster
  * cluster, proper or improper
  * orbit of cluster
"""


import itertools
import numpy as np
from collections import Counter

from .sites import Site, PeriodicSite
from .structure import Structure
from .lattice import Lattice
from .analyzer import SpacegroupAnalyzer
from .util.io_utils import read_into_lines, zopen, load_scmatrix
from .util.mathtool import multi_index_fac, remove_duplicates, MixedIndexForImproper, relativePosition
from .coord_utils import pbc_shortest_vectors
from .util.tool import list2index_unique


class ClusterSite(PeriodicSite):
    """
    lattice points, a periodic site plus
        ijkl ==[ijk, l]
          ijk: 3 integers for lattice points and
          l: index in primitive cell
          example [[0 0 0] 1]

        prim: symmetrized primitive cell

    """

    def __init__(self, ijkl, f_crd, prim):
        """
        set ijkl or f_crd to none to auto-generate the other
        """
        if not isinstance(prim, Structure):
            raise TypeError("set_prim need Structure, not", prim.__class__)
        self.prim = prim

        if ijkl is None and f_crd is None:
            raise TypeError("Supply at least one of ijkl or f_crd")
        if ijkl is None:
            self._ijkl = prim.frac2ijkl(f_crd)
        else:
            self._ijkl = np.array(ijkl, dtype=np.int)
            f_crd = prim.ijkl2frac(self._ijkl)
        self.ijkl = self._ijkl.tolist()
        self._orb_id = prim.orbit_of_l[self.ijkl[3]]
        PeriodicSite.__init__(self, '', f_crd, prim.lattice)


    def __hash__(self):
        """
        Possibly dangerous:
        i,j,k assumed to be between -500, 499, so each taking 3 digits (i+500)
        l assumed to be between 0 and 9999, so 4 digits
        """
        code=0
        n= self.ijkl[3]
        if n < 0 or n > 9999:
            raise ValueError("l cannot be hashed", self.l)
        code = code*10000 + n
        for x in self.ijkl[:3]:
            n = 500 + x
            if n < 0 or n > 999:
                raise ValueError("ijk cannot be hashed", self.ijk)
            code = code*1000 + n
        return code


    def __repr__(self):
        return str(self.ijkl)

    def __str__(self):
        return self.__repr__()

    def copy(self):
        """
        A copy of the cartesian coordinates of the site as a numpy array.
        """
        return ClusterSite(self._ijkl[:], self.f_crd[:], self.prim)


    def __lt__(self, other):
        """
        Sets a default sort order.
        1. by l, the index im primitive cell
        2. then by ijk
        """
        if self.ijkl[3] != other.ijkl[3]:
            return self.ijkl[3] < other.ijkl[3]
        return self.ijkl[:3] < other.ijkl[:3]

    def __eq__(self, other):
        """
        compare ijk and l
        """
        return self.ijkl == other.ijkl



class Cluster():
    """
    a collection of cluster sites. Can be either proper or improper
    """

    def __init__(self, vertices, prim):
        """
        Create a cluster.

        Args:
            vertices: cluster points
            prim: pointer to symmetrized primitive cell
        """
        self.vertices = vertices
        self.prim = prim
        self.order = len(vertices)
        self.orders = [self.order]
        self._diameter = None
        self._l_list = [v.ijkl[3] for v in vertices]
        self.update_books()
        self.order_uniq = len(self.uniq)
        self.factorial = multi_index_fac(self.multi_index)
        self._orb_id = sorted([v._orb_id for v in vertices])

    def update_books(self):
        counts = dict(Counter(self.vertices))
        self.uniq = list(counts.keys())
        self.multi_index = list(counts.values())

    def copy(self):
        return self.from_ijkl(self.ijkls, self.prim)

    @property
    def diameter(self):
        if self._diameter is None:
            diameter = 0
            for i in range(self.order_uniq):
                for j in range(i+1, self.order_uniq):
                    if not isinstance(self.uniq[i], int):
                        rij = np.linalg.norm(self.uniq[i].coords - self.uniq[j].coords)
                    else:
                        rij = np.linalg.norm(self.prim.coords[self.uniq[i]] - self.prim.coords[self.uniq[j]])
                    if diameter < rij:
                        diameter = rij
            self._diameter = diameter
        return self._diameter

    @property
    def center(self):
        # return self.prim.lattice.get_cartesian_coords()
        return np.average(self.frac_coords, axis=0)

# useful for filtering clusters
    def bond_counts(self, cutoff):
        """
        :param cutoff: max bond length
        :return: number of bonds between distinct atoms
        """
        nbond = 0
        for i in range(self.order_uniq):
            for j in range(i+1, self.order_uniq):
                rij = Site.dist(self.uniq[i], self.uniq[j])
                if rij <= cutoff:
                    nbond+=1
        return nbond


# is this the smallest cluster within ONE OF the supercells?
    def is_small_in(self, scmat,maxord=2):
        """
        Is the pair within the Wigner-Seitz cell of the supercell?
        :scmat1: 3x3 integer matrix or text file containing the matrix
        """
        sclist= [scmat] if isinstance(scmat, str) else scmat
        return any([self.is_small_in_1cell(sc, maxord) for sc in sclist])


# for each supercell
    def is_small_in_1cell(self, scmat1,maxord=2):
        """
        Is the pair within the Wigner-Seitz cell of the supercell?
        :scmat1: 3x3 integer matrix or text file containing the matrix
        """
# by default filter pairs only
        if self.order>maxord or self.order_uniq <2:
            return True
        #from .coord_utils import pbc_images
        scmat= load_scmatrix(scmat1, self.prim)
        latt = Lattice(np.dot(scmat, self.prim.lattice.matrix))
        #images = pbc_images(nogamma=True)
        #if self.uniq[0].ijkl==[0,0,0,4] and self.uniq[1].ijkl[-1]==11 and self.uniq[1].ijkl[0]==0 and self.uniq[1].ijkl[2]==0 and (self.diameter < 1.6 or self.diameter > 8):
         #   print('debug', self.ijkls)
        #elif self.order_uniq>=2:
        #    return False
        for i in range(self.order_uniq-1):
            for j in range(i+1, self.order_uniq):
               if not latt.vector_in_wigner_seitz(self.uniq[i].coords - self.uniq[j].coords, -1e-7):
                   return False
        return True


    def __repr__(self):
        outs = []
        for v, n in zip(self.uniq, self.multi_index):
            outs.append(v.__repr__() +  "*"+str(n) if n>1 else v.__repr__())
        return '['+ " ".join(outs) + ']'


    @staticmethod
    def cluster_name(n):
        """
        Translate to English names
        :param n:
        :return: name
        """
        names={0:"empty", 1:"point", 2:"pair", 3:"triplet", 4:"quadruplet", 5:"quintuplet",
               6:"sextuplet", 7:"septuplet", 8:"octuplet"}
        if n>=0 and n<=8:
            return names[n]
        else:
            return str(n)+"-body"

    @property
    def coords(self):
        return np.array([v.coords for v in self.vertices])

    @property
    def frac_coords(self):
        return np.array([v.frac_coords for v in self.vertices])

    def __str__(self):
        outs = ["cluster: order "+str(self.order), Cluster.cluster_name(self.order_uniq)]
        if self.order>1:
            outs.append("diameter %.4f" %(self.diameter))
        return " ".join(outs)

    @staticmethod
    def read_cluster(prim, clusf):
        """
        Format of cluster file from corrdump routine of Axel's ATAT
        Example
        ***************************
        1
        0.00000
        0

        1
        0.00000
        1
        1.00000 1.00000 1.00000 0 0

        1
        0.00000
        1
        0.50000 0.50000 0.50000 0 0

        6
        0.50000
        2
        0.50000 0.50000 0.50000 0 0
        0.00000 1.00000 0.00000 0 0

        ***************************
        """
        lines = read_into_lines(zopen(clusf, "r").read())
        lines = [x for x in lines if x != '']
        line = 0
        raw = []
        while line<len(lines):
            ord= int(lines[line+2])
            xyz = [list(map(float, aline.split()[:3])) for aline in lines[line+3:line+3+ord]]
            raw.append(Cluster.from_coords(xyz, prim))
            line+= 3 + ord
        return raw


    @staticmethod
    def remove_equivalent_clusters(clusters):
        return remove_duplicates(clusters, Cluster.compare_clusters_by_symmetry)

    @staticmethod
    def improper_from_proper(clus, n1, prim):
        """
        generate improper clusters by duplicating points in an proper one
        :param clus: proper
        :param ord:
        :return: list of generated improper ones
        """
        n0 = clus.order
        improp = []
        if n1<=n0 or n0<=0:
            raise ValueError("cannot extend order-",n0," cluster to order-", n1)

        # distribute n1 slots in the improper cluster to n0 sites in input proper cluster
        # or insert n0-1 bars among n1-n0 stars
        for bars in itertools.combinations(range(1, n1), n0 - 1):
            pos= [0] + list(bars) + [n1]
            cimp= []
            for i in range(n0):
#                print(i+1, "-th point expanded to", pos[i + 1]- pos[i])
                for _ in range(pos[i + 1]- pos[i]):
                    cimp.append(ClusterSite(clus.vertices[i]._ijkl, None, prim))
            improp.append( Cluster(cimp, prim))
#        print("improper generated=", improp)
        improp = Cluster.remove_equivalent_clusters(improp)
        return improp




    def _must_differ(self, other):
        """
        Compare if clusters are unambiguously different
        :param other:
        :return:
        """
        if self.order != other.order:
            return True
        elif self.order_uniq != other.order_uniq:
            return True
        elif abs(self.diameter - other.diameter) > 1E-6:
            return True
        elif self.factorial != other.factorial:
            return True
        elif self._orb_id != other._orb_id:
            return True
        else:
            return False

    def __lt__(self, other):
        """
        Sets a default sort order.
        1. by order
        2. then by order_uniq
        3. then by diameter
        4. then by -factorial, i.e. (aaab) before (aabb)
        """
        if self.order != other.order:
            return self.order < other.order
        if self.order_uniq != other.order_uniq:
            return self.order_uniq < other.order_uniq
        if abs(self.diameter - other.diameter) > 1E-6:
            return self.diameter < other.diameter
        return self.factorial > other.factorial

    @staticmethod
    def _sorted(v):
        return sorted(v)

    def equivalent_by_lattranslation(self, c2, tr_vec=False):
        """
        :param c2 is not an instantiated cluster, but a list of ijkl:
        :return: [bool, translated c2 or permuted c1]
        """
        if self.order != len(c2):
            return [False, []]
        trc2 = [v if isinstance(v, list) else v.tolist() for v in c2]

        sortc1 = self._sorted(self.ijkls)
        if sortc1 == self._sorted(trc2):
            return [True, trc2] +([[0,0,0]] if tr_vec else [])

        l1 = self.vertices[0].ijkl[3]
        ijk1 = self.vertices[0]._ijkl[:3]
        trijk = np.zeros(4, dtype=np.int)
        for i in range(self.order):
            if l1 != c2[i][3]:
                continue
            #  only the same INDEX of atoms might be related
            # print("self", self.__class__, "c2", c2)
            trijk[:3] = ijk1 - c2[i][:3]
            trc2 = [(trijk + v).tolist() for v in c2]
            # print("trc2", trc2)
            # print(sortc1 , sorted([[list(v[0]), v[1]] for v in trc2]))
            if sortc1 == self._sorted(trc2):
                return [True, trc2] +([trijk[:3]] if tr_vec else [])

        return [False, []]


    def operated_by(self, ig):
        """

        :param ig: index of symmetry operator
        :return:
        """
        return [self.prim.operate_ijkl(ig, v._ijkl) for v in self.vertices]



    def equivalent_by_symmetry(self, other):
        if self._must_differ(other):
            return False
        for ig in range(self.prim.nsymm):
            gclus = other.operated_by(ig)
            if self.equivalent_by_lattranslation(gclus)[0]:
                return True
        return False

    @staticmethod
    def compare_clusters_by_symmetry(c1, c2):
        return c1.equivalent_by_symmetry(c2)

    @staticmethod
    def _compare_clusters_ig_by_lattrans(c1, c2):
        return c1[1].equivalent_by_lattranslation(c2[1].ijkls)[0]

    @property
    def ijkls(self):
        return [v.ijkl for v in self.vertices]

    @property
    def l_list(self):
        return [i[-1] for i in self.ijkls]

    @property
    def _ijkls(self):
        return [v._ijkl for v in self.vertices]

    @property
    def _ijkls_np(self):
        return np.array([v._ijkl for v in self.vertices], dtype=int)

    @property
    def _ls(self):
        return [v.ijkl[-1] for v in self.vertices]

    @property
    def atomic_numbers_uniq(self):
        return list(set(np.array(self.prim.atomic_numbers)[self._ls].tolist()))


    @classmethod
    def from_ijkl(cls, ijkls, prim, *type):
        """
        Alternative Cluster constructor
        :param ijkls: ijkl of each vertex
        :return: Cluster
        """
        return cls([ClusterSite(ijkl, None, prim) for ijkl in ijkls], prim)

    def as_index(self):
        return list2index_unique(self.vertices, self.uniq)

    def append_site(self, ijkl):
        return Cluster(self.vertices + [ClusterSite(ijkl, None, self.prim)], self.prim)

    @classmethod
    def from_coords(cls, xyz, prim, frac_coords=True):
        """
        Alternative Cluster constructor
        :param xyz: coordinates
        :return: Cluster
        """
        coor = xyz if frac_coords else prim.lattice.get_fractional_coords(xyz)
        return cls([ClusterSite(None, np.array(c_f), prim) for c_f in coor], prim)

    def move_to_prim(self, dest=[0,0,0], dijk=None):
        """
        move first atom to destination cell
        """
        dijkl = np.zeros(4, dtype=np.int)
        if self.order > 0:
            if dijk is None:
                dijkl[:3] = dest - self.vertices[0]._ijkl[:3]
            else:
                dijkl[:3] = np.array(dijk)
        for i in range(self.order):
            self.vertices[i] = ClusterSite(self.vertices[i]._ijkl + dijkl, None, self.prim)
        self.update_books()

    def move_closest_to(self, pt):
        """

        :param pt: fractional coordinate
        :return:
        """
        lat = self.prim.lattice
        center = self.center
        min_vec = lat.get_fractional_coords(pbc_shortest_vectors(lat, center, pt))
        self.move_to_prim(dijk=((pt-center)-min_vec).astype(int))


    def index_cluster_in_supercell(self, sc):
        """
        indexing a single cluster
        :param sc: supercell
        :return:
        """
        SCref = sc.sc_ref.tolist()
        invSCmat = sc.inv_sc_mat
        n_cell = sc.n_cell
        scmat = sc.sc_mat
        return [SCref.index(np.dot(np.mod(np.dot(pt[:3], invSCmat), 1), scmat).astype(np.int).tolist()) +
                n_cell*pt[3] for pt in self.ijkls]


    def reduce_improper_fct_output(self, fct):
        """
        set to zero the rest of the FCT elements equivalent by symmetric index of IMPROPER clusters
        example: cluster {a,a}, the FCT will be redefined as
        fxx 2fxy 2fxz
        0    fyy 2fyz
        0     0   fzz
        fct: the FCT
        """
        if self.factorial <= 1:
            return fct
        # this subroutine is bugged. For safety do nothing
        return fct
        idxlist = MixedIndexForImproper(self.vertices, 3)[1]
        fctnew = fct
        for i in idxlist:
            if len(i) > 1:
                fctnew[i[0]] *= len(i)
                for j in i[1:]:
                     fctnew[j] = 0
   # (*Print["in, out=",Transpose[{fct, fctnew}]];*)
        return fctnew

    @staticmethod
    def clusters_to_orbit_isotropy(prim, clusters):
        orbits = []
        for clus in clusters:
            iso = []
            orbit = [[0,clus]]
####            # assume identity operator is the very first!
            for ig in range(prim.nsymm):
                gclus = clus.operated_by(ig)
                if clus.order>0 and len(clus.orders)<=1 and hasattr(prim, 'spacegroup'):
                    # move to center of mass
                    ijk0 = np.append(np.mean(np.array(gclus)[:,:3], axis=0).astype(np.int), 0)
                    # print(ijk0, clus, gclus)
                    gclus = [v - ijk0 for v in gclus]
                [flag, gp] = clus.equivalent_by_lattranslation(gclus)
                if flag:
                    iso.append([ig, relativePosition(clus.ijkls, gp, dim=len(clus.orders))])
                else:
                    orbit.append([ig, clus.from_ijkl(gclus, prim, clus)])

            orbit = remove_duplicates(orbit, Cluster._compare_clusters_ig_by_lattrans)
            # check for Lagrange's theorem
            assert len(iso)*len(orbit) == prim.nsymm, TypeError('ERROR is, orbit nsymm=%d %d %d'%(len(iso), len(orbit), prim.nsymm))
            if hasattr(prim, 'spacegroup'):
                pointgroup_symbol = SpacegroupAnalyzer.get_point_group_symbol([prim.spacegroup[ig].rot_f for ig, pi in iso])
            else:
                pointgroup_symbol = '1  (TBD)'
            orbits.append(ClusterOrbit(orbit, iso, pointgroup_symbol))
            # print("orbit", orbit)
            # print("iso", iso)
        # print("  multiplicity", [len(orbiso[0]) for orbiso in orbits])
        return orbits


    def index_of_uniq(self):
        return [self.uniq.index(v) for v in self.vertices]

    def permutations(self):
        from itertools import permutations
        uniq_ids = self.index_of_uniq()
        permutated_ids = []
        perms=[]
        for iper in permutations(range(self.order)):
            id_permutated = [uniq_ids[i] for i in iper]
            if id_permutated not in permutated_ids:
                perms.append(iper)
                permutated_ids.append(id_permutated)
        return perms


class ClusterOrbit():
    """
    orbit of an cluster, consisting of representative cluster and the set of
    clustered the spacegroup symmops map the cluster into
    """
    def __init__(self, orbit, iso, pointgroup_symbol='TBD'):
        self.clusters_ig, self.clusters = list(zip(*orbit))
        self.cluster = self.clusters[0]
        self.isotropy = iso
        self.multiplicity = len(self.clusters)
        self.pointgroup_symbol = pointgroup_symbol
        self.ncorr_full = -1
        self.order = self.cluster.order
