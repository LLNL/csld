#!/usr/bin/env python3

"""
This module provides classes used to define a periodic structure.
"""
# adapted from original version in pymatgen version from pymatgen


import math
import collections
import itertools
from abc import ABCMeta, abstractmethod, abstractproperty
import warnings

import numpy as np

from fractions import gcd
from .lattice import Lattice
from .util.periodic_table import Element, Specie, get_el_sp
from .sites import PeriodicSite
from .composition import Composition
from .coord_utils import get_angle, all_distances
from .util.units import Mass, Length
from .coord_utils import supercell_latticepoints
from .util.tool import non_1to1


class SiteCollection(collections.Sequence):
    """
    Basic SiteCollection. Essentially a sequence of Sites or PeriodicSites.
    This serves as a base class for Molecule (a collection of Site, i.e., no
    periodicity) and Structure (a collection of PeriodicSites, i.e.,
    periodicity). Not meant to be instantiated directly.
    """
    __metaclass__ = ABCMeta

    #Tolerance in Angstrom for determining if sites are too close.
    DISTANCE_TOLERANCE = 0.01

    @abstractproperty
    def sites(self):
        """
        Returns a tuple of sites.
        """
        return

    @abstractmethod
    def get_distance(self, i, j):
        """
        Returns distance between sites at index i and j.

        Args:
            i (int): Index of first site
            j (int): Index of second site

        Returns:
            (float) Distance between sites at index i and index j.
        """
        return

    @property
    def distance_matrix(self):
        """
        Returns the distance matrix between all sites in the structure. For
        periodic structures, this is overwritten to return the nearest image
        distance.
        """
        return all_distances(self.cart_coords, self.cart_coords)

    @property
    def species(self):
        """
        Only works for ordered structures.
        Disordered structures will raise an AttributeError.

        Returns:
            ([Specie]) List of species at each site of the structure.
        """
        return [site.specie for site in self]

    @property
    def elements(self):
        return [site.specie.__str__() for site in self]

    @property
    def species_and_occu(self):
        """
        List of species and occupancies at each site of the structure.
        """
        return [site.species_and_occu for site in self]

    @property
    def ntypesp(self):
        """Number of types of atoms."""
        return len(self.types_of_specie)

    @property
    def types_of_specie(self):
        """
        List of types of specie. Only works for ordered structures.
        Disordered structures will raise an AttributeError.
        """
        # Cannot use set since we want a deterministic algorithm.
        types = []
        for site in self:
            if site.specie not in types:
                types.append(site.specie)
        return types

    @property
    def types_of_elements(self):
        return [i.symbol for i in self.types_of_specie]


    def group_by_types(self):
        """Iterate over species grouped by type"""
        for t in self.types_of_specie:
            for site in self:
                if site.specie == t:
                    yield site

    def indices_from_symbol(self, symbol):
        """
        Returns a tuple with the sequential indices of the sites
        that contain an element with the given chemical symbol.
        """
        indices = []
        for i, specie in enumerate(self.species):
            if specie.symbol == symbol:
                indices.append(i)
        return tuple(indices)

    @property
    def symbol_set(self):
        """
        Tuple with the set of chemical symbols.
        Note that len(symbol_set) == len(types_of_specie)
        """
        return tuple([specie.symbol for specie in self.types_of_specie])

    @property
    def atomic_masses(self):
        """List of atomic masses."""
        return [site.specie.atomic_mass for site in self]

    @property
    def atomic_numbers(self):
        """List of atomic numbers."""
        return [site.specie.number for site in self]

    @property
    def site_properties(self):
        """
        Returns the site properties as a dict of sequences. E.g.,
        {"magmom": (5,-5), "charge": (-4,4)}.
        """
        props = collections.defaultdict(list)
        for site in self:
            for k, v in site.properties.items():
                props[k].append(v)
        return props

    def __contains__(self, site):
        return site in self.sites

    def __iter__(self):
        return self.sites.__iter__()

    def __getitem__(self, ind):
        return self.sites[ind]

    def __len__(self):
        return len(self.sites)

    def __hash__(self):
        #for now, just use the composition hash code.
        return self.composition.__hash__()

    @property
    def num_sites(self):
        """
        Number of sites.
        """
        return len(self)

    @property
    def cart_coords(self):
        """
        Returns a list of the cartesian coordinates of sites in the structure.
        """
        return np.array([site.coords for site in self])

    @property
    def formula(self):
        """
        (str) Returns the formula.
        """
        return self.composition.formula

    @property
    def composition(self):
        """
        (Composition) Returns the composition
        """
        elmap = collections.defaultdict(float)
        for site in self:
            for species, occu in site.species_and_occu.items():
                elmap[species] += occu
        return Composition(elmap)

    @property
    def charge(self):
        """
        Returns the net charge of the structure based on oxidation states. If
        Elements are found, a charge of 0 is assumed.
        """
        charge = 0
        for site in self:
            for specie, amt in site.species_and_occu.items():
                charge += getattr(specie, "oxi_state", 0) * amt
        return charge

    @property
    def is_ordered(self):
        """
        Checks if structure is ordered, meaning no partial occupancies in any
        of the sites.
        """
        return all((site.is_ordered for site in self))

    def get_angle(self, i, j, k):
        """
        Returns angle specified by three sites.

        Args:
            i (int): Index of first site.
            j (int): Index of second site.
            k (int): Index of third site.

        Returns:
            (float) Angle in degrees.
        """
        v1 = self[i].coords - self[j].coords
        v2 = self[k].coords - self[j].coords
        return get_angle(v1, v2, units="degrees")

    def get_dihedral(self, i, j, k, l):
        """
        Returns dihedral angle specified by four sites.

        Args:
            i (int): Index of first site
            j (int): Index of second site
            k (int): Index of third site
            l (int): Index of fourth site

        Returns:
            (float) Dihedral angle in degrees.
        """
        v1 = self[k].coords - self[l].coords
        v2 = self[j].coords - self[k].coords
        v3 = self[i].coords - self[j].coords
        v23 = np.cross(v2, v3)
        v12 = np.cross(v1, v2)
        return math.degrees(math.atan2(np.linalg.norm(v2) * np.dot(v1, v23),
                            np.dot(v12, v23)))

    def is_valid(self, tol=DISTANCE_TOLERANCE):
        """
        True if SiteCollection does not contain atoms that are too close
        together. Note that the distance definition is based on type of
        SiteCollection. Cartesian distances are used for non-periodic
        Molecules, while PBC is taken into account for periodic structures.

        Args:
            tol (float): Distance tolerance. Default is 0.01A.

        Returns:
            (bool) True if SiteCollection does not contain atoms that are too
            close together.
        """
        if len(self.sites) == 1:
            return True
        all_dists = self.distance_matrix[np.triu_indices(len(self), 1)]
        return bool(np.min(all_dists) > tol)


class IStructure(SiteCollection):
    """
    Basic immutable Structure object with periodicity. Essentially a sequence
    of PeriodicSites having a common lattice. IStructure is made to be
    (somewhat) immutable so that they can function as keys in a dict. To make
    modifications, use the standard Structure object instead. Structure
    extends Sequence and Hashable, which means that in many cases,
    it can be used like any Python sequence. Iterating through a
    structure is equivalent to going through the sites in sequence.
    """

    def __init__(self, lattice, species, coords, validate_proximity=False,
                 to_unit_cell=False, coords_are_cartesian=False,
                 site_properties=None,site_properties_T=None,
                 intensive_properties={}, extensive_properties={}):
        """
        Create a periodic structure.

        Args:
            lattice (Lattice/3x3 array): The lattice, either as a
                :class:`pymatgen.core.lattice.Lattice` or
                simply as any 2D array. Each row should correspond to a lattice
                vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
                lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].
            species ([Specie]): Sequence of species on each site. Can take in
                flexible input, including:

                i.  A sequence of element / specie specified either as string
                    symbols, e.g. ["Li", "Fe2+", "P", ...] or atomic numbers,
                    e.g., (3, 56, ...) or actual Element or Specie objects.

                ii. List of dict of elements/species and occupancies, e.g.,
                    [{"Fe" : 0.5, "Mn":0.5}, ...]. This allows the setup of
                    disordered structures.
            fractional_coords (Nx3 array): list of fractional coordinates of
                each species.
            validate_proximity (bool): Whether to check if there are sites
                that are less than 0.01 Ang apart. Defaults to False.
            coords_are_cartesian (bool): Set to True if you are providing
                coordinates in cartesian coordinates. Defaults to False.
            site_properties (dict): Properties associated with the sites as a
                dict of sequences, e.g., {"magmom":[5,5,5,5]}. The sequences
                have to be the same length as the atomic species and
                fractional_coords. Defaults to None for no properties.

            site_properties_T (list): alternative way to specify site_properties
                essentially site_properties transposed
        """
        if len(species) != len(coords):
            raise StructureError("The list of atomic species must be of the"
                                 "same length as the list of fractional"
                                 " coordinates.")

        if isinstance(lattice, Lattice):
            self._lattice = lattice
        else:
            self._lattice = Lattice(lattice)

        sites = []
        for i in range(len(species)):
            prop = None
            if site_properties:
                prop = {k: v[i] for k, v in site_properties.items()}
            elif site_properties_T:
                prop = site_properties_T[i]
            sites.append(
                PeriodicSite(species[i], coords[i], self._lattice,
                             to_unit_cell,
                             coords_are_cartesian=coords_are_cartesian,
                             properties=prop))
        self._sites = tuple(sites)
        if validate_proximity and not self.is_valid():
            raise StructureError(("Structure contains sites that are ",
                                  "less than 0.01 Angstrom apart!"))
        self.intensive_properties=intensive_properties
        self.extensive_properties=extensive_properties


    @classmethod
    def from_sites(cls, sites, validate_proximity=False,
                   to_unit_cell=False):
        """
        Convenience constructor to make a Structure from a list of sites.

        Args:
            sites: Sequence of PeriodicSites. Sites must have the same
                lattice.
            validate_proximity (bool): Whether to check if there are sites
                that are less than 0.01 Ang apart. Defaults to False.
            to_unit_cell (bool): Whether to translate sites into the unit
                cell.

        Returns:
            (Structure) Note that missing properties are set as None.
        """
        prop_keys = []
        props = {}
        lattice = None
        for i, site in enumerate(sites):
            if not lattice:
                lattice = site.lattice
            elif site.lattice != lattice:
                raise ValueError("Sites must belong to the same lattice")
            for k, v in site.properties.items():
                if k not in prop_keys:
                    prop_keys.append(k)
                    props[k] = [None] * len(sites)
                props[k][i] = v
        for k, v in props.items():
            if any((vv is None for vv in v)):
                warnings.warn("Not all sites have property %s. Missing values "
                              "are set to None." % k)
        return cls(lattice, [site.species_and_occu for site in sites],
                   [site.frac_coords for site in sites],
                   site_properties=props,
                   validate_proximity=validate_proximity,
                   to_unit_cell=to_unit_cell)

    @property
    def distance_matrix(self):
        """
        Returns the distance matrix between all sites in the structure. For
        periodic structures, this should return the nearest image distance.
        """
        return self.lattice.get_all_distances(self.frac_coords,
                                              self.frac_coords)

    @property
    def distance_matrix_noself(self):
        return self.lattice.get_all_distances(self.frac_coords,
                                              self.frac_coords, nogamma=True)

    @property
    def sites(self):
        """
        Returns an iterator for the sites in the Structure.
        """
        return self._sites

    @property
    def lattice(self):
        """
        Lattice of the structure.
        """
        return self._lattice

    @property
    def reciprocal_lattice(self):
        """
        Reciprocal lattice of the structure.
        """
        return self._lattice.reciprocal_lattice

    def lattice_vectors(self, space="r"):
        """
        Returns the vectors of the unit cell in Angstrom.

        Args:
            space: "r" for real space vectors, "g" for reciprocal space basis
                vectors.
        """
        if space.lower() == "r":
            return self.lattice.matrix
        if space.lower() == "g":
            return self.lattice.reciprocal_lattice.matrix
        raise ValueError("Wrong value for space: %s " % space)

    @property
    def density(self):
        """
        Returns the density in units of g/cc
        """
        m = Mass(self.composition.weight, "amu")
        return m.to("g") / (self.volume * Length(1, "ang").to("cm") ** 3)

    def __eq__(self, other):
        if other is None:
            return False
        if len(self) != len(other):
            return False
        if self._lattice != other._lattice:
            return False
        for site in self:
            if site not in other:
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        # For now, just use the composition hash code.
        return self.composition.__hash__()

    @property
    def frac_coords(self):
        """
        Fractional coordinates as a Nx3 numpy array.
        """
        return np.array([site.frac_coords for site in self._sites])

    @property
    def volume(self):
        """
        Returns the volume of the structure.
        """
        return self._lattice.volume

    def get_distance(self, i, j, jimage=None):
        """
        Get distance between site i and j assuming periodic boundary
        conditions. If the index jimage of two sites atom j is not specified it
        selects the jimage nearest to the i atom and returns the distance and
        jimage indices in terms of lattice vector translations if the index
        jimage of atom j is specified it returns the distance between the i
        atom and the specified jimage atom.

        Args:
            i (int): Index of first site
            j (int): Index of second site
            jimage: Number of lattice translations in each lattice direction.
                Default is None for nearest image.

        Returns:
            distance
        """
        return self[i].distance(self[j], jimage)

    def distance_and_image(self, i, j, jimage=None):
        return self[i].distance_and_image(self[j], jimage)

    def get_sites_in_sphere(self, pt, r, include_index=False):
        """
        Find all sites within a sphere from the point. This includes sites
        in other periodic images.

        Algorithm:

        1. place sphere of radius r in crystal and determine minimum supercell
           (parallelpiped) which would contain a sphere of radius r. for this
           we need the projection of a_1 on a unit vector perpendicular
           to a_2 & a_3 (i.e. the unit vector in the direction b_1) to
           determine how many a_1"s it will take to contain the sphere.

           Nxmax = r * length_of_b_1 / (2 Pi)

        2. keep points falling within r.

        Args:
            pt (3x1 array): cartesian coordinates of center of sphere.
            r (float): Radius of sphere.
            include_index (bool): Whether the non-supercell site index
                is included in the returned data

        Returns:
            [(site, dist) ...] since most of the time, subsequent processing
            requires the distance.
        """
        site_fcoords = np.mod(self.frac_coords, 1)
        neighbors = []
        for fcoord, dist, i in self._lattice.get_points_in_sphere(
                site_fcoords, pt, r):
            nnsite = PeriodicSite(self[i].species_and_occu,
                                  fcoord, self._lattice,
                                  properties=self[i].properties)
            neighbors.append((nnsite, dist) if not include_index
                             else (nnsite, dist, i))
        return neighbors

    def get_neighbors(self, site, r, include_index=False):
        """
        Get all neighbors to a site within a sphere of radius r.  Excludes the
        site itself.

        Args:
            site:
                site, which is the center of the sphere.
            r:
                radius of sphere.
            include_index:
                boolean that determines whether the non-supercell site index
                is included in the returned data

        Returns:
            [(site, dist) ...] since most of the time, subsequent processing
            requires the distance.
        """
        nn = self.get_sites_in_sphere(site.coords, r,
                                      include_index=include_index)
        return [d for d in nn if site != d[0]]

    def get_all_neighbors(self, r, include_index=False):
        """
        Get neighbors for each atom in the unit cell, out to a distance r
        Returns a list of list of neighbors for each site in structure.
        Use this method if you are planning on looping over all sites in the
        crystal. If you only want neighbors for a particular site, use the
        method get_neighbors as it may not have to build such a large supercell
        However if you are looping over all sites in the crystal, this method
        is more efficient since it only performs one pass over a large enough
        supercell to contain all possible atoms out to a distance r.
        The return type is a [(site, dist) ...] since most of the time,
        subsequent processing requires the distance.

        Args:
            r (float): Radius of sphere.
            include_index (bool): Whether to include the non-supercell site
                in the returned data

        Returns:
            A list of a list of nearest neighbors for each site, i.e.,
            [[(site, dist, index) ...], ..]
            Index only supplied if include_index = True.
            The index is the index of the site in the original (non-supercell)
            structure. This is needed for ewaldmatrix by keeping track of which
            sites contribute to the ewald sum.
        """

        # Use same algorithm as get_sites_in_sphere to determine supercell but
        # loop over all atoms in crystal
        recp_len = self.lattice.reciprocal_lattice.abc
        sr = r + 0.15
        nmax = [sr * l / (2 * math.pi) for l in recp_len]
        site_nminmax = []
        floor = math.floor
        inds = (0, 1, 2)
        for site in self:
            pcoords = site.frac_coords
            inmax = [int(floor(pcoords[i] + nmax[i])) for i in inds]
            inmin = [int(floor(pcoords[i] - nmax[i])) for i in inds]
            site_nminmax.append(zip(inmin, inmax))

        nmin = [min([i[j][0] for i in site_nminmax]) for j in inds]
        nmax = [max([i[j][1] for i in site_nminmax]) for j in inds]

        all_ranges = [range(nmin[i], nmax[i] + 1) for i in inds]

        neighbors = [list() for i in range(len(self._sites))]
        all_fcoords = np.mod(self.frac_coords, 1)

        site_coords = np.array(self.cart_coords)
        latt = self._lattice
        frac_2_cart = latt.get_cartesian_coords
        n = len(self)
        indices = np.array(range(n))
        for image in itertools.product(*all_ranges):
            for (j, fcoord) in enumerate(all_fcoords):
                fcoords = fcoord + image
                coords = frac_2_cart(fcoords)
                submat = np.tile(coords, (n, 1))
                dists = np.power(site_coords - submat, 2)
                dists = np.sqrt(dists.sum(axis=1))
                withindists = (dists <= r) * (dists > 1e-8)
                sp = self[j].species_and_occu
                props = self[j].properties
                for i in indices[withindists]:
                    nnsite = PeriodicSite(sp, fcoords, latt,
                                          properties=props)
                    item = (nnsite, dists[i], j) if include_index else (
                        nnsite, dists[i])
                    neighbors[i].append(item)
        return neighbors

    def get_neighbors_in_shell(self, origin, r, dr):
        """
        Returns all sites in a shell centered on origin (coords) between radii
        r-dr and r+dr.

        Args:
            origin (3x1 array): Cartesian coordinates of center of sphere.
            r (float): Inner radius of shell.
            dr (float): Width of shell.

        Returns:
            [(site, dist) ...] since most of the time, subsequent processing
            requires the distance.
        """
        outer = self.get_sites_in_sphere(origin, r + dr)
        inner = r - dr
        return [(site, dist) for (site, dist) in outer if dist > inner]

    def get_sorted_structure(self):
        """
        Get a sorted copy of the structure.
        Sites are sorted by the electronegativity of the species.
        """
#        sites = sorted(self)
# WARNING: Sorting by electronegativity was NOT implemented?????
        sites= self.sites
        return self.__class__.from_sites(sites)

    def get_reduced_structure(self, reduction_algo="niggli"):
        """
        Get a reduced structure.

        Args:
            reduction_algo (str): The lattice reduction algorithm to use.
                Currently supported options are "niggli" or "LLL".
        """
        if reduction_algo == "niggli":
            reduced_latt = self._lattice.get_niggli_reduced_lattice()
        elif reduction_algo == "LLL":
            reduced_latt = self._lattice.get_lll_reduced_lattice()
        else:
            raise ValueError("Invalid reduction algo : {}"
                             .format(reduction_algo))

        return self.__class__(reduced_latt, self.species_and_occu,
                              self.cart_coords,
                              coords_are_cartesian=True, to_unit_cell=True)

    def copy(self, site_properties=None, sanitize=False):
        """
        Convenience method to get a copy of the structure, with options to add
        site properties.

        Args:
            site_properties (dict): Properties to add or override. The
                properties are specified in the same way as the constructor,
                i.e., as a dict of the form {property: [values]}. The
                properties should be in the order of the *original* structure
                if you are performing sanitization.
            sanitize (bool): If True, this method will return a sanitized
                structure. Sanitization performs a few things: (i) The sites are
                sorted by electronegativity, (ii) a LLL lattice reduction is
                carried out to obtain a relatively orthogonalized cell,
                (iii) all fractional coords for sites are mapped into the
                unit cell.

        Returns:
            A copy of the Structure, with optionally new site_properties and
            optionally sanitized.
        """
        props = self.site_properties
        if site_properties:
            props.update(site_properties)
        if not sanitize:
            return self.__class__(self._lattice,
                                  self.species_and_occu,
                                  self.frac_coords,
                                  site_properties=props)
        else:
            reduced_latt = self._lattice.get_lll_reduced_lattice()
            new_sites = []
            for i, site in enumerate(self):
                frac_coords = reduced_latt.get_fractional_coords(site.coords)
                site_props = {}
                for p in props:
                    site_props[p] = props[p][i]
                new_sites.append(PeriodicSite(site.species_and_occu,
                                              frac_coords, reduced_latt,
                                              to_unit_cell=True,
                                              properties=site_props))
            new_sites = sorted(new_sites)
            return self.__class__.from_sites(new_sites)

    def interpolate(self, end_structure, nimages=10,
                    interpolate_lattices=False, pbc=True):
        """
        Interpolate between this structure and end_structure. Useful for
        construction of NEB inputs.

        Args:
            end_structure (Structure): structure to interpolate between this
                structure and end.
            nimages (int): No. of interpolation images. Defaults to 10 images.
            interpolate_lattices (bool): Whether to interpolate the lattices.
                Interpolates the lengths and angles (rather than the matrix)
                so orientation may be affected.
            pbc (bool): Whether to use periodic boundary conditions to find
                the shortest path between endpoints.

        Returns:
            List of interpolated structures. The starting and ending
            structures included as the first and last structures respectively.
            A total of (nimages + 1) structures are returned.
        """
        #Check length of structures
        if len(self) != len(end_structure):
            raise ValueError("Structures have different lengths!")

        if interpolate_lattices:
            #interpolate lattices
            lstart = np.array(self.lattice.lengths_and_angles)
            lend = np.array(end_structure.lattice.lengths_and_angles)
            lvec = lend - lstart

        #Check that both structures have the same lattice
        elif not self.lattice == end_structure.lattice:
            raise ValueError("Structures with different lattices!")

        #Check that both structures have the same species
        for i in range(0, len(self)):
            if self[i].species_and_occu != end_structure[i].species_and_occu:
                raise ValueError("Different species!\nStructure 1:\n" +
                                 str(self) + "\nStructure 2\n" +
                                 str(end_structure))

        start_coords = np.array(self.frac_coords)
        end_coords = np.array(end_structure.frac_coords)
        vec = end_coords - start_coords
        if pbc:
            vec -= np.round(vec)
        sp = self.species_and_occu
        structs = []
        for x in range(nimages+1):
            if interpolate_lattices:
                l_a = lstart + x / nimages * lvec
                l = Lattice.from_lengths_and_angles(*l_a)
            else:
                l = self.lattice
            fcoords = start_coords + x / nimages * vec
            structs.append(self.__class__(l, sp, fcoords,
                           site_properties=self.site_properties))
        return structs

    def get_primitive_structure(self, tolerance=0.25):
        """
        This finds a smaller unit cell than the input. Sometimes it doesn"t
        find the smallest possible one, so this method is recursively called
        until it is unable to find a smaller cell.

        The method works by finding possible smaller translations
        and then using that translational symmetry instead of one of the
        lattice basis vectors if more than one vector is found (usually the
        case for large cells) the one with the smallest norm is used.

        Things are done in fractional coordinates because its easier to
        translate back to the unit cell.

        NOTE: if the tolerance is greater than 1/2 the minimum inter-site
        distance, the algorithm may find 2 non-equivalent sites that are
        within tolerance of each other. The algorithm will reject this
        lattice.

        Args:
            tolerance (float): Tolerance for each coordinate of a particular
                site. For example, [0.5, 0, 0.5] in cartesian coordinates
                will be considered to be on the same coordinates as
                [0, 0, 0] for a tolerance of 0.5. Defaults to 0.5.

        Returns:
            The most primitive structure found. The returned structure is
            guaranteed to have len(new structure) <= len(structure).
        """
        original_volume = self.volume

        #get the possible symmetry vectors
        sites = sorted(self._sites, key=lambda site: site.species_string)
        grouped_sites = [list(a[1]) for a
                         in itertools.groupby(sites,
                                              key=lambda s: s.species_string)]

        num_fu = reduce(gcd, map(len, grouped_sites))
        min_vol = original_volume * 0.5 / num_fu

        min_site_list = min(grouped_sites, key=lambda group: len(group))

        min_site_list = [site.to_unit_cell for site in min_site_list]
        org = min_site_list[0].coords
        possible_vectors = [min_site_list[i].coords - org
                            for i in range(1, len(min_site_list))]

        #Let's try to use the shortest vector possible first. Allows for faster
        #convergence to primitive cell.
        possible_vectors = sorted(possible_vectors,
                                  key=lambda x: np.linalg.norm(x))

        # Pre-create a few varibles for faster lookup.
        all_coords = [site.coords for site in sites]
        all_sp = [site.species_and_occu for site in sites]
        new_structure = None

        #all lattice points need to be projected to 0 under new basis
        l_points = np.array([[0, 0, 1], [0, 1, 0], [0, 1, 1], [1, 0, 0],
                             [1, 0, 1], [1, 1, 0], [1, 1, 1]])
        l_points = self._lattice.get_cartesian_coords(l_points)

        for v, repl_pos in itertools.product(possible_vectors, range(3)):
            #Try combinations of new lattice vectors with existing lattice
            #vectors.
            latt = self._lattice.matrix
            latt[repl_pos] = v

            #Exclude coplanar lattices from consideration.
            if abs(np.dot(np.cross(latt[0], latt[1]), latt[2])) < min_vol:
                continue
            latt = Lattice(latt)

            #Convert to fractional tol
            tol = tolerance / np.array(latt.abc)

            #check validity of new basis
            new_l_points = latt.get_fractional_coords(l_points)
            f_l_dist = np.abs(new_l_points - np.round(new_l_points))
            if np.any(f_l_dist > tol[None, None, :]):
                continue

            all_frac = latt.get_fractional_coords(np.array(all_coords))

            #calculate grouping of equivalent sites, represented by
            #adjacency matrix
            fdist = all_frac[None, :, :] - all_frac[:, None, :]
            fdist = np.abs(fdist - np.round(fdist))
            groups = np.all(fdist < tol[None, None, :], axis=2)

            #check that all group sizes are the same
            sizes = np.unique(np.sum(groups, axis=0))
            if len(sizes) > 1:
                continue

            #check that reduction in number of sites was by the same
            #amount as the volume reduction
            if round(self._lattice.volume / latt.volume) != sizes[0]:
                continue

            new_sp = []
            new_frac = []
            #this flag is set to ensure that all sites in a group are
            #the same species, it is set to false if a group is found
            #where this is not the case
            correct = True

            added = np.zeros(len(groups), dtype='bool')
            for i, g in enumerate(groups):
                if added[i]:
                    continue
                indices = np.where(g)[0]
                i0 = indices[0]
                sp = all_sp[i0]
                added[indices] = 1
                if not all([all_sp[i] == sp for i in indices]):
                    correct = False
                    break
                new_sp.append(all_sp[i0])
                new_frac.append(all_frac[i0])

            if correct:
                new_structure = self.__class__(
                    latt, new_sp, new_frac, to_unit_cell=True)
                break

        if new_structure and len(new_structure) != len(self):
            # If a more primitive structure has been found, try to find an
            # even more primitive structure again.
            return new_structure.get_primitive_structure(tolerance=tolerance)
        else:
            return self

    def __repr__(self):
        outs = ["Structure Summary", repr(self.lattice)]
        for s in self:
            outs.append(repr(s))
        return "\n".join(outs)

    def __str__(self):
        outs = ["Structure Summary ({s})".format(s=str(self.composition)),
                "Reduced Formula: {}"
                .format(self.composition.reduced_formula)]
        to_s = lambda x: "%0.6f" % x
        outs.append("abc   : " + " ".join([to_s(i).rjust(10)
                                           for i in self.lattice.abc]))
        outs.append("angles: " + " ".join([to_s(i).rjust(10)
                                           for i in self.lattice.angles]))
        outs.append("Sites ({i})".format(i=len(self)))
        for i, site in enumerate(self):
            outs.append(" ".join([str(i + 1), site.species_string,
                                  " ".join([to_s(j).rjust(12)
                                            for j in site.frac_coords])]))
        return "\n".join(outs)

    @property
    def to_dict(self):
        """
        Json-serializable dict representation of Structure
        """
        latt_dict = self._lattice.to_dict
        del latt_dict["@module"]
        del latt_dict["@class"]

        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__,
             "lattice": latt_dict, "sites": []}
        for site in self:
            site_dict = site.to_dict
            del site_dict["lattice"]
            del site_dict["@module"]
            del site_dict["@class"]
            d["sites"].append(site_dict)
        return d

    @classmethod
    def from_dict(cls, d):
        """
        Reconstitute a Structure object from a dict representation of Structure
        created using to_dict.

        Args:
            d (dict): Dict representation of structure.

        Returns:
            Structure object
        """
        lattice = Lattice.from_dict(d["lattice"])
        sites = [PeriodicSite.from_dict(sd, lattice) for sd in d["sites"]]
        return cls.from_sites(sites)

    def dot(self, coords_a, coords_b, space="r", frac_coords=False):
        """
        Compute the scalar product of vector(s) either in real space or
        reciprocal space.

        Args:
            coords (3x1 array): Array-like object with the coordinates.
            space (str): "r" for real space, "g" for reciprocal space.
            frac_coords (bool): Whether the vector corresponds to fractional or
                cartesian coordinates.

        Returns:
            one-dimensional `numpy` array.
        """
        lattice = {"r": self.lattice,
                   "g": self.reciprocal_lattice}[space.lower()]
        return lattice.dot(coords_a, coords_b, frac_coords=frac_coords)

    def norm(self, coords, space="r", frac_coords=True):
        """
        Compute the norm of vector(s) either in real space or reciprocal space.

        Args:
            coords (3x1 array): Array-like object with the coordinates.
            space (str): "r" for real space, "g" for reciprocal space.
            frac_coords (bool): Whether the vector corresponds to fractional or
                cartesian coordinates.

        Returns:
            one-dimensional `numpy` array.
        """
        return np.sqrt(self.dot(coords, coords, space=space,
                                frac_coords=frac_coords))



class Structure(IStructure, collections.MutableSequence):
    """
    Mutable version of structure. Much easier to use for editing,
    but cannot be used as a key in a dict.
    """
    __hash__ = None

    def __init__(self, lattice, species, coords, validate_proximity=False,
                 to_unit_cell=False, coords_are_cartesian=False,
                 site_properties=None,site_properties_T=None,
                 intensive_properties={}, extensive_properties={}):
        """
        Create a periodic structure.

        Args:
            scale: scaling factor, real number
            R: lattice vectors in rows. Note R*scale == lattice!!!
            lattice: The lattice, either as a pymatgen.core.lattice.Lattice or
                simply as any 2D array. Each row should correspond to a lattice
                vector. E.g., [[10,0,0], [20,10,0], [0,0,30]] specifies a
                lattice with lattice vectors [10,0,0], [20,10,0] and [0,0,30].
            species: List of species on each site. Can take in flexible input,
                including:

                i.  A sequence of element / specie specified either as string
                    symbols, e.g. ["Li", "Fe2+", "P", ...] or atomic numbers,
                    e.g., (3, 56, ...) or actual Element or Specie objects.

                ii. List of dict of elements/species and occupancies, e.g.,
                    [{"Fe" : 0.5, "Mn":0.5}, ...]. This allows the setup of
                    disordered structures.
            fractional_coords: list of fractional coordinates of each species.
            validate_proximity (bool): Whether to check if there are sites
                that are less than 0.01 Ang apart. Defaults to False.
            coords_are_cartesian (bool): Set to True if you are providing
                coordinates in cartesian coordinates. Defaults to False.
            site_properties (dict): Properties associated with the sites as a
                dict of sequences, e.g., {"magmom":[5,5,5,5]}. The sequences
                have to be the same length as the atomic species and
                fractional_coords. Defaults to None for no properties.
        """
        IStructure.__init__(
            self, lattice, species, coords,
            validate_proximity=validate_proximity, to_unit_cell=to_unit_cell,
            coords_are_cartesian=coords_are_cartesian,
            site_properties=site_properties,site_properties_T=site_properties_T,
            intensive_properties=intensive_properties,extensive_properties=extensive_properties)

        self._sites = list(self._sites)

    def set_extra(self, scale, R):
        self.lattice.set_R(scale, R)

    def __setitem__(self, i, site):
        """
        Modify a site in the structure.

        Args:
            i (int): Index
            site (PeriodicSite/Specie/Sequence): Three options exist. You
                can provide a PeriodicSite directly (lattice will be
                checked). Or more conveniently, you can provide a
                specie-like object or a tuple of up to length 3. Examples:
                s[0] = "Fe"
                s[0] = Element("Fe")
                both replaces the species only.
                s[0] = "Fe", [0.5, 0.5, 0.5]
                Replaces site and *fractional* coordinates. Any properties
                are inherited from current site.
                s[0] = "Fe", [0.5, 0.5, 0.5], {"spin": 2}
                Replaces site and *fractional* coordinates and properties.
        """
        if isinstance(site, PeriodicSite):
            if site.lattice != self._lattice:
                raise ValueError("PeriodicSite added must have same lattice "
                                 "as Structure!")
            self._sites[i] = site
        else:
            if isinstance(site, str) or (not isinstance(site, collections.Sequence)):
                sp = site
                frac_coords = self._sites[i].frac_coords
                properties = self._sites[i].properties
            else:
                sp = site[0]
                frac_coords = site[1] if len(site) > 1 else self._sites[i]\
                    .frac_coords
                properties = site[2] if len(site) > 2 else self._sites[i]\
                    .properties

            self._sites[i] = PeriodicSite(sp, frac_coords, self._lattice,
                                          properties=properties)

    def __delitem__(self, i):
        """
        Deletes a site from the Structure.
        """
        self._sites.__delitem__(i)

    def append(self, species, coords, coords_are_cartesian=False,
               validate_proximity=False, properties=None):
        """
        Append a site to the structure.

        Args:
            species: Species of inserted site
            coords (3x1 array): Coordinates of inserted site
            coords_are_cartesian (bool): Whether coordinates are cartesian.
                Defaults to False.
            validate_proximity (bool): Whether to check if inserted site is
                too close to an existing site. Defaults to False.

        Returns:
            New structure with inserted site.
        """
        return self.insert(len(self), species, coords,
                           coords_are_cartesian=coords_are_cartesian,
                           validate_proximity=validate_proximity,
                           properties=properties)

    def insert(self, i, species, coords, coords_are_cartesian=False,
               validate_proximity=False, properties=None):
        """
        Insert a site to the structure.

        Args:
            i (int): Index to insert site
            species (species-like): Species of inserted site
            coords (3x1 array): Coordinates of inserted site
            coords_are_cartesian (bool): Whether coordinates are cartesian.
                Defaults to False.
            validate_proximity (bool): Whether to check if inserted site is
                too close to an existing site. Defaults to False.

        Returns:
            New structure with inserted site.
        """
        if not coords_are_cartesian:
            new_site = PeriodicSite(species, coords, self._lattice,
                                    properties=properties)
        else:
            frac_coords = self._lattice.get_fractional_coords(coords)
            new_site = PeriodicSite(species, frac_coords, self._lattice,
                                    properties=properties)

        if validate_proximity:
            for site in self:
                if site.distance(new_site) < self.DISTANCE_TOLERANCE:
                    raise ValueError("New site is too close to an existing "
                                     "site!")

        self._sites.insert(i, new_site)


    def add_site_property(self, property_name, values):
        """
        Adds a property to all sites.

        Args:
            property_name (str): The name of the property to add.
            values: A sequence of values. Must be same length as number of
                sites.
        """
        if values is None:
            return
        if len(values) != len(self._sites):
            raise ValueError("Values must be same length as sites.")
        for i in range(len(self._sites)):
            site = self._sites[i]
            props = site.properties
            if not props:
                props = {}
            props[property_name] = values[i]
            self._sites[i] = PeriodicSite(site.species_and_occu,
                                          site.frac_coords, self._lattice,
                                          properties=props)

    def replace_species(self, species_mapping):
        """
        Swap species in a structure.

        Args:
            species_mapping (dict): Dict of species to swap. Species can be
                elements too. e.g., {Element("Li"): Element("Na")} performs
                a Li for Na substitution. The second species can be a
                sp_and_occu dict. For example, a site with 0.5 Si that is
                passed the mapping {Element('Si): {Element('Ge'):0.75,
                Element('C'):0.25} } will have .375 Ge and .125 C.
        """
        latt = self._lattice
        species_mapping = {get_el_sp(k): v
                           for k, v in species_mapping.items()}

        def mod_site(site):
            new_atom_occu = collections.defaultdict(int)
            for sp, amt in site.species_and_occu.items():
                if sp in species_mapping:
                    if isinstance(species_mapping[sp], collections.Mapping):
                        for new_sp, new_amt in species_mapping[sp].items():
                            new_atom_occu[get_el_sp(new_sp)] \
                                += amt * new_amt
                    else:
                        new_atom_occu[get_el_sp(
                            species_mapping[sp])] += amt
                else:
                    new_atom_occu[sp] += amt
            return PeriodicSite(new_atom_occu, site.frac_coords, latt,
                                properties=site.properties)

        self._sites = map(mod_site, self._sites)

    def replace(self, i, species, coords=None, coords_are_cartesian=False,
                properties=None):
        """
        Replace a single site. Takes either a species or a dict of species and
        occupations.

        Args:
            i (int): Index of the site in the _sites list.
            species (species-like): Species of replacement site
            coords (3x1 array): Coordinates of replacement site. If None,
                the current coordinates are assumed.
            coords_are_cartesian (bool): Whether coordinates are cartesian.
                Defaults to False.
            validate_proximity (bool): Whether to check if inserted site is
                too close to an existing site. Defaults to False.
        """
        if coords is None:
            frac_coords = self[i].frac_coords
        elif coords_are_cartesian:
            frac_coords = self._lattice.get_fractional_coords(coords)
        else:
            frac_coords = coords

        new_site = PeriodicSite(species, frac_coords, self._lattice,
                                properties=properties)
        self._sites[i] = new_site

    def remove_species(self, species):
        """
        Remove all occurrences of several species from a structure.

        Args:
            species: Sequence of species to remove, e.g., ["Li", "Na"].
        """
        new_sites = []
        species = map(get_el_sp, species)

        for site in self._sites:
            new_sp_occu = {sp: amt for sp, amt in site.species_and_occu.items()
                           if sp not in species}
            if len(new_sp_occu) > 0:
                new_sites.append(PeriodicSite(
                    new_sp_occu, site.frac_coords, self._lattice,
                    properties=site.properties))
        self._sites = new_sites

    def remove_sites(self, indices):
        """
        Delete sites with at indices.

        Args:
            indices: Sequence of indices of sites to delete.
        """
        self._sites = [self._sites[i] for i in range(len(self._sites))
                       if i not in indices]

    def apply_operation(self, symmop):
        """
        Apply a symmetry operation to the structure and return the new
        structure. The lattice is operated by the rotation matrix only.
        Coords are operated in full and then transformed to the new lattice.

        Args:
            symmop (SymmOp): Symmetry operation to apply.
        """
        self._lattice = Lattice([symmop.apply_rotation_only(row)
                                 for row in self._lattice.matrix])

        def operate_site(site):
            new_cart = symmop.operate(site.coords)
            new_frac = self._lattice.get_fractional_coords(new_cart)
            return PeriodicSite(site.species_and_occu, new_frac, self._lattice,
                                properties=site.properties)
        self._sites = map(operate_site, self._sites)

    def modify_lattice(self, new_lattice):
        """
        Modify the lattice of the structure.  Mainly used for changing the
        basis.

        Args:
            new_lattice (Lattice): New lattice
        """
        self._lattice = new_lattice
        new_sites = []
        for site in self._sites:
            new_sites.append(PeriodicSite(site.species_and_occu,
                                          site.frac_coords,
                                          self._lattice,
                                          properties=site.properties))
        self._sites = new_sites

    def apply_strain(self, strain):
        """
        Apply a strain to the lattice.

        Args:
            strain (float or list): Amount of strain to apply. Can be a float,
                or a sequence of 3 numbers. E.g., 0.01 means all lattice
                vectors are increased by 1%. This is equivalent to calling
                modify_lattice with a lattice with lattice parameters that
                are 1% larger.
        """
        s = (1 + np.array(strain)) * np.eye(3)
        self.modify_lattice(Lattice(np.dot(self._lattice.matrix.T, s).T))

    def translate_sites(self, indices, vector, frac_coords=True,
                        to_unit_cell=True):
        """
        Translate specific sites by some vector, keeping the sites within the
        unit cell.

        Args:
            indices: Integer or List of site indices on which to perform the
                translation.
            vector: Translation vector for sites.
            frac_coords (bool): Whether the vector corresponds to fractional or
                cartesian coordinates.
            to_unit_cell (bool): Whether new sites are transformed to unit
                cell
        """
        if not isinstance(indices, collections.Iterable):
            indices = [indices]

        for i in indices:
            site = self._sites[i]
            if frac_coords:
                fcoords = site.frac_coords + vector
            else:
                fcoords = self._lattice.get_fractional_coords(site.coords
                                                              + vector)
            new_site = PeriodicSite(site.species_and_occu, fcoords,
                                    self._lattice, to_unit_cell=to_unit_cell,
                                    coords_are_cartesian=False,
                                    properties=site.properties)
            self._sites[i] = new_site

    def perturb(self, distance, to_unit_cell=True):
        """
        Performs a random perturbation of the sites in a structure to break
        symmetries.

        Args:
            distance (float): Distance in angstroms by which to perturb each
                site.
        """
        def get_rand_vec():
            #deals with zero vectors.
            vector = np.random.randn(3)
            vnorm = np.linalg.norm(vector)
            return vector / vnorm * distance if vnorm != 0 else get_rand_vec()

        for i in range(len(self._sites)):
            self.translate_sites([i], get_rand_vec(), frac_coords=False,to_unit_cell=to_unit_cell)


    def displace_by(self, dr, to_unit_cell=True):
        for i in range(len(self._sites)):
            self.translate_sites([i], dr[i], frac_coords=False,to_unit_cell=to_unit_cell)

    def from_displacement(self, dr, to_unit_cell=True):
        s1= Structure(self._lattice, self.species_and_occu, self.frac_coords)
        s1.displace_by(dr, to_unit_cell=to_unit_cell)
        return s1

    def add_oxidation_state_by_element(self, oxidation_states):
        """
        Add oxidation states to a structure.

        Args:
            oxidation_states (dict): Dict of oxidation states.
                E.g., {"Li":1, "Fe":2, "P":5, "O":-2}
        """
        try:
            for i, site in enumerate(self._sites):
                new_sp = {}
                for el, occu in site.species_and_occu.items():
                    sym = el.symbol
                    new_sp[Specie(sym, oxidation_states[sym])] = occu
                new_site = PeriodicSite(new_sp, site.frac_coords,
                                        self._lattice,
                                        coords_are_cartesian=False,
                                        properties=site.properties)
                self._sites[i] = new_site

        except KeyError:
            raise ValueError("Oxidation state of all elements must be "
                             "specified in the dictionary.")

    def add_oxidation_state_by_site(self, oxidation_states):
        """
        Add oxidation states to a structure by site.

        Args:
            oxidation_states (list): List of oxidation states.
                E.g., [1, 1, 1, 1, 2, 2, 2, 2, 5, 5, 5, 5, -2, -2, -2, -2]
        """
        try:
            for i, site in enumerate(self._sites):
                new_sp = {}
                for el, occu in site.species_and_occu.items():
                    sym = el.symbol
                    new_sp[Specie(sym, oxidation_states[i])] = occu
                new_site = PeriodicSite(new_sp, site.frac_coords,
                                        self._lattice,
                                        coords_are_cartesian=False,
                                        properties=site.properties)
                self._sites[i] = new_site

        except IndexError:
            raise ValueError("Oxidation state of all sites must be "
                             "specified in the dictionary.")

    def remove_oxidation_states(self):
        """
        Removes oxidation states from a structure.
        """
        for i, site in enumerate(self._sites):
            new_sp = collections.defaultdict(float)
            for el, occu in site.species_and_occu.items():
                sym = el.symbol
                new_sp[Element(sym)] += occu
            new_site = PeriodicSite(new_sp, site.frac_coords,
                                    self._lattice,
                                    coords_are_cartesian=False,
                                    properties=site.properties)
            self._sites[i] = new_site

    def generate_supercell(self, scaling_matrix, scref=None):
        """
        Create a supercell.

        Args:
            scaling_matrix: A scaling matrix for transforming the lattice
                vectors. Has to be all integers. Several options are possible:

                a. A full 3x3 scaling matrix defining the linear combination
                   the old lattice vectors. E.g., [[2,1,0],[0,3,0],[0,0,
                   1]] generates a new structure with lattice vectors a' =
                   2a + b, b' = 3b, c' = c where a, b, and c are the lattice
                   vectors of the original structure.
                b. An sequence of three scaling factors. E.g., [2, 1, 1]
                   specifies that the supercell should have dimensions 2a x b x
                   c.
                c. A number, which simply scales all lattice vectors by the
                   same factor.
        """
        scmat = np.array(scaling_matrix, np.int16)
        if scmat.shape != (3, 3):
            scmat= np.array(scmat* np.eye(3), np.int16)
        n_cell=int(round(np.linalg.det(scmat)))
        old_lattice = self._lattice
        new_lattice = Lattice(np.dot(scmat, old_lattice.matrix))
        tvects = supercell_latticepoints(scmat)
        inv=np.linalg.inv(scmat)
        if scref is None:
            sc_ref= supercell_latticepoints(scmat)
        else:
            sc_ref= scref
        return Structure(Lattice(np.dot(scmat, self._lattice.matrix)),
           [s.species_and_occu for s in self for _ in range(n_cell)], 
           (self.frac_coords[:,None,:]+sc_ref[None,:,:]).reshape((-1,3)).dot(inv),
           coords_are_cartesian=False, to_unit_cell=True,
           site_properties_T=[s.properties for s in self for _ in range(n_cell)],
           intensive_properties=self.intensive_properties,extensive_properties=
                           {k:v*self.n_cell for k,v in self.extensive_properties.items()})

    def optimize_supercell(nsc1, maxIter=2000):
        """
        search for optimal supercell shape (as cubic like as possible)
        Args:
            nsc1: positive means number of supercells targeted
                 negative means a certain number of atoms desired
            maxIter: number of iterations
        """
        nsc=nsc1
        if nsc<0:
            nsc=int(round(nsc/self.num_sites))
        volsc = nsc*self.volume
        invprim = np.linalg.inv(self.lattice)
        ntry=0
        bestsc=np.identity(3, dtype=int)
        bestlen=999999.0
        for i in range(maxIter):
            scmat=1



    def scale_lattice(self, volume):
        """
        Performs a scaling of the lattice vectors so that length proportions
        and angles are preserved.

        Args:
            volume (float): New volume of the unit cell in A^3.
        """
        self.modify_lattice(self._lattice.scale(volume))


    def ijkl2frac(self, ijkl):
        """
        same but ijkl in one array
        """
        return np.array(ijkl[:3]) + self._sites[int(ijkl[3])].frac_coords


    def ijkl2cart(self, ijkl):
        """
        :return: cartesian coordinates
        """
        return np.dot(self.ijkl2frac(ijkl), self.lattice._matrix)

    def frac2ijkl(self, coords, frac_coords=True, tolerance= 1E-3):
        """
        Identify which atom corresponds to the specified coordinates

        Args:
            coords (3x1 array): Array-like object with the coordinates.
            frac_coords (bool): Whether the vector corresponds to fractional or
                cartesian coordinates.

        Returns:
            integer index for the identified site
        """
        assert frac_coords == True
        return Structure.frac2ijkl_fromapos(coords, self.frac_coords, tolerance)

    @staticmethod
    def frac2ijkl_fromapos(c_f, apos, tolerance):
        for l in range(len(apos)):
            rvec= c_f - apos[l]
            rvec_frac = rvec - np.round(rvec)
            if np.linalg.norm(rvec_frac) < tolerance:
                return np.append(np.round(rvec), [l]).astype(int)
        print('debug from apos', c_f)
        raise ValueError("could not find [%f, %f, %f] in cell" % (c_f[0], c_f[1], c_f[2]))

    @staticmethod
    def ijkl_in_supercell(scmat, ijkl):
        """

        :param scmat: supercell scaling matrix
        :param refpts: list of lattice points within the supercell
        :return: new ijkl in the supercell.
        Asuming ijkl.inv(scmat) give new ijk, and refpts should be within the supercell spanned by scmat
        """
        inv = np.linalg.inv(scmat)
        nsc = abs(int(round(np.linalg.det(scmat))))
        return _ijkl_in_supercell(nsc, inv, supercell_latticepoints(scmat).dot(inv), ijkl)

    @staticmethod
    def _ijkl_in_supercell(nsc, invscmat, refpts, ijkl):
        newfrac= np.dot(ijkl[:3], invscmat)
        newijkl = Structure.frac2ijkl_fromapos(newfrac, refpts, 1E-3)
        newijkl[3] += ijkl[3]*nsc
        return newijkl


    def nbtable(self, r):
        """
        r: cutoff distance
        return a table: neighbor list of each atom
        """
        if not hasattr(self, '_nbtable'):
            self._nbtable = {}
        if r in self._nbtable:
            return self._nbtable[r]
        recp_len = np.array(self.lattice.reciprocal_lattice.abc)
        sr = r + 0.15
        nmax = sr * recp_len / (2 * math.pi)
        floor = math.floor
        n = self.num_sites
        fcoords = self.frac_coords
        indices = np.array(range(n))
        nbtable_per_atom = []
        pmin = np.amin(fcoords, axis=0)
        pmax = np.amax(fcoords, axis=0)

        arange = np.arange(int(floor(pmin[0] - nmax[0])),
                           int(floor(pmax[0] + nmax[0])) + 1)
        brange = np.arange(int(floor(pmin[1] - nmax[1])),
                           int(floor(pmax[1] + nmax[1])) + 1)
        crange = np.arange(int(floor(pmin[2] - nmax[2])),
                           int(floor(pmax[2] + nmax[2])) + 1)
        # print("debug arange=", arange.shape)

        arange = arange[:, None] * np.array([1, 0, 0])[None, :]
        brange = brange[:, None] * np.array([0, 1, 0])[None, :]
        crange = crange[:, None] * np.array([0, 0, 1])[None, :]
        # print("debug arange=", arange.shape, arange)

        images = arange[:, None, None] + brange[None, :, None] + crange[None, None, :]
        images = images.reshape((-1,3))
        shifted_coords = fcoords[:, None, :] + images[None, :, :]
        shifted_coords= shifted_coords.reshape((-1,3))
        coords = self.lattice.get_cartesian_coords(shifted_coords)
        ijkls = np.array([[img[0], img[1], img[2], l] for l in range(n) for img in images])
        for i in range(n):
            pct = self.cart_coords[i]
            dists = np.array([np.sqrt(np.sum((p-pct) ** 2)) for p in coords])
            within_r = np.where(dists <= r)
            nbtable_per_atom.append(ijkls[within_r])
        self._nbtable[r] = nbtable_per_atom
        return nbtable_per_atom

    def find_nb_cluster(self, ijkls, cut):
        """
        cut: cutoff distance
        Find atoms within cutoff of EVERY ijkl
        """
        # print(ijkls)
        nbtable = self.nbtable(cut)
        nsite = ijkls.shape[0]
        # print("nsite", nsite)
        if nsite <=0:
            raise ValueError("At least 1 atom in cluster needed")

        # print("ijkls", ijkls)
        atoms = []
        for _atom in nbtable[ijkls[0][3]]:
            atom = _atom.copy()
            atom[:3] += ijkls[0][:3]
            # print("testing", atom)
            within = True
            for j in range(1, nsite):
                atom_wrt_j = atom.copy()
                atom_wrt_j[:3] -= ijkls[j,:3]
                within = False
                for x in nbtable[ijkls[j][3]]:
                    if (atom_wrt_j == x).all():
                        within = True
                        break
                if not within:
                    # print("atom_wrt_j", atom_wrt_j, "not in", ijkls[j,:])
                    break
            if within:
                atoms.append(atom)

        return atoms

    def get_scmat(self, sc_R):
        """
        Given primitive cell (self) and supercell lattice vectors
        :param sc_R: lattice vectors of supercell
        :return:     integer scaling matrix
        """
        return np.dot(sc_R.lattice.matrix if isinstance(sc_R, Structure) else sc_R, self.lattice.inv_matrix).round().astype(int)


    def map_to_prim(self, prim):
        """
        Given supercell (self) and primitive cell, get supercell without displacement
        :param prim: primitive
        :return:     supercell without displacement
        """
        scmat = prim.get_scmat(self.lattice.matrix)
        sc=prim.generate_supercell(scmat)
        return self.map_to_reference(sc)


    def map_to_reference(self, sc):
        """
        Given supercell (self) and ideal supercell, get supercell without displacement
        Assume that
        :param sc: supercell
        :return:     supercell without displacement
        """
        assert self.num_sites == sc.num_sites
        dist_mat = self.lattice.get_all_distances(self.frac_coords, sc.frac_coords)
        jmin = np.argmin(dist_mat, axis=1)
        bad_i, bad_j = non_1to1(jmin)
        if len(bad_i) > 0:
            print("** WARNING** found %d conflicting mappings" % (len(bad_i)))
            print(self.frac_coords[bad_i], "==>", sc.frac_coords[bad_j])
            # try to resolve conflict
            from itertools import permutations
            min_dist = 1E99
            solve_j = bad_j
            for try_j in permutations(bad_j):
                dist = np.sum([dist_mat[bad_i[i], try_j[i]] for i in range(len(bad_i))])
                if dist < min_dist:
                    min_dist = dist
                    solve_j = [try_j]
            jmin[bad_i] = solve_j
            print(bad_j, solve_j)
            print("resolved", self.frac_coords[bad_i], "==>", sc.frac_coords[solve_j])

        for i in range(self.num_sites):
            # print("%d %.4f" % (i, np.linalg.norm(self.frac_coords[i] - sc.frac_coords[jmin[i]])))
            self[i].set_coords(np.round(self.frac_coords[i] - sc.frac_coords[jmin[i]]) + sc.frac_coords[jmin[i]], cart=False)
            # self[i].set_coords(sc.frac_coords[jmin[i]], cart=False)
            # print("%d %.4f" % (i, np.linalg.norm(self.frac_coords[i] - sc.frac_coords[jmin[i]])))
        return self


    def set_coords(self, c, cart=True):
#        print('debug set_coords', c.shape)
        for i in range(self.num_sites):
            self[i].set_coords(c[i], cart=cart)

    def get_order_wrt(self, p1, inverse=False, tol=1E-4):
        from _c_util import get_structure_ordering
        if p1 is None:
            return list(range(self.num_sites))

        if isinstance(p1, Structure):
            assert (np.abs(self.lattice._matrix-p1.lattice._matrix)<1E-6).all(), "ERROR difference lattice"
        pos= p1.frac_coords if isinstance(p1, Structure) else p1
        if inverse:
            ordering= get_structure_ordering(pos, self.frac_coords, 1, tol).tolist()
        else:
            ordering= get_structure_ordering(self.frac_coords, pos, 1, tol).tolist()

        # if isinstance(p1, Structure):
        #     inSC = p1
        # else:
        #     inSC = Structure(self.lattice, ['']*self.num_sites, p1)
        # ordering= [inSC.frac2ijkl(pf)[3] for pf in self.frac_coords]
        assert sorted(ordering) == list(range(self.num_sites))
        return ordering

    def to_spglib(self):
        """
        To the 'cell' format (tuple) required by spglib
        """
        unique_species = []
        zs = []
        magmoms = []

        for species, g in itertools.groupby(self,
                                            key=lambda s: s.species_and_occu):
            if species in unique_species:
                ind = unique_species.index(species)
                zs.extend([ind + 1] * len(tuple(g)))
            else:
                unique_species.append(species)
                zs.extend([len(unique_species)] * len(tuple(g)))

        for site in self:
            if hasattr(site, 'magmom'):
                magmoms.append(site.magmom)
            elif site.is_ordered and hasattr(site.specie, 'spin'):
                magmoms.append(site.specie.spin)
            else:
                magmoms.append(0)

        return self.lattice._matrix, self.frac_coords, zs, magmoms

    def from_spglib(self, latt, fcoord, id):
        spe= self.types_of_specie
        return Structure(latt, [spe[i-1] for i in id], fcoord)

    def standardize_cell(self, to_primitive=False, no_idealize=False, symprec=1e-5):
        """
        call standardize_cell of spglib
        """
        import spglib
        cell=self.to_spglib()
        return self.from_spglib(*spglib.standardize_cell(cell, to_primitive, no_idealize, symprec))

    def match(self, p2, tol_match=0.15, tol_distinct=1.0):
        """
        self: ideal (supercell) structure
        p2: structure with distortion, defect and/or disorder
        """
        from f_util import f_util
        return f_util.match_structure(self.lattice.matrix.T, self.frac_coords.T, self.atomic_numbers, p2.frac_coords.T, p2.atomic_numbers, tol_match, tol_distinct)


class SupercellStructure(Structure):
    """
    This class represents a Supercell related to SymmetrizedStructure by sc_mat matrix

    Args:
        prim: primitive cell
        sc_mat: defining supercell lattice vectors as R_sc = sc_mat . R_prim

    .. attribute: equivalent_indices

        indices of structure grouped by equivalency
    """

    def __init__(self, prim, sc_mat, strc, tol=1E-4, match_prim=True):
        self.sc_mat = np.array(sc_mat, dtype=int)
        n_cell = abs(np.linalg.det(self.sc_mat))
        assert abs(n_cell - np.round(n_cell).astype(int))< 1E-10, "ncell= %f "%(n_cell)
        self.n_cell = np.round(n_cell).astype(int)
        self.prim = prim
        self.inv_sc_mat = np.linalg.inv(sc_mat)

        if True or match_prim:
            strc.map_to_prim(prim)
            strc.perturb(0, to_unit_cell=True) # move to_unit_cell
            ijkl= np.array([prim.frac2ijkl(c) for c in strc.frac_coords.dot(sc_mat)])
            self.ijk_ref= np.array([i[:3] for i in ijkl if i[3]==0])
            self._ijkl = ijkl
            self._ijkl_list = ijkl.tolist()
            species=[prim[i[3]].species_and_occu for i in ijkl]
            site_properties_T=[prim[i[3]].properties for i in ijkl]
        else:
            species=[s.species_and_occu for s in strc]
            site_properties_T=[s.properties for s in strc]
        self.sc_ref=self.ijk_ref
        Structure.__init__(self, Lattice(np.dot(self.sc_mat, prim.lattice._matrix)), species,
                strc.frac_coords, site_properties_T=site_properties_T,
                intensive_properties=prim.intensive_properties,extensive_properties=
                           {k:v*self.n_cell for k,v in prim.extensive_properties.items()})

    @classmethod
    def from_scmat(cls, prim, scmat, scref= None):
        return cls(prim, scmat, prim.generate_supercell(scmat, scref))
        n_cell = abs(np.linalg.det(scmat))
        assert abs(n_cell - np.round(n_cell).astype(int))< 1E-10, "ncell= %f "%(n_cell)
        n_cell=int(n_cell)
        newlattice = Lattice(np.dot(scmat, prim.lattice.matrix))
        if scref is None:
            sc_ref= supercell_latticepoints(scmat)
        else:
            sc_ref= scref
        assert n_cell == len(sc_ref)
        new_sites = []
        for site in prim:
            for t in sc_ref:
                fcoords = site.frac_coords + t
                coords = prim.lattice.get_cartesian_coords(fcoords)
                new_site = PeriodicSite(
                    site.species_and_occu, coords, newlattice,
                    coords_are_cartesian=True, properties=site.properties,
                    to_unit_cell=True)
                new_sites.append(new_site)
        strc=Structure.__init__(newlattice, [s.species_and_occu for s in prim for j in range(n_cell)],
                newfrac, site_properties_T=[s.properties for s in prim for j in range(n_cell)],
                intensive_properties=prim.intensive_properties,extensive_properties=
                           {k:v*n_cell for k,v in prim.extensive_properties.items()})
        return cls(prim, scmat, strc)


    @classmethod
    def from_structure(cls, prim, strc):
        return cls(prim, prim.get_scmat(strc), strc)


    @classmethod
    def from_file(cls, prim, f):
        from .interface_vasp import Poscar
        strc= Poscar.from_file(f).structure
        return cls(prim, prim.get_scmat(strc), strc)    


    def compatible_kpoints(self):
        """
        return K-points compatible with this supercell
        """
        return np.dot(supercell_latticepoints(self.sc_mat.T), self.inv_sc_mat.T)


    def to_unperturbed(self):
        """
        Map positions to unperturbed (prim) positions by the best match,
        i.e. the --task 2 of polaron_main
        """
        self.map_to_prim(self.prim)
        

class StructureError(Exception):
    """
    Exception class for Structure.
    Raised when the structure has problems, e.g., atoms that are too close.
    """
    pass
