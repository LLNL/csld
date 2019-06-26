#!/usr/bin/env python3

"""
Classes for reading POSCAR, OUTCAR
"""

# adapted from pymatgen.io_util


import os
import re
import itertools
import warnings
import logging

import numpy as np
from numpy.linalg import det


from .structure import Structure
from .lattice import Lattice
from .util.periodic_table import Element
from .util.io_utils import clean_lines, zopen


logger = logging.getLogger(__name__)


class Poscar():
    """
    Object for representing the data in a POSCAR or CONTCAR file.
    Please note that this current implementation. Most attributes can be set
    directly.

    Args:
        structure (Structure):  Structure object.
        comment (str): Optional comment line for POSCAR. Defaults to unit
            cell formula of structure. Defaults to None.
        selective_dynamics (Nx3 array): bool values for selective dynamics,
            where N is number of sites. Defaults to None.
        true_names (bool): Set to False is the names in the POSCAR are not
            well-defined and ambiguous. This situation arises commonly in
            vasp < 5 where the POSCAR sometimes does not contain element
            symbols. Defaults to True.
        velocities (Nx3 array): Velocities for the POSCAR. Typically parsed
            in MD runs or can be used to initialize velocities.
        predictor_corrector (Nx3 array): Predictor corrector for the POSCAR.
            Typically parsed in MD runs.

    .. attribute:: structure

        Associated Structure.

    .. attribute:: comment

        Optional comment string.

    .. attribute:: true_names

        Boolean indication whether Poscar contains actual real names parsed
        from either a POTCAR or the POSCAR itself.

    .. attribute:: selective_dynamics

        Selective dynamics attribute for each site if available. A Nx3 array of
        booleans.

    .. attribute:: velocities

        Velocities for each site (typically read in from a CONTCAR). A Nx3
        array of floats.

    .. attribute:: predictor_corrector

        Predictor corrector coordinates for each site (typically read in from a
        MD CONTCAR).

    .. attribute:: temperature

        Temperature of velocity Maxwell-Boltzmann initialization. Initialized
        to -1 (MB hasn"t been performed).
    """

    def __init__(self, structure, comment=None,
                 true_names=True):
        self.structure = structure
        self.true_names = true_names
        self.comment = structure.formula if comment is None else comment

    @property
    def site_symbols(self):
        """
        Sequence of symbols associated with the Poscar. Similar to 6th line in
        vasp 5+ POSCAR.
        """
        syms = [site.specie.symbol for site in self.structure]
        return [a[0] for a in itertools.groupby(syms)]

    @property
    def natoms(self):
        """
        Sequence of number of sites of each type associated with the Poscar.
        Similar to 7th line in vasp 5+ POSCAR or the 6th line in vasp 4 POSCAR.
        """
        syms = [site.specie.symbol for site in self.structure]
        return [len(tuple(a[1])) for a in itertools.groupby(syms)]

    def __setattr__(self, name, value):
        super(Poscar, self).__setattr__(name, value)

    @staticmethod
    def from_file(filename, check_for_POTCAR= False, read_CE=False, warn_vasp4=True):
        """
        Reads a Poscar from a file.

        The code will try its best to determine the elements in the POSCAR in
        the following order:
        1. If check_for_POTCAR is True, the code will try to check if a POTCAR
        is in the same directory as the POSCAR and use elements from that by
        default. (This is the VASP default sequence of priority).
        2. If the input file is Vasp5-like and contains element symbols in the
        6th line, the code will use that if check_for_POTCAR is False or there
        is no POTCAR found.
        3. Failing (2), the code will check if a symbol is provided at the end
        of each coordinate.

        If all else fails, the code will just assign the first n elements in
        increasing atomic number, where n is the number of species, to the
        Poscar. For example, H, He, Li, ....  This will ensure at least a
        unique element is assigned to each site and any analysis that does not
        require specific elemental properties should work fine.

        Args:
            filename (str): File name containing Poscar data.
            check_for_POTCAR (bool): Whether to check if a POTCAR is present
                in the same directory as the POSCAR. Defaults to True.

        Returns:
            Poscar object.
        """
        names = None
        with zopen(filename, "r") as f:
            return Poscar.from_string(f.read(), read_CE=read_CE, warn_vasp4=warn_vasp4)

    @staticmethod
    def from_string(data, default_names=None, read_CE=False, warn_vasp4=True):
        """
        Reads a Poscar from a string.

        The code will try its best to determine the elements in the POSCAR in
        the following order:
        1. If default_names are supplied and valid, it will use those. Usually,
        default names comes from an external source, such as a POTCAR in the
        same directory.
        2. If there are no valid default names but the input file is Vasp5-like
        and contains element symbols in the 6th line, the code will use that.
        3. Failing (2), the code will check if a symbol is provided at the end
        of each coordinate.

        If all else fails, the code will just assign the first n elements in
        increasing atomic number, where n is the number of species, to the
        Poscar. For example, H, He, Li, ....  This will ensure at least a
        unique element is assigned to each site and any analysis that does not
        require specific elemental properties should work fine.

        Args:
            data (str): String containing Poscar data.
            default_names ([str]): Default symbols for the POSCAR file,
                usually coming from a POTCAR in the same directory.

        Returns:
            Poscar object.
        """
        #"^\s*$" doesn't match lines with no whitespace
        chunks = re.split("\n\s*\n", data.rstrip(), flags=re.MULTILINE)
        if chunks[0] == "":
            chunks.pop(0)
            chunks[0] = "\n" + chunks[0]
        #Parse positions
        lines = tuple(clean_lines(chunks[0].split("\n"), False))
        comment = lines[0]
        scale = float(lines[1])
        R = np.array([list(map(float, line.split()))
                            for line in lines[2:5]])
        if scale < 0:
            # In vasp, a negative scale factor is treated as a volume. We need
            # to translate this to a proper lattice vector scaling.
            vol = abs(det(R))
            scale= (-scale / vol) ** (1 / 3)

        lattice = R* scale

        vasp5_symbols = False
        try:
            natoms = list(map(int, lines[5].split()))
            ipos = 6
        except ValueError:
            vasp5_symbols = True
            symbols = lines[5].split()
            natoms = list(map(int, lines[6].split()))
            atomic_symbols = list()
            for i in range(len(natoms)):
                atomic_symbols.extend([symbols[i]] * natoms[i])
            ipos = 7

        postype = lines[ipos].split()[0]

        sdynamics = False
        # Selective dynamics
        if postype[0] in "sS":
            sdynamics = True
            ipos += 1
            postype = lines[ipos].split()[0]

        cart = postype[0] in "cCkK"
        nsites = sum(natoms)

        # If default_names is specified (usually coming from a POTCAR), use
        # them. This is in line with Vasp"s parsing order that the POTCAR
        # specified is the default used.
        if default_names:
            try:
                atomic_symbols = []
                for i in range(len(natoms)):
                    atomic_symbols.extend([default_names[i]] * natoms[i])
                vasp5_symbols = True
            except IndexError:
                pass

        if not vasp5_symbols:
            ind = 3 if not sdynamics else 6
            try:
                #check if names are appended at the end of the coordinates.
                atomic_symbols = [l.split()[ind]
                                  for l in lines[ipos + 1:ipos + 1 + nsites]]
                #Ensure symbols are valid elements
                if not all([Element.is_valid_symbol(sym)
                            for sym in atomic_symbols]):
                    raise ValueError("Non-valid symbols detected.")
                vasp5_symbols = True
            except (ValueError, IndexError):
                #Defaulting to false names.
                atomic_symbols = []
                for i in range(len(natoms)):
                    sym = Element.from_Z(i + 1).symbol
                    atomic_symbols.extend([sym] * natoms[i])
                if warn_vasp4:
                    warnings.warn("Elements in POSCAR cannot be determined. "
                              "Defaulting to false names {}."
                              .format(" ".join(atomic_symbols)))

        # read the atomic coordinates
        coords = []
        selective_dynamics = list() if sdynamics else None
        alloy_species= [] if read_CE else None

        for i in range(nsites):
            toks = lines[ipos + 1 + i].split()
            crd_scale = scale if cart else 1
            coords.append([float(j)*crd_scale for j in toks[:3]])
            if sdynamics:
                selective_dynamics.append([tok.upper()[0] == "T"
                                           for tok in toks[3:6]])
            if read_CE:
                alloy_species.append(' '.join(toks[3:]))


#        print("debug> lattice vec=", atomic_symbols, Z_atom, cart)
        struct = Structure(lattice, atomic_symbols, coords, False, False, cart)
        struct.set_extra(scale, R)
        if read_CE:
            struct.add_site_property("alloy_species", alloy_species)

        return Poscar(struct, comment, vasp5_symbols)

    def get_string(self, direct=True, vasp4_compatible=False,
                   significant_figures=16):
        """
        Returns a string to be written as a POSCAR file. By default, site
        symbols are written, which means compatibility is for vasp >= 5.

        Args:
            direct (bool): Whether coordinates are output in direct or
                cartesian. Defaults to True.
            vasp4_compatible (bool): Set to True to omit site symbols on 6th
                line to maintain backward vasp 4.x compatibility. Defaults
                to False.
            significant_figures (int): No. of significant figures to
                output all quantities. Defaults to 6. Note that positions are
                output in fixed point, while velocities are output in
                scientific format.

        Returns:
            String representation of POSCAR.
        """

        # This corrects for VASP really annoying bug of crashing on lattices
        # which have triple product < 0. We will just invert the lattice
        # vectors.
        latt = self.structure.lattice
        if np.linalg.det(latt.matrix) < 0:
            latt = Lattice(-latt.matrix)

        lines = [self.comment, "1.0", str(latt)]
        if self.true_names and not vasp4_compatible:
            lines.append(" ".join(self.site_symbols))
        lines.append(" ".join([str(x) for x in self.natoms]))
#        if self.selective_dynamics:
#            lines.append("Selective dynamics")
        lines.append("direct" if direct else "cartesian")

        format_str = "{{:.{0}f}}".format(significant_figures)

        for (i, site) in enumerate(self.structure):
            coords = site.frac_coords if direct else site.coords
            line = " ".join([format_str.format(c) for c in coords])
            # if self.selective_dynamics is not None:
            #     sd = ["T" if j else "F" for j in self.selective_dynamics[i]]
            #     line += " %s %s %s" % (sd[0], sd[1], sd[2])
            line += " " + site.species_string
            lines.append(line)


        return "\n".join(lines) + "\n"

    def __str__(self):
        """
        String representation of Poscar file.
        """
        return self.get_string()

    def write_file(self, filename, **kwargs):
        """
        Writes POSCAR to a file. The supported kwargs are the same as those for
        the Poscar.get_string method and are passed through directly.
        """
        if isinstance(filename, str):
            with open(filename, "w") as f:
                f.write(self.get_string(**kwargs))
        else:
            filename.write(self.get_string(**kwargs))

    @property
    def to_dict(self):
        return {"@module": self.__class__.__module__,
                "@class": self.__class__.__name__,
                "structure": self.structure.to_dict,
                "true_names": self.true_names,
                "comment": self.comment}

    @classmethod
    def from_dict(cls, d):
        return Poscar(Structure.from_dict(d["structure"]),
                      comment=d["comment"],
                      true_names=d["true_names"])
