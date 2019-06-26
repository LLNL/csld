#!/usr/bin/env python3

"""
This module provides base class used to define a collection of atoms that 
represent a certain interaction, most notably a cluster in a
crystalline lattice, e.g.
  * proper cluster as in cluster expansion, Ising/Heisenberg model
  * improper/proper cluster, in lattice dynamics
"""


from abc import ABCMeta, abstractmethod, abstractproperty


class AtomicModel():
    """
    A model such as cluster expansion, lattice dynamics
    """
    __metaclass__ = ABCMeta

    @abstractproperty
    def n_raw_parameters(self):
        """
        Returns the number of ALL parameters (except lattice translation).
        """
        return

    @abstractproperty
    def n_independent_parameters(self):
        """
        Returns the number of independent parameters.
        """
        return

    @abstractmethod
    def evaluate(self, structure):
        """
         return the property given the input structure
        """
        return
