from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

__author__ = 'Ehud Tsivion'

from molecule import Molecule
from verbosity import verbprint


class Complex(object):
    """
    A Complex is a collection of monomers

    """

    def __init__(self, name="default", verbosity=1):
        self.name = name
        self.verbosity = 1
        self.monomer_list = []

    def from_molecule(self, molecule, size_list):
        """
        split a Molecule object into a Complex object
        (which contains several Molecule object,
        representing monomers)

        :param: size_list: a list of the number of atoms of
        each of the monomers.

        example:

        size_list = [5, 6, 1, ...] means:

        monomer 0 is the first five atoms
        monomer 1 is the next 6 atoms
        monomer 2 is the next 1 atom
        and so on

        :return: A complex object
        """

        # count the atoms to make sure all atoms
        # are assigned to a fragment
        atom_counter = 0
        i = 0

        for frag_size in size_list:

            # I don't know why, but the "atom_list=[]"
            # is essential
            j = i + frag_size

            monomer = Molecule(atoms=molecule.atoms[i:j])

            i = j

            self.add_monomer(monomer)

            atom_counter += frag_size

        if atom_counter == molecule.atom_count:
            pass
        else:
            raise Exception('Failed to parse molecule to monomers: '
                            'specified {}, but there are {} atoms'
                            .format(atom_counter, molecule.atom_count))

        verbprint(2, self.verbosity, '\nMonomer list:')
        verbprint(2, self.verbosity, self)

        return self

    def add_monomer(self, monomer_molecule):
        self.monomer_list.append(monomer_molecule)

    def __str__(self):
        counter = 0
        text = "\nComplex Object:\n"
        for momoner in self.monomer_list:
            counter += 1
            text += "Monomer number {}:\n".format(counter)
            text += str(momoner)
        return text








