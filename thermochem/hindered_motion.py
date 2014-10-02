from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

__author__ = 'Ehud Tsivion'

from verbosity import verbprint
from molden_parser import MoldenIO
from complex import Complex


class HNMA(object):
    """
    A class for managing the "Hindered motions analysis"
    which corrects of hindered motions that are lost in
    the standard process of extraction of thermodynamical
    properties from a frequency analysis.

    how the fragments are organized?

    fragments_list - is a list of fragments: [[frag1], [frag2], ..]
    each fragment, is a list of atoms that belong to that fragment: [1, 4, 5]

    """

    def __init__(self, verbosity=1, molecule=None):
        self.verbosity = verbosity
        self.molecule = molecule
        self.complex = None
        self.normal_nodes = None

    def parse_molden(self, molden_file, molden_job=1):
        """
        obtain relevant data from molden file

        :return:
        """

        verbprint(1, self.verbosity, 'Reading information from Molden Format output')

        # open the file, and automatically parse data
        mf = MoldenIO(molden_file, molden_job, self.verbosity)

        # get references for needed information
        self.molecule = mf._atoms
        self.normal_nodes = mf._normal_modes

    def get_thermodynamics(self, temperature=297.15, pressure=1):
        """

        :param temperature:
        :param pressure:
        :return:
        """
        verbprint(1, self.verbosity, '\n{:*^40}'.format(' Vibrational analysis '))

        # if fragments not set, set all the atoms on a single fragment
        if not self.monomer_list:
            verbprint(0, self.verbosity, "Fragments are not set,"
                                         " this is just a standard vibrational"
                                         " analysis")

            verbprint(0, self.verbosity, "setting molecule as single fragment")
            self.fragments_number = 1
            self.monomer_list = [range(self.atoms.atom_count)]

        for mode in self.normal_modes.modes_list:
            for frag in self.monomer_list:
                pass

    def set_complex(self, atom_list):

        if self.molecule:
            self.complex = Complex().from_molecule(self.molecule, atom_list)

        else:

            print("There is no Molecule object. Cannot split to complex")
            return None


if __name__ == "__main__":
    job = HNMA(verbosity=1)
    job.parse_molden('../example_outputs/catechol-Mg-3Me-opt-freq.qchem', molden_job=2)
    job.set_complex([24, 4])
    print(job.normal_nodes)