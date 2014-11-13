from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

__author__ = 'Ehud Tsivion'

from verbosity import verbprint
from molden_parser import MoldenIO
from molecule import Molecule
import numpy as np


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

    def __init__(self, verbosity=1, complex=None):
        self.verbosity = verbosity
        self.molecule = complex
        self.complex = complex
        self.adsorbate = None
        self.sorbent = None


    def parse_molden(self, molden_file, molden_job=1):
        """
        obtain relevant data from molden file

        :return:
        """

        verbprint(1, self.verbosity, 'Reading information from Molden Format output')

        # open the file, and automatically parse data
        mf = MoldenIO(molden_file, molden_job, self.verbosity)

        # get references for needed information
        self.molecule = mf.atoms
        self.molecule._normal_modes = mf.normal_modes

    def get_thermodynamics(self, temperature=297.15, pressure=1):
        """

        :param temperature:
        :param pressure:
        :return:
        """

        entropy = 0
        internal_energy = 0
        zero_point_energy = 0
        ads_vib_count = 0


        verbprint(1, self.verbosity, '\n{:*^50}'.format('Starting Hindered Vibrational analysis'))

        if not self.adsorbate:
            raise Exception("Adsorbate molecule is not defined")

        ads = self.adsorbate

        mode_counter = 0
        vib_counter = 0
        freq_list = ads.normal_modes.frequency_list
        last_vib_freq = 0

        # two lists to hold the coefficients of the rotations
        # and translations of each normal mode
        t_coef_list = list()
        r_coef_list = list()

        ads.align_to_prinaxes()
        atoms_com_distance = np.linalg.norm(ads.positions, axis=1)

        print(ads.to_molden_format())
        exit()

        # replace small values
        low_values_indices = atoms_com_distance < 0.1
        atoms_com_distance[low_values_indices] = 0.1

        for mode in self.adsorbate.normal_modes:

            mode_counter += 1
            is_pure_vib = None
            freq = freq_list[mode_counter - 1]

            motion_max = np.sqrt(max(mode.motions[:, 0])**2 +
                                 max(mode.motions[:, 1])**2 +
                                 max(mode.motions[:, 2])**2)

            # calculate the motion of the center of mass
            trans_coefs = np.dot(ads.mass_array,  mode.motions) / ads.mass_array.sum()

            # the "trans_tot" variable should give an idea
            # about how much the ceneter-of-mass of the molecule
            # is linearly moving
            trans_tot = np.linalg.norm(trans_coefs)

            # detect whether this is a pure vibration
            if np.abs(motion_max - trans_tot) > 0.3 and trans_tot < 0.005:

                # identify cases where there are two modes
                # with almost exactly the same adsorbate vibration.
                if abs(freq - last_vib_freq) < 10:
                    if np.linalg.norm(mode.motions - ads.normal_modes[mode_counter - 2].motions) > 0.3:
                        is_pure_vib = True
                    else:
                        verbprint(1, self.verbosity, "douplicated vibration found "
                                                     "in normal mode {}".format(mode_counter))

                else:
                    is_pure_vib = True

            if is_pure_vib:
                vib_counter += 1
                last_vib_freq = freq
                verbprint(2, self.verbosity,"{} purevib {: 5f} {: 5f}".format(mode_counter, np.abs(motion_max - trans_tot), trans_tot))

            # if it is not a pure vibration, then is must be a hindered rotation/translation
            if not is_pure_vib:

                # add translational coefficients to list
                t_coef_list.append(trans_coefs.tolist())

                # purify the COM linear motion - so all is left is rotations
                # assuming no significant vibrational motion remains.
                non_trans_motion = mode.motions - trans_coefs

                # project the rotational motion
                # you can look at the following wikipedia page for
                # http://en.wikipedia.org/w/index.php?title=Angular_velocity&oldid=632430097

                angular_velcity_vector = np.cross(non_trans_motion, ads.positions) \
                                         / np.array([(atoms_com_distance ** 2)]).T

                print(mode_counter, angular_velcity_vector)
                # if mode_counter == 13:
                #     exit()
                # TODO: this shouldn't be max, but rather absolute max value
                # TODO: this shouldn't be max, but rather absolute max value
                # TODO: this shouldn't be max, but rather absolute max value
                # TODO: this shouldn't be max, but rather absolute max value
                # TODO: this shouldn't be max, but rather absolute max value
                # TODO: this shouldn't be max, but rather absolute max value
                # TODO: this shouldn't be max, but rather absolute max value
                # TODO: this shouldn't be max, but rather absolute max value

                r_coef_list.append(np.max(non_trans_motion, axis=0).tolist())

        t_array = np.array(t_coef_list)
        t_array /= np.linalg.norm(t_array, axis=0)

        r_array = np.array(r_coef_list)
        r_array /= np.linalg.norm(r_array, axis=0)

        print(r_array**2)

        verbprint(1, self.verbosity, "Found {} pure vibrations,\nexpected 3*N-6 = 3*{}-6 = {}".format(vib_counter, self.adsorbate.atom_count, self.adsorbate.atom_count*3 - 6))

    def set_adsorbate(self, atom_list):
        """
        specify the atom numbers that compose
        the adsorbed molecule

        :param atom_list: a list of atoms
        """

        # this conversion is necessary, because user numbering is
        # different from internal numbering

        atom_list = np.array(atom_list, dtype=int)
        atom_list -= 1
        atom_list = atom_list.tolist()

        self.adsorbate = self.molecule.get_fragment(atom_list)



if __name__ == "__main__":
    hm = HNMA(verbosity=1)

    np.set_printoptions(precision=2, suppress=True)

    hm.parse_molden('../example_outputs/cat-ca-1Me-freq.qchem', molden_job=1)
    adsorbate_atoms = range(14, 19)

    # hm.parse_molden('/home/tsivion/Dropbox/abinitio/metane_storage/catechol-Mg/no_solvent/cat-mg-1Me-opt-freq.qchem', molden_job=2)
    # adsorbate_atoms = range(14, 19)

    # hm.parse_molden('../example_outputs/benzene_opt_freq.qchem', molden_job=2)
    # adsorbate_atoms = range(1, 13)

    hm.set_adsorbate(adsorbate_atoms)

    hm.get_thermodynamics()