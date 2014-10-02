from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

__author__ = 'Ehud Tsivion'

from verbosity import verbprint
import constants as const
import numpy as np


class StdThermo(object):
    """

    """

    def __init__(self, molecule, normal_modes, verbosity=1, linear=False, temperature=297.15):

        self._molecule = molecule
        self._normal_modes = normal_modes
        self._verbosity = verbosity
        self._linear = linear
        self._temperature = temperature
        self._internal_energy = None
        self._entropy = None

        verbprint(1, verbosity, '\n{:*^40}'.format(' Std Thermodynamic properties '))

        # validate input correctness
        # Is this a linear molecule, what is the expected number of modes?

        self._linear = molecule.is_linear

        if self._linear:
            verbprint(1, verbosity, "This is a linear molecule")
            expected_modes = molecule.atom_count * 3 - 5

        else:
            verbprint(1, verbosity, "This os a non-linear molecule")
            expected_modes = molecule.atom_count * 3 - 6

        if expected_modes == normal_modes.count and not self._linear:
            verbprint(1, verbosity, "Validation of 3N-6 normal modes: ... OK")

        elif expected_modes == normal_modes.count and self._linear:
            verbprint(1, verbosity, "Validation of 3N-5 normal modes: ... OK")

        else:
            if not self._linear:
                raise Exception('There should be 3N-6 normal modes (N={}). '
                                'However, there are only {}.\n'
                                .format(molecule.atom_count, normal_modes.count))
            else:
                raise Exception('There should be 3N-5 normal modes (N={}). '
                                'However, there are only {}.\n'
                                .format(molecule.atom_count, normal_modes.count))

    @property
    def temperature(self):
        """
        The temperature which is set of the calculation
        """
        return self._temperature

    @temperature.setter
    def temperature(self, temperature):
        self._temperature = temperature

    @property
    def s_rotations(self):
        return None

    @property
    def s_translations(self):
        """

        :return:
        """

        p_cons = const.Planck_constant
        blz_cons = const.Boltzmann_constant
        tmprt = self._temperature

        # obtain principle moments and axes
        prm_moments, prm_axes = self._molecule.inertia_tensor

        # the rotational temperature for each principle axis (a, b, c)
        theta_a = p_cons**2 / 8 * np.pi * blz_cons * prm_moments[0]
        theta_b = p_cons**2 / 8 * np.pi * blz_cons * prm_moments[1]
        theta_c = p_cons**2 / 8 * np.pi * blz_cons * prm_moments[2]

        q_rot = np.sqrt(np.pi)/sigma * np.sqrt(tmprt**3/(theta_a*theta_c*theta_b))

        return None

    @property
    def s_vibrations(self):
        """

        :return: the internal energy of vibrations in [kJ/mol] units
        """

        # obtain an array of the frequencies of the normal modes
        freq_list = self._normal_modes.frequency_array

        summator = 0

        T = self._temperature

        for freq in freq_list:

            # theta is the characteristic vibrational temperature
            # units of theta are [J]/([J]/[K]) = [K]
            theta = const.inverse_cm_to_joule * freq / const.Boltzmann_constant


            # summator doesn't have units
            summator += (theta / T) * (1/(np.exp(theta/T) - 1)) - np.log(1-np.exp(-1*theta/T))

        # units are [J]/[Mol][K]
        rotational_entropy = const.molar_gas_constant * summator

        # convert to [kj]/[Mol][K]
        rotational_entropy = rotational_entropy * 1E-3

        return rotational_entropy
    @property
    def zpe(self):
        return None

    @property
    def h_rotations(self):
        """

        :return: the internal energy of rotations, which is 1/2RT for
        each degree of freedom. Units are kJ/Mol
        """

        temperature = self._temperature
        gas_const = const.molar_gas_constant

        # units are [J]/[Mol]
        h_rot = 1.5 * gas_const * temperature

        # convert to [kJ]/[mol]
        h_rot *= 1E-3

        return h_rot

    @property
    def h_translations(self):
        return None

    @property
    def h_vibrations(self):
        """

        \langle E_{vib}\rangle=R\sum_{i=1}^{3N-6}\theta_{v_i} \left[\frac{1}{2}+\frac{e^{-\theta_{v_i}/T}}{1-e^{-\theta_{v_i}/T}}\right


        :return: the internal energy of vibrations in [kJ/mol] units
        """

        # obtain an array of the frequencies of the normal modes
        freq_list = self._normal_modes.frequency_array

        summator = 0

        for freq in freq_list:

            # theta is the characteristic vibrational temperature
            # units of theta are [J]/([J]/[K]) = [K]
            theta = const.inverse_cm_to_joule * freq / const.Boltzmann_constant

            # a temporary term to make the whole thing more readable
            temp = np.exp(-1 * theta / self._temperature)

            # units of summator are [K] - same as theta
            summator += theta * (0.5 + temp / (1 - temp))

        # units are [K] * [J]/[Mol][K] = [J]/[Mol]
        internal_energy = const.molar_gas_constant * summator

        # convert to [kj]/[Mol]
        internal_energy *= 1E-3

        return internal_energy

    @property
    def enthalpy(self):
        return None

    @property
    def entropy(self):
        return None

    @property
    def gibbs(self):
        return None

    @property
    def helmholtz(self):
        return None


    def __str__(self):
        return self.gibbs

if __name__ == "__main__":

    from molden_parser import MoldenIO
    mol = MoldenIO("../example_outputs/acetylene_opt_freq.qchem", 2)
    calc = StdThermo(molecule=mol._atoms,
                     normal_modes=mol._normal_modes)
    print("Internal energy of vibrations: {}".format(calc.s_vibrations))


