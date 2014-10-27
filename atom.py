# This is a free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ASE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with ASE.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

__author__ = 'Ehud Tsivion'

import numpy as np
import constants


class Atom(object):

    def __init__(self, sym=None, coords=[], coord_units="bohr"):
        """

        default constructor of atom object

        :param sym: chemical symbol of atom
        :param coords: cartesian coordinates of the atom, in units
        :param coord_units: the units in which coordinates are provided, can
                be bohr or angstroms, angs, Angs etc.. anthing with first letter 'A'
        """

        # make sure chemical symbol has the correct
        # letter-case form
        sym = correct_symbol_case(sym)
        self.symbol = str(sym)

        # set mass
        self.mass = mass_for_sym(sym)

        # verify correctness of input
        if len(coords) == 3:

            # if units are given in Angstroms, convert to bohr
            if coord_units[0].lower() == "a":
                self._xyz = np.array(coords).astype(np.float64)*constants.angstrom_to_bohr

            else:
                self._xyz = np.array(coords).astype(np.float64)

        else:
            raise ValueError("error reading the coordinates: \"{}\"".format(coords))

    @property
    def coords(self):
        return self._xyz

    @coords.setter
    def coords(self, new_coord):
        self._xyz = new_coord

    @property
    def atomic_number(self):
        return number_for_sym(self.symbol)

    def __str__(self):
        return "{} {} {:< 10.8f} {:< 10.8f} {:< 10.8f}".format(self.symbol, self.mass, self._xyz[0], self._xyz[1],
                                                               self._xyz[2])


chemical_symbols = ['X', 'H', 'He', 'Li', 'Be',
                    'B', 'C', 'N', 'O', 'F',
                    'Ne', 'Na', 'Mg', 'Al', 'Si',
                    'P', 'S', 'Cl', 'Ar', 'K',
                    'Ca', 'Sc', 'Ti', 'V', 'Cr',
                    'Mn', 'Fe', 'Co', 'Ni', 'Cu',
                    'Zn', 'Ga', 'Ge', 'As', 'Se',
                    'Br', 'Kr', 'Rb', 'Sr', 'Y',
                    'Zr', 'Nb', 'Mo', 'Tc', 'Ru',
                    'Rh', 'Pd', 'Ag', 'Cd', 'In',
                    'Sn', 'Sb', 'Te', 'I', 'Xe',
                    'Cs', 'Ba', 'La', 'Ce', 'Pr',
                    'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
                    'Tb', 'Dy', 'Ho', 'Er', 'Tm',
                    'Yb', 'Lu', 'Hf', 'Ta', 'W',
                    'Re', 'Os', 'Ir', 'Pt', 'Au',
                    'Hg', 'Tl', 'Pb', 'Bi', 'Po',
                    'At', 'Rn', 'Fr', 'Ra', 'Ac',
                    'Th', 'Pa', 'U', 'Np', 'Pu',
                    'Am', 'Cm', 'Bk', 'Cf', 'Es',
                    'Fm', 'Md', 'No', 'Lr']

# this should have the same size as "chemical_symbols" list
atomic_names = [
    '', 'Hydrogen', 'Helium', 'Lithium', 'Beryllium', 'Boron',
    'Carbon', 'Nitrogen', 'Oxygen', 'Fluorine', 'Neon', 'Sodium',
    'Magnesium', 'Aluminium', 'Silicon', 'Phosphorus', 'Sulfur',
    'Chlorine', 'Argon', 'Potassium', 'Calcium', 'Scandium',
    'Titanium', 'Vanadium', 'Chromium', 'Manganese', 'Iron',
    'Cobalt', 'Nickel', 'Copper', 'Zinc', 'Gallium', 'Germanium',
    'Arsenic', 'Selenium', 'Bromine', 'Krypton', 'Rubidium',
    'Strontium', 'Yttrium', 'Zirconium', 'Niobium', 'Molybdenum',
    'Technetium', 'Ruthenium', 'Rhodium', 'Palladium', 'Silver',
    'Cadmium', 'Indium', 'Tin', 'Antimony', 'Tellurium',
    'Iodine', 'Xenon', 'Caesium', 'Barium', 'Lanthanum',
    'Cerium', 'Praseodymium', 'Neodymium', 'Promethium',
    'Samarium', 'Europium', 'Gadolinium', 'Terbium',
    'Dysprosium', 'Holmium', 'Erbium', 'Thulium', 'Ytterbium',
    'Lutetium', 'Hafnium', 'Tantalum', 'Tungsten', 'Rhenium',
    'Osmium', 'Iridium', 'Platinum', 'Gold', 'Mercury',
    'Thallium', 'Lead', 'Bismuth', 'Polonium', 'Astatine',
    'Radon', 'Francium', 'Radium', 'Actinium', 'Thorium',
    'Protactinium', 'Uranium', 'Neptunium', 'Plutonium',
    'Americium', 'Curium', 'Berkelium', 'Californium',
    'Einsteinium', 'Fermium', 'Mendelevium', 'Nobelium',
    'Lawrencium', 'Unnilquadium', 'Unnilpentium', 'Unnilhexium']

# this should have the same size as "chemical_symbols" list
# the atomic masses are given in atomic mass units
atomic_masses = np.array([
    0.00000,  # X
    1.00794,  # H
    4.00260,  # He
    6.94100,  # Li
    9.01218,  # Be
    10.81100,  # B
    12.01100,  # C
    14.00670,  # N
    15.99940,  # O
    18.99840,  # F
    20.17970,  # Ne
    22.98977,  # Na
    24.30500,  # Mg
    26.98154,  # Al
    28.08550,  # Si
    30.97376,  # P
    32.06600,  # S
    35.45270,  # Cl
    39.94800,  # Ar
    39.09830,  # K
    40.07800,  # Ca
    44.95590,  # Sc
    47.88000,  # Ti
    50.94150,  # V
    51.99600,  # Cr
    54.93800,  # Mn
    55.84700,  # Fe
    58.93320,  # Co
    58.69340,  # Ni
    63.54600,  # Cu
    65.39000,  # Zn
    69.72300,  # Ga
    72.61000,  # Ge
    74.92160,  # As
    78.96000,  # Se
    79.90400,  # Br
    83.80000,  # Kr
    85.46780,  # Rb
    87.62000,  # Sr
    88.90590,  # Y
    91.22400,  # Zr
    92.90640,  # Nb
    95.94000,  # Mo
    np.nan,  # Tc
    101.07000,  # Ru
    102.90550,  # Rh
    106.42000,  # Pd
    107.86800,  # Ag
    112.41000,  # Cd
    114.82000,  # In
    118.71000,  # Sn
    121.75700,  # Sb
    127.60000,  # Te
    126.90450,  # I
    131.29000,  # Xe
    132.90540,  # Cs
    137.33000,  # Ba
    138.90550,  # La
    140.12000,  # Ce
    140.90770,  # Pr
    144.24000,  # Nd
    np.nan,  # Pm
    150.36000,  # Sm
    151.96500,  # Eu
    157.25000,  # Gd
    158.92530,  # Tb
    162.50000,  # Dy
    164.93030,  # Ho
    167.26000,  # Er
    168.93420,  # Tm
    173.04000,  # Yb
    174.96700,  # Lu
    178.49000,  # Hf
    180.94790,  # Ta
    183.85000,  # W
    186.20700,  # Re
    190.20000,  # Os
    192.22000,  # Ir
    195.08000,  # Pt
    196.96650,  # Au
    200.59000,  # Hg
    204.38300,  # Tl
    207.20000,  # Pb
    208.98040,  # Bi
    np.nan,  # Po
    np.nan,  # At
    np.nan,  # Rn
    np.nan,  # Fr
    226.02540,  # Ra
    np.nan,  # Ac
    232.03810,  # Th
    231.03590,  # Pa
    238.02900,  # U
    237.04820,  # Np
    np.nan,  # Pu
    np.nan,  # Am
    np.nan,  # Cm
    np.nan,  # Bk
    np.nan,  # Cf
    np.nan,  # Es
    np.nan,  # Fm
    np.nan,  # Md
    np.nan,  # No
    np.nan]).astype(np.float)  # Lw

def correct_symbol_case(chem_symbol):
    """
    correct the chemical symbol to be in the form:

    Upperlower
    example:

    correct_symbol_case('he') -> 'He'
    correct_symbol_case('HE') -> 'He'
    correct_symbol_case('He') -> 'He'

    correct_symbol_case('C') -> 'C'
    correct_symbol_case('c') -> 'c'


    :param chem_symbol: the chemical symbole to correct
    :return: a correct case chemical symbol
    """

    # make sure chemical symbol is not longer than 2.
    sym_length = len(chem_symbol)

    if sym_length > 2:
        raise Exception('Found error while reading the chemical symbol \"{}\"'.format(chem_symbol))

    elif sym_length == 2:
        return chem_symbol[0].upper() + chem_symbol[1].lower()

    elif sym_length == 1:
        return chem_symbol.upper()

def mass_for_sym(symbol):
    """

    :param symbol:
    :return: the atomic mass of the chemical symbol
    """

    if len(chemical_symbols) != len(atomic_masses):
        raise Exception('\"chem_symbol\" list and \"atomic_masses\" numpy array are not '
                        'the same size.')

    try:
        index = chemical_symbols.index(symbol)

    except ValueError():
        print('Chemical symbol \"{}\" does not exist'.format(symbol))
        return 0

    mass = constants.dict_of_atomic_masses[symbol]

    return mass

def number_for_sym(symbol):
    """

    :param symbol:
    :return: the atomic number of the chemical symbol
    """

    try:
        number = chemical_symbols.index(symbol)

    except ValueError:
        print("{} is not a symbol of an atom".format(symbol))

    return number









