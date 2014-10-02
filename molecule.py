from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

__author__ = 'Ehud Tsivion'

import numpy as np
from verbosity import verbprint
import subprocess as sp
import tempfile
import copy
import os
import molden_parser


class Molecule(object):
    """

        A class for holding the information and operations
         related to collections of atoms, such as molecule,
         dimer or any complex.

    """

    """
    some notes to the programmer, about atom numbering:

    like everything in python, the atoms are numbered starting 0
    however, to make things more friendly to the user, he refers
    to the atoms as starting from 1. The numbers entered by the user
    can then be converted to the usual python convention.

    for instance, if there are 15 atoms, and the user wishes to split
    the molecule after the 10'th atom, he would use:

    molecule.split(10) (1 .. numbering)

    and not molecule.split(9) (0 .. numbering)
    """

    def __init__(self, verbosity=1, atoms=None, normal_modes=None):

        self.verbosity = verbosity

        # all the *_list variables are managed by getter and setters
        self._atoms = atoms

        # threshold for a moment of inertia purification
        # and other stuff
        self.mom_thresh = 1E-5

    def add_atom(self, new_atom):
        """

        :param new_atom:
        :return: adds the atom to the atom list of the molecule
        """

        if not self._atoms:
            self._atoms = []

        self._atoms.append(copy.deepcopy(new_atom))

    def rotate(self, atom_num1, atom_num2, angle, atom_list=None, radians=False):
        """

        counter-clockwise rotation of the molecule around the
        axis formed by connecting atom1 and atom2

        :param angle: how much to rotate around the axis formed
        by atoms 1 and 2. Degrees is the defaults

        :param radians: False is angle is given in degrees, True if
        given in radians

        :param atom_list: a list of atoms to be rotated

        :return: Changes the internal state (more specifically: the location)
        of the location of the atoms that make the molecule
        """

        if not atom_list:
            atom_list = np.arange(self.atom_count)

        else:
            # correct user input to start from 0 (not user 1..)
            atom_list = np.array(atom_list, dtype=int)
            atom_list -= 1

        max_num = max(atom_num1, atom_num2)
        if max_num > self.atom_count:
            raise ValueError("atom index {} exceeds the number of atoms in molecule."
                             .format(max_num))

        min_num = min(atom_num1, atom_num2)
        if min_num < 1:
            raise ValueError("atom index {} is smaller than 1".format(min_num))

        verbprint(2, self.verbosity,
                  "rotate molecule around the axis "
                  "which goes through: "
                  "\natom {}: {}\natom {}: {}\nby {} degrees/radians"
                  .format(atom_num1, self.atoms[atom_num1],
                          atom_num2, self.atoms[atom_num1],
                          angle))

        # correct the atom numbering to start from 0
        # since python uses C numbering, which starts at 0 and not 1
        atom_num1 -= 1
        atom_num2 -= 1

        if not radians:
            # convert to radians
            angle = angle * np.pi / 180

        coordinates = self.positions

        # get rotation vector, and normalize
        rot_vec = coordinates[atom_num1] - coordinates[atom_num2]
        rot_vec /= np.linalg.norm(rot_vec)

        # shift the coordinate system, such that (0,0,0) will be located
        # on one of the atoms (atom_num1 or atom_num2) which define
        # the rotation angle
        coord_shift = np.copy(coordinates[atom_num1])
        coordinates -= coord_shift


        # form the rotation matrix
        # see: http://goo.gl/lUCF3P for wikipedia article
        rot_mat = np.zeros((3, 3), dtype=np.float)

        cos_theta = np.cos(angle)
        sin_theta = np.sin(angle)

        rot_mat[0, 0] = cos_theta + rot_vec[0] ** 2 * (1 - cos_theta)
        rot_mat[0, 1] = rot_vec[0] * rot_vec[1] * (1 - cos_theta) - rot_vec[2] * sin_theta
        rot_mat[0, 2] = rot_vec[0] * rot_vec[2] * (1 - cos_theta) + rot_vec[1] * sin_theta

        rot_mat[1, 0] = rot_vec[1] * rot_vec[0] * (1 - cos_theta) + rot_vec[2] * sin_theta
        rot_mat[1, 1] = cos_theta + rot_vec[1] ** 2 * (1 - cos_theta)
        rot_mat[1, 2] = rot_vec[1] * rot_vec[2] * (1 - cos_theta) - rot_vec[0] * sin_theta

        rot_mat[2, 0] = rot_vec[2] * rot_vec[0] * (1 - cos_theta) - rot_vec[1] * sin_theta
        rot_mat[2, 1] = rot_vec[2] * rot_vec[1] * (1 - cos_theta) + rot_vec[0] * sin_theta
        rot_mat[2, 2] = cos_theta + rot_vec[2] ** 2 * (1 - cos_theta)

        coordinates[atom_list] = np.dot(rot_mat, coordinates[atom_list].T).T

        # shift the coordinate system back to original center
        coordinates += coord_shift

        # update the locations of the molecule object
        self.positions = coordinates

    @property
    def mass_array(self):
        """

        Return a numpy float array of the masses of the molecules atoms
        with the same ordering.

        If N is the number of atoms, the size of the array is N x 1

        """
        if not self.atoms:
            raise Exception('There are no atoms set for this molecule')

        mass_list = []

        for atom in self.atoms:
            mass_list.append(atom.mass)

        return np.array(mass_list).astype(float)

    @property
    def positions(self):
        """
        Return a numpy float array of the positions of the molecules atoms
        with the same ordering.

        If N is the number of atoms, the size of the array is N x 3
        """

        position_list = []

        for atom in self.atoms:
            position_list.append(atom._xyz)

        positions = np.array(position_list).astype(float)

        return positions

    @positions.setter
    def positions(self, new_positions):
        if not new_positions.shape[0] == self.atom_count:
            raise IndexError("number new atomic positions do not match number of atoms")

        counter = 0

        for atom in self.atoms:
            atom.coords = new_positions[counter]
            counter += 1

    @property
    def center_of_mass(self):
        """

            Get the center of mass
            return a vector of the form [x_position, y_position, z_position]

            units are Bohr (atomic)

        """

        com = np.dot(self.mass_array, self.positions) / self.mass_array.sum()

        return com

    @property
    def inertia_tensor(self):
        """Get the moments of inertia along the principal axes.

        Units of the moments of inertia are amu*bohr**2.

        :returns a 1 x 3 vector of momenta values, a 3 x 3 matrix of which
        the columns are the principle axes of inertia.

        """
        com = self.center_of_mass
        positions = self.positions
        positions -= com  # translate center of mass to origin
        masses = self.mass_array

        # initialize elements of the inertial tensor
        i_xx = 0.0
        i_yy = 0.0
        i_zz = 0.0
        i_xy = 0.0
        i_xz = 0.0
        i_yz = 0.0

        for i in range(self.atom_count):
            x, y, z = positions[i]
            m = masses[i]

            i_xx += m * (y ** 2 + z ** 2)
            i_yy += m * (x ** 2 + z ** 2)
            i_zz += m * (x ** 2 + y ** 2)
            i_xy += -m * x * y
            i_xz += -m * x * z
            i_yz += -m * y * z

        inertia_tensor = np.array([[i_xx, i_xy, i_xz],
                                   [i_xy, i_yy, i_yz],
                                   [i_xz, i_yz, i_zz]])

        moments, principle_axes = np.linalg.eigh(inertia_tensor)

        sort_perm = moments.argsort()

        # sort them
        moments.sort()
        principle_axes = principle_axes[sort_perm]

        # purify moments if smaller than certain thresh
        moments[moments <= self.mom_thresh] = 0

        return moments, principle_axes

    def split(self, atom_number):
        """
        splits the molecule into several molecules

        :param atom_number: a list number of atoms, which will  be the spliting
        points of the molecule. The molecule is split after the giver atom.

        for instance, if molecule M has 18 atoms
        M.split([5, 6, 12]) would return in the following:
        [M1, M2, M3, M4]

        M1 contains the atoms 0..5
        M2 contains the atom 6
        M3 contains the atoms 7..12
        M4 contains the atoms 13..17

        :return: list of molecules
        """

        if not isinstance(atom_number, list):
            atom_number = [atom_number]

        try:
            point_list = np.array(atom_number, dtype=int)

        except ValueError:
            print('there are non-integer values')

        point_list = point_list[point_list > 0]
        point_list = np.unique(point_list)
        point_list = np.sort(point_list)

        if point_list[-1] >= self.atom_count:
            raise ValueError('Cannot split molecule at atom number {} - '
                             'not enough atoms'.format(point_list[-1]))

        verbprint(2, self.verbosity, "splitting the molecule at "
                                     "atoms number: {}".format(point_list))

        molecule_counter = 0

        molecule_list = list()

        molecule_list.append(Molecule())

        for i in range(self.atom_count):

            molecule_list[molecule_counter].add_atom(self.atoms[i])

            if i in point_list - 1:  # - 1 to correct numbering to start from 1
                # and not from 0 like C lang
                molecule_list.append(Molecule())
                molecule_counter += 1

        return molecule_list

    @property
    def symmetry_number(self):
        return int(1)

    @property
    def atom_count(self):
        return len(self._atoms)

    @property
    def atoms(self):
        return self._atoms

    @atoms.setter
    def atoms(self, new_atoms):
        self._atoms = new_atoms

    @property
    def is_linear(self):
        """
        Returns True value if molecule is linear
         as determined by the Inertia moments
        """

        moments, axes = self.inertia_tensor

        if moments[0] == 0.0 and (moments[1] - moments[2]) < self.mom_thresh:
            linear = True

        else:
            linear = False

        return linear


    def to_jmol(self):
        """
        use JMol to visualize molecule
        """

        # an xyz file that contains all the information about the molecule
        jmol_script = tempfile.NamedTemporaryFile(suffix='.molden', prefix="pyqchemVis_")

        jmol_script.write(self.to_molden_format())

        jmol_script.flush()

        try:
            sp.call(["jmol", "-i", "{}".format(jmol_script.name)])

        except OSError:
            print('cannot find JMOL executable')

        jmol_script.close()


    def to_molden_format(self):
        """
        write the molecular data in xyz format

        :return: the content of an xyz formatted file
        """
        output = ""

        output += "[Molden Format]\n"
        output += "[Atoms] (AU)\n".format(self.atom_count)

        counter = 1

        for atom in self.atoms:
            output += "{:>5} {:>5} {:>5} {:10f} {:10f} {:10f}\n" \
                .format(atom.symbol,
                        counter,
                        atom.atomic_number,
                        atom.coords[0],
                        atom.coords[1],
                        atom.coords[2])
            counter += 1

        return output


    def to_xyz(self):
        """
        Generate an xyz format text with molecular specifications

        :return: xyz format text with molecular specifications
        """

        xyz_output = ""

        xyz_output += "{:>5}\n{}\n".format(self.atom_count, "pyQChem output")

        for atom in self.atoms:
            xyz_output += "{:>5} {:10f} {:10f} {:10f}\n" \
                .format(atom.symbol,
                        atom.coords[0],
                        atom.coords[1],
                        atom.coords[2])

        return xyz_output


    def to_file(self, file_name="pyQChem_output", file_format="xyz"):
        """
        writes a file with molecular specification

        :param file_name: the name of the file to write

        :param file_format: the file format. Currently supported formats are:
                xyz and Molden

        :return: None
        """

        file_format = file_format.lower()

        supported_formats = ["xyz", "molden"]

        if file_format not in supported_formats:
            raise NotImplementedError("format {} is not supported "
                                      "(yet?)".format(file_format))

        if file_format == "xyz":
            file_text = self.to_xyz()
            file_name += ".xyz"

        elif file_format == "molden":
            file_text = self.to_molden_format()
            file_name += ".molden"

        else:
            raise NotImplementedError("format {} is not supported "
                                      "(yet?)".format(file_format))

        f = open(file_name, 'w')
        f.write(file_text)
        f.close()

        verbprint(1, self.verbosity, "write the file {} to "
                                     "directory {}".format(file_name, os.getcwd()))


    def __str__(self):
        output = ""

        for atom in self.atoms:
            output += "{:4} {:.2e} {: e} {: e} {: e}\n" \
                .format(atom.symbol,
                        atom.mass,
                        atom.coords[0],
                        atom.coords[1],
                        atom.coords[2])

        return output


    def __add__(self, other_molecule):
        """
        an operator to add two molecule together

        :return: a molecule object with the atoms of the two molecules
        joined
        """

        joined_molecule = Molecule()

        if self.atoms:
            for atom in self.atoms:
                joined_molecule.add_atom(atom)

        if other_molecule.atoms:
            for atom in other_molecule.atoms:
                joined_molecule.add_atom(atom)

        return joined_molecule


    def __getitem__(self, item):
        raise Exception("__getitem__ not implemented!")


def join(molecule_list):
    new_mol = Molecule()

    for mol in molecule_list:
        new_mol += mol

    return new_mol


def from_molden(file_name, molden_job_num=1, verbosity=1):
    """
    parse a molecule object out of a molden format

    :param: file_name: the file containing the molden format

    :param molden_job_num: the number of molden output. Relevant
                        if a job contains several "molden format"
                        groups, use the output from molden_job.

    :param verbosity: how much output: 0 none, 1 some (default),
                                            2 (or higher) debug.

    :return: a molecule object
    """

    mp = molden_parser.MoldenIO(file_name)

    return mp.molecule


if __name__ == "__main__":
    mp = molden_parser.MoldenIO('./example_outputs/benzene_opt_freq.qchem')
    molecule = mp.molecule
    # molecule.to_jmol()
    molecule.rotate(3, 5, 38, [1, 11, 2, 12])
    molecule.to_jmol()
    # molecules = molecule.split([23])
    # molecules[1].rotate(2, 5, 90)
    # join(molecules).to_jmol()


