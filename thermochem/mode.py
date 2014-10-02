from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

__author__ = 'Ehud Tsivion'

import numpy as np


class NMode(object):
    """
    a molecular motion is collaborative movement of its atoms.


    """

    def __init__(self, coords_list, frequency):
        self.coordinates = np.array(coords_list).astype(np.float)
        self.frequency = frequency

    @property
    def velocity(self):
        """
        The velocity of a cluster of atoms (usually a molecule)
        is a collective movement in a certain direction, assuming
        that the molecule behaves as a rigid body, where
        the distance between all the atoms remains constant.

        This procedure extract the collective motion off all the modes
        by finding the minimal motion of each of the normal modes along
        each of the cartesian axes.

        :return: motion_array. A 1 x 3 vector in the cartesian representation
        which indicates the velocity (speed and direction) common to all atoms
        for this normal mode:

            motion_vector = [speed along x axis, speed along y axis, speed along z axis]
        """

        motion_list = self.coordinates

        # the velocity vector of this normal mode
        motion_vector = np.array([abs_min(motion_list[:, 0]),
                                  abs_min(motion_list[:, 1]),
                                  abs_min(motion_list[:, 2])])

        return motion_vector

    @property
    def rotation_projection(self, principle_axes):
        """
        The projection of the rotation of a cluster of atoms (usually a molecule)
        is a collective rotation_projection along the principle axes of inertia,
        assuming that the molecule behaves as a rigid body, where
        the distance between all the atoms remains constant.

        This procedure extract the collective rotation off all the modes
        by finding the minimal motion of each of the normal modes along
        each of principle axes of the moment of inertia.

        :param principle_axes: a matrix (3 x 3) of which the columns are
        the normalized principle axes of inertia.

        :return: rotations_vector: an array of 1 x 3, which indicates the rotation
        along the principle axes which common to all atoms

        for this normal mode:
            rotations_vector = [rotation prin. ax. 1, rotation prin. ax. 2, rotation prin. ax. 3]
        """

        motion_list = self.coordinates

        projected_rotation = np.dot(motion_list, principle_axes)

        rotations_vector = np.array([abs_min(projected_rotation[:, 0]),
                                     abs_min(projected_rotation[:, 1]),
                                     abs_min(projected_rotation[:, 2])])

        return rotations_vector

    @property
    def size(self):
        return len(self.coordinates)


    def __str__(self):
        return np.array_str(self.coordinates, precision=8)


def abs_min(some_array):
    some_array_all_positive = (some_array >= 0.0).all()
    some_array_all_negative = (some_array < 0.0).all()

    # if all term in array have positive sign
    # return the minimum (smallest positive)
    if some_array_all_positive:
        return some_array.min()

    # if all term in array have negative sign
    # return the maximum ((absolute) smallest negative)
    elif some_array_all_negative:
        return some_array.max()

    # if some_array has mixed sings return 0.0
    else:
        return np.float(0.0)


if __name__ == "__main__":
    pass