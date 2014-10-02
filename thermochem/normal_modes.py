from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

__author__ = 'Ehud Tsivion'

import numpy as np


class NormalModes(object):
    """
    A class to handle the collection of motions

    if N_modes is the number of normal modes
    and N_atoms is the number of atoms

    then the _modes_list is an array with dimension
    _modes_list N_nodes x N_atoms x 3

    for example:

    _modes_list[4, 5, 0] is the x_axis value of the 5th atom in the
    4th normal node.

    """

    def __init__(self):
        self._modes_list = []

    def add_mode(self, mode):
        self._modes_list.append(mode)

    @property
    def count(self):
        return len(self._modes_list)

    @property
    def modes(self):
        return self._modes_list

    @property
    def frequency_list(self):
        """

        :return: a list of the frequencies of the normal modes
        """

        freq_list = []

        for mode in self._modes_list:
            freq_list.append(mode.frequency)

        return freq_list

    @property
    def frequency_array(self):
        """

        :return: a numpy float array of the frequencies of the normal modes
        """
        return np.array(self.frequency_list, dtype=np.float)

    def __str__(self):
        counter = 0
        text = ""
        for mode in self._modes_list:
            counter += 1
            text += "NMode {}\n".format(counter)
            text += "{}\n".format(mode)
        return text

