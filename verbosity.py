from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import

__author__ = 'Ehud Tsivion'


def verbprint(intrinsic_verbosity, requested_verbosity, print_this):
    """"
    a customized print function:

    print only if the requested verbosity degree is larger or
    equal from the intrinsic verbosity
    """

    if requested_verbosity >= intrinsic_verbosity:
        # return the standard print function (Python3 style)
        return print(print_this)

    else:
        # do noting
        return lambda *a, **k: None