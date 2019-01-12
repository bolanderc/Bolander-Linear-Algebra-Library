#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Calculates machine precision.

Routine Listings
-----------------
smaceps()    Calculates and returns machine epsilon in single precision.
dmaceps()    Calculates and returns machine epsilon in double precision.

Notes
------
A description of the algorithm is given in Uri and Greif on page 23 and
on Wikipedia (https://en.wikipedia.org/wiki/Machine_epsilon#cite_note-1).

@author: Christian Bolander
Math 5610
HW 1 Task 1

"""


def smaceps():
    """Calculates and returns machine epsilon in single precision.

    Returns
    -------
    mach_eps : float
        Single precision machine epsilon.

    """
    mach_eps = 2**-(24 - 1)/2
    return mach_eps


def dmaceps():
    """Calculates and returns machine epsilon in double precision.

    Returns
    -------
    mach_eps : float
        Double precision machine epsilon.

    """
    mach_eps = 2**-(53 - 1)/2
    return mach_eps
