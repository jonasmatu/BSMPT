"""This modules provides the interface to BSMPT to calculate
the nucleation temperature of a cosmological first order phase transition."""
import ctypes
import pathlib

import numpy as np

# Import the library
libname = pathlib.Path().absolute() / "libBSMPTinterface.so"
c_lib = ctypes.CDLL(libname)
c_lib.getNuclTemp.restype = ctypes.c_double

def getNuclTemp(A, lambda3, lambda4, T0, verbose=False):
    """Calculate the nucleation temperature of the phase transition.

    Parameters
    ----------


    Returns
    ----------
    Tnucl : float
        The nucleation temperature of the phase transition."""

    # Prepare the input for C++ function call
    _A = ctypes.c_double(A)
    _l3 = ctypes.c_double(lambda3)
    _l4 = ctypes.c_double(lambda4)
    _T0 = ctypes.c_double(T0)

    # Call the function getNuclTemp from libBSMPTinterface.so
    Tnucl = c_lib.getNuclTemp(_A, _l3, _l4, _T0)

    # Check for sensible output
    if np.isnan(Tnucl) or Tnucl == 0.0:
        if verbose:
            print("Could not calculate the nucleation temperature!")
            print("This could be because of a 2nd order transition \
            or no phase transition at all.")
            print("Returning 0.")
        return 0.0
    return Tnucl


# Example input parameters for which we get a
# first order phase transition
if __name__=="__main__":
    A = 1
    lambda3 = 0.5
    lambda4 = 2
    T0 = 10

    Tnucl = getNuclTemp(A, lambda3, lambda4, T0)
    print("Found a first order phase transition with a nucleation temperature of Tnucl = {:2.3f}".format(Tnucl))
