"""This modules provides the interface to BSMPT to calculate
the nucleation temperature of a cosmological first order phase transition."""
import ctypes
import pathlib

import numpy as np

# Import the library
libname = pathlib.Path().absolute() / "libBSMPTinterface.so"
c_lib = ctypes.CDLL(libname)
c_lib.getNuclTemp.restype = ctypes.c_double

def getNuclTemp(D, A, lam, T0, verbose=False):
    """Calculate the nucleation temperature of the phase transition.

    Parameters
    ----------


    Returns
    ----------
    Tnucl : float
        The nucleation temperature of the phase transition."""

    # Prepare the input for C++ function call
    _D = ctypes.c_double(D)
    _A = ctypes.c_double(A)
    _l = ctypes.c_double(lam)
    _T0 = ctypes.c_double(T0)

    # Call the function getNuclTemp from libBSMPTinterface.so
    Tnucl = c_lib.getNuclTemp(_D, _A, _l, _T0)

    # Check for sensible output
    if np.isnan(Tnucl) or Tnucl == 0.0:
        if verbose:
            print("Could not calculate the nucleation temperature!")
            print("This could be because of a 2nd order transition" + \
                  "or no phase transition at all.")
            print("Returning 0.")
        return 0.0
    return Tnucl


# Example input parameters for which we get a
# first order phase transition
if __name__=="__main__":
    D = 1
    A = 1.5
    lam = 2
    T0 = 10

    Tnucl = getNuclTemp(D, A, lam, T0, verbose=True)
    print("Found a first order phase transition with a nucleation temperature of Tnucl = {:2.3f}".format(Tnucl))
