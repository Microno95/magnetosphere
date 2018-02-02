# magnetosphere
A set of python modules for analysing the output of the Gorgon MHD code.

It contains the following functionality:

* Importing Gorgon `pvti` files into Numpy arrays, in `gorgon/gorgon_import.py`
* Saving scalar 3D Numpy arrays to `vti`, in `gorgon/gorgon_import.py`
* A streamtracer, written in Fortran and parallelised with OpenMP, in `gorgon/streamline.py`
* A magnetic field tracer with connectivity mapping, in `gorgon/fieldlines.py` and `gorgon/connectivity.py`
* Null point determination (untested), in `gorgon/null_points.py`
* Empirical Earth Bow shock and Magnetopause models, in `models/bowshock.py` and `models/magnetopause.py`

## Fortran subroutines and F2Py

The Fortran source files are contained in `gorgon/fortran/`. This repo contains compiled versions for Windows and Linux, on numpy v1.13.3.
If you need it to run on other platforms, they need to be compiled using F2Py by running the shell script `gorgon/fortran/compile.sh`.
There have been issues using the pre-compiled versions with other versions of numpy, so either update/downgrade numpy, or recompile.