#!/usr/bin/env python
""" Simple examples of calling C functions through ctypes module. """
import ctypes
import sys
import pathlib

if __name__ == "__main__":
    libname = pathlib.Path("C:\MinGW\msys\1.0\home\Alain\lh_library\collocation_cinter\demo_realpython\").absolute()		
print("libname: ", libname)

#***stackoverflow
#BASE_DIR = Path('.')
#EXTENSIONS = {'.xls', '.txt'}

#for path in BASE_DIR.glob(r'**/*'):
#    if path.suffix in EXTENSIONS:
#        print(path)
#***

    # Load the shared library into c types.
#    if sys.platform.startswith("win"):
#        c_lib = ctypes.CDLL(libname / "cmult.dll")
#    else:
#        c_lib = ctypes.CDLL(libname / "libcmult.so")
    c_lib = ctypes.CDLL(libname / "libcmult.so")
    # Sample data for our call:
    x, y = 6, 2.3

    # This will not work:
    # answer = c_lib.cmult(x, y)

    # This produces a bad answer:
    answer = c_lib.cmult(x, ctypes.c_float(y))
    print(f"    In Python: int: {x} float {y:.1f} return val {answer:.1f}")
    print()

    # You need tell ctypes that the function returns a float
    c_lib.cmult.restype = ctypes.c_float
    answer = c_lib.cmult(x, ctypes.c_float(y))
    print(f"    In Python: int: {x} float {y:.1f} return val {answer:.1f}")
