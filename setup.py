from setuptools import setup, Extension

setup(
        name="idcolloc",
        version="0.0.1",
        description="idcolloc module",
        ext_modules=[Extension("idcolloc", sources=["idcollocmodule.c"],
           include_dirs=["opt/OpenBLAS/include/lapacke.h"], 
           library_dirs=["/usr/lib/x86_64-linux-gnu","/opt/OpenBLAS/lib"],
           libraries=["lapack", "openblas"])]
)
