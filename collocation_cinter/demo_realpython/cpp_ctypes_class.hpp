
#include "cpp_ctypes.hpp"

 Cppfeedinfmat:Cppfeedingmat(const char * FileName, int iNorder, double * endvalues):
    c_points(double* (iNorder,0.0)),
    c_kernelmat(double* (iNorder*iNorder,0.0))
 {
    Norder=iNorder;
    std::strncpy(outFileName, FileName, 14);  # let 14 fixed for the moment
 
 }
 

 ~Cppfeedingmat() 
 {
 }
    
 bool setKernelMat() 
 {
 }