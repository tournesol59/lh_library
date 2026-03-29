/*********TO DO
* Suivre tuto ctypes de realpython
* Ecrire ici un eclasse en C++ que l'on appellera depuis python comme librairie
****************************************/
#ifndef __CPP_CTYPES__
#define __CPP_CTYPES__

#include <fstream>
#include <cstdio>
#include <iostream>
#include <math.h>
#include <cstring>

class  Cppfeedingmat {
  
  public:
   Cppfeedingmat(const char * FileName, int iNorder, double * endvalues);
   ~Cppfeedingmat();
   setKernelMat();

  private:
   char * outFileName;
   int Norder;
   double * c_points;
   double * c_kernelmat;
}

#endif
