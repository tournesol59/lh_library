/********************************************************************************************
* The aim is to test the collocation method from python 
*    y'' + 1/2*y' + 16*y = 0  (1) with y(0)=0 y'(0)=1
*  TODO: first of all a basic wrapper (copy of parts from collocation/main_example.cpp) with only file parsing test
*      further difficulty: link as a library
*********************************************************************************************/
#include "Python.h"

#include "../include/ident05_coll.hpp"
#include "../include/ident05_fftdata.hpp"
#include <memory>  // for Smart Ptr: unique_ptr.h

//#include "lapacke.h"  // prefer set INCLUDE DIR in setup.py for portability

static PyObject *ExtLapackError = NULL;
static PyObject * idcolloc_system(PyObject *self, PyObject *args);


using namespace lhlib;

static PyMethodDef IdCollocMethods[] = {
   {"system", idcolloc_system, METH_VARARGS,
	    "numeric matrix multiplication"},
   {NULL,NULL,0,NULL}
};

static struct PyModuleDef idcollocmodule = {
   PyModuleDef_HEAD_INIT,
   "idcolloc", // name of the module
    NULL, // module documentation, may be NULL (matmult_doc)
   -1,        // size of per-interpreter state of the module
              // or -1 if the module keeps state in global variables
   IdCollocMethods  // defined above
};


PyMODINIT_FUNC PyInit_idcolloc(void) {
   PyObject *m;

   m = PyModule_Create(&idcollocmodule);
   if (m==NULL) 
      return NULL;
   ExtLapackError = PyErr_NewException("idcolloc.lapack.error", NULL, NULL);
   Py_XINCREF(ExtLapackError);
   if (PyModule_AddObject(m, "error", ExtLapackError) < 0) {
       Py_XDECREF(ExtLapackError);
       Py_CLEAR(ExtLapackError);
       Py_DECREF(m);
       return NULL;
   }
   return m;
}

static PyObject * idcolloc_system(PyObject *self, PyObject *args) {
 
//int main(int argc, char ** argv) {
  const char * strFileName;
  const char * strCodeName;
  int * num_totalpoints;
  int lenFileName=14;
  int i, sts;

  if (!PyArg_ParseTuple(args, "(ssd)", &strFileName, &strCodeName, &num_totalpoints))
      return NULL;
 
  std::cout << "Have string " << strFileName << " of length " << lenFileName << " as arg[1] and " << strCodeName << " as arg[2]\n";

  // list of import export vectors
  li_doubles datafft;  // li_doubles is a list<double> container
  li_doubles datatimesol;  // idem

  // solver class
  std::cout << "Instanciate the collocation class \n"; 
  IDENT05_COLL cnClassInst=IDENT05_COLL(6,strFileName,strCodeName);

  // data solution container and class list<T>
  char fftFileName[14];
  strncpy(fftFileName,"fft.out", 8); //..  
  IDENT05_IODATA fftClassInst=IDENT05_IODATA(256, 0.1, fftFileName); // size, sampling, file to read from colloc

  char expFileName[14];
  strncpy(expFileName,"out.bat", 8); //..
  IDENT05_IODATA expClassInst=IDENT05_IODATA(num_totalpoints, 0.1, expFileName); // two signals
  
  // read problem formulation:
  std::cout << "Start the reading of data file " << strFileName << "\n";
  cnClassInst.read_parse_file();
  std::cout << "Start the reading of codes file " << strCodeName << "\n";  
  cnClassInst.read_parse_code();

  std::cout << "complete first test for the moment \n";
  //std::cout << "Start the collocation on each interval \n";
  /*
   * AFTER
// If Testing: Decompose:
//  cnClassInst.ExpandSeriesLinearSys_ref1();
//  cnClassInst.SolveSeriesLinearSys_ref1();
// or complete auto solving:
  cnClassInst.SolveNumRangesSys_ref1();

    // data exchange between classes
  cnClassInst.pass_dataarray_col( 40, datatimesol);
  //Index select[2]; //={1,2}
  expClassInst.exportToDisk(datatimesol);
// import fft after being verified
//  fftClassInst.read_extern_output(256, datafft); // TBC

  *
  */
  return 0;
}

