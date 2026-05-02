#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "../include/ident05_coll.hpp"
#include "../include/ident05_fftdata.hpp"
#include <memory>

using namespace lhlib;

static PyObject* ident05_system(PyObject *self, PyObject *args)
{
    const char *argFileName;
    char strFileName[14];
    const char *argCodeName;
    char strCodeName[14];
    const char *argSpecName;
    char strSpecName[14];

    int lenFileName=14;
    int num_totalpoints;

    if (!PyArg_ParseTuple(args, "s", &argFileName))
	    return NULL;
    if (!PyArg_ParseTuple(args, "s", &argCodeName))
	    return NULL;
    if (!PyArg_ParseTuple(args, "s", &argSpecName))
	    return NULL;
    if (!PyArg_ParseTuple(args, "i", num_totalpoints))
	    return NULL;
    // comment: il faut plutot grouper les expressions suivantes

}

