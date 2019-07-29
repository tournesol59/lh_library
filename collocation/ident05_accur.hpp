/*
**************************************************
 Module part of lh_library collocation "accuracy
 test" in order to test the solution, whichever 
 collocation method (tau, inverse) was used
 
 LIEBHERR TOULOUSE
**************************************************
*/

#ifndef __IDENT05_ACCUR_HPP__
#define __IDENT05_ACCUR_HPP__

#include <cassert>
#include <iostream>
#include <istream>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <array>
#include <cstring>

using namespace lhlib;

class IDENT05_ACCUR {
	public:
		// constructor
	bool IDENT05_ACCUR(const char* TestNameOut, const char* SpecNameOut);
	// TestNameOut: out file of ident05_coll.SolveNumRangesSys_ref1(), contains time series data
	// SpecNameOut: out file(2) with spectral values, completed after by re-writing from Matlab (fft)
	// destructor:
	bool ~IDENT05_ACCUR();

	bool LSQ_SpecError();

	bool LSQ_SeriesError();
}

#endif

