/*******************************************************
 * Class IDENT05_IODATA defined here
 * the methods are previewed to import format-specific input 
 * files data for collocation inputs (e.g. fft for verification),
 * and also to EXPORT a list-previewed output data in a file
 *
 * Note that the list container is a list of double, which is
 * filled alternately by a x value then a y value ..
 *
 * LIEBHERR TOULOUSE
 *******************************************************/
#include "../include/lhTypes.hpp"
//#include "../MathFunctions/MathFunctions.hpp"
#include <cassert>
#include <iostream>
#include <istream>
#include <sstream>
#include <fstream>
#include <malloc.h>
#include <cstdio>
#include <vector>  // for read_parse_file
#include <list>
#include <array>
#include <cstring>
#include <math.h>

#ifndef __IDENT05_FFTDATA_HPP__
#define __IDENT05_FFTDATA_HPP__

using namespace lhlib;

//using namespace std;

typedef std::list<dpair> li_doubles;
//typedef list<Number[6]> li_6doubles;

class IDENT05_IODATA {

	public:
// constructor
  IDENT05_IODATA(Number size, char* giveFileName) ;

  ~IDENT05_IODATA();

  Index sizefft;
  Number samplefft;

  Index  sizeevalsol; // length of following table
    
   char FileName[14]; // name of the file for export (exportToDisk)

    	// vector imported from y text file
	// previewed for FFT externally performed whose output file format
	// is in a one-dimensional vector [x,y]
   bool read_extern_output(Index dim, li_doubles &li_reals);

	// for solutions of collocation or ident algorithms to export
	// on "FileName" (arg of instantiation) file (e.g. ".bat")
   bool exportToDisk(li_doubles &li_reals);

	private:
};


#endif
