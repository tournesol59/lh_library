#include <iostream>
//#include <generator>
#include <random>
#include <vector>
#include <iterator>
#include <numeric>
#include <algorithm>
#include <cstring>

#include "../include/ident05_fftdata.hpp"

using namespace std;

int generate_rdom_example( vector<double> &sig, int n, double ca, double cb, double fvar) {

  int i;

  double x;

//  double mean=0.0;
//  double stddev = 0.5;
//  std::default_random_engine generator;
//  std::normal_distribution<double> dist(mean, stddev);

  std::random_device rd;
  std::mt19937 gen( rd());
  uniform_int_distribution<int> dist(0,99);
  for (i=0; i<n; i++) {
    x=ca+cb*1*i + (dist(gen)-50)*fvar/100.0;
    sig.push_back(x);
  }
  for (auto it=sig.begin(); it != sig.end(); it++) {
     std::cout << " " << (*it);
  }
  std::cout<< "\n";

  return 1;
}

int generate_arN_example(std::vector<double> &sig, int n,
		     int order, double x0, std::vector<double> &ar, double fvar)
{
  int i,j;

  double x,y;

  std::random_device rd;
  std::mt19937 gen( rd());
  uniform_int_distribution<int> dist(0,99);
  for (i=0; i<order; i++) {
    x=x0;
    sig.push_back(x);
  }
  for (i=order; i<n; i++) {
    x=0.0;
    for (j=0; j<order; j++) {
       y=sig.at(i-j-1);
       x+=ar[j]*y;  // fred: CORRECTION VITE
    }
    x=x + (dist(gen)-50)*fvar/100.0;
    sig.push_back(x);
  }
  for (auto it=sig.begin(); it != sig.end(); it++) {
     std::cout << " " << (*it);
  }
  std::cout<< "\n";  

  return 1;
}


int generate_from_file(std::vector<double> &sig, int &n,
	       	const char* fileName)
{
   li_doubles li_reals;
  // use the method of the related class IODATA declared in ../include/ident05_fftdata.hpp
   IDENT05_IODATA instIODATA = IDENT05_IODATA(100, 1, fileName) ;
   instIODATA.read_extern_output(2, li_reals);
   std::cout << "size " << li_reals.size() << "\n";  // DEBUG
  // recopy .. it seems not practical to have a vector instead a li_doubles however
   n=1;
   for (auto it=li_reals.begin(); it != li_reals.end(); it++) {
      std::pair<double, double> singleton = *it;
      std::cout << singleton.second << " : " << singleton.second << "\n";
      sig.push_back(singleton.second);

      n++;
   }
   return 1;
}

int generate_more_from_file(std::vector<double> &sig, int &Npty,
	       	const char* fileName, int n, int m)
{
//   std::pair<double,double> singleton; // TBRD
   //double line_reals[8]; // maximum 8
   //std::cout << fileName << "\n";  // DEBUG
   std::ifstream Inpfile(fileName);
   std::istream_iterator<double> my_it(Inpfile);
   int k=1;
   for ( ; my_it != std::istream_iterator<double>(); my_it++)
   {
      //std::cout << "it at " << k << "\t" << *my_it << "\n";   
      sig.push_back(*my_it);  // m-1 means m-th from 0..
      k++;
      //if (k >= Npty) {
	//      break;
      //}
   }
   std::cout << "exit procedure read more from file\n";
   Inpfile.close();
   return 1;
}

int generate_more_from_file_second(std::vector<double> &sig, int dsize, const char* fileName, int m) {

    // assume the file has exactly three columns, and more readable way
    std::string line;
    std::ifstream Inpfile(fileName);
    double* linedou = new double[3];
    while (getline(Inpfile, line)) {
	if (line == "") 
		break;
        std::stringstream ssline(line);
        ssline >> linedou[0] >> linedou[1] >> linedou[2];
	sig.push_back(linedou[0]); // we are in C++, index begins by 0
    }
    Inpfile.close();

    return 1;
}

