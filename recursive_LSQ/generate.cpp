#include <iostream>
//#include <generator>
#include <random>
//#include <vector>

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
