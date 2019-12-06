#include <iostream>
//#include <generator>
#include <random>
//#include <vector>

using namespace std;

int generate_rdom_example( vector<double> &sig, int n, double ca, double cb) {

  int i;
  double mean=0.0;
  double stddev = 0.5;
  double x;

//  std::default_random_engine generator;
//  std::normal_distribution<double> dist(mean, stddev);
  std::random_device rd;
  std::mt19937 mt( rd());
  uiform_int_distribution<int> dist(0,99);
  for (i=0; i<n; i++) {
    x=ca+cb*1*i + (dist(generator)-50)*1.0/100;
    sig.push_back(x);
  }

  std::copy(begin(sig), end(sig), std::ostream_iterator<double>(std::cout, " "));
  std::cout<< "\n";

  return 1;
}
