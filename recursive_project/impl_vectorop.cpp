/* ------- declaration of an own class vector (just three constructors and at, pushback ---------*/
#include "decl_vectorop.h"

vectorprec::vectorprec(int dsize) : size(dsize) {
   values = new double[dsize];
}

vectorop::~vectorop() {
}

/*-------------------------------*/
vectorop::vectorop(int dsize) : vectorprec(dsize) {
   values = new double[dsize];
}

vectorop::vectorop(const vectorop & v) : vectorprec(v) {
   size=v.size;
   values = new double[size];
   for (int i=0; i<size; i++) {
       values[i]=v.values[i];
   }
}

vectorop &vectorop::operator=(const vectorop & v) {
  if (size!=v.size) {
   size=v.size;
   delete values;
   values = new double[size];
  }
  for (int i=0; i<size; i++) {
     values[i]=v.values[i];
  }
}

vectorop::~vectorop() {
}

double &vectorop::operator() (int i) {
   return values[i];
}

// or equivalently the "at" submethod
double &vectorop::at(int i) {
   return values[i];
}

double vectorop::operator() (int i) const {
   return values[i];
}
// or equivalently the "at" submethod
double vectorop::at(int i) {
   return values[i];
}

void vectorop::push_back(double val) {
   size=size+1;
   double* temp=new double[size];
   temp[0]=val;
   for (int i=1; i<size; i++) {
      temp[i]=values[i-1];
   }
   delete values;
   double* values = temp;
}

double vectorop:: scal_prod(const vectorop &bvec) {
  int sizemin=(size<bvec.size? size : bvec.size);
  double sum=0.;
  for (int i=0; i<size; i++) {
     sum += values[i]*bvec.values[i];
  }
  return sum;
}

void vectorop::addvector(const vectorop & bvec) {
   if (size==bvec.size) {
      for (int i=0; i<size; i++) {
         values[i] += bvec.values[i];
      }
   }
}

vectorop &vectorop::operator+=(const vectorop & bvec) {
   this->addvector(bvec);
   return *this;
}

void vectorop::decvector(const vectorop & bvec) {
   if (size==bvec.size) {
      for (int i=0; i<size; i++) {
         values[i] -= bvec.values[i];
      }
   }
}

vectorop &vectorop::operator-=(const vectorop & bvec) {
   this->decvector(bvec);
   return *this;
}

void vectorop::multscal(double scalar) {
   for (int i=0; i<size; i++) {
      values[i] *= scalar;
   }
}


