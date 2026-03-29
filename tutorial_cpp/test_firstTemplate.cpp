#include <iostream>
#include <stdio.h>

// nombre int modulo un autre nombre
class argroup {
   public:
      argroup(int dmod);
      argroup(int dmod, int di);
      ~argroup(void);
      int getimod();
      int getorder();

      bool operator==(const argroup &B) const;
      bool operator<(const argroup &B) const;
      /*
      argroup operator+=(const argroup &B);
      argroup operator-=(const argroup &B);
      argroup operator*=(const argroup &B);
*/
   private:
      int imod;
      int order;
};

argroup::argroup(int dord) :
   order(dord) {}

argroup::argroup(int dord, int di) :
   imod(di),
   order(dord) 
{
   imod=imod % order;
   if (imod<0) imod=imod+order;
 }

argroup::~argroup() 
{}

int argroup::getimod() {
   return imod;
}

int argroup::getorder() {
   return order;
}

// equality comparator
bool argroup::operator==(const argroup &B ) const
{//error: passing 'const argroup' as 'this' argument discards qualifiers [-fpermissive]
 // if (order!=B.getorder())
   if (order!=B.order) {
      std::cout << "groups have not the same order" << std::endl;
   }
   return ((imod % order)==(B.imod % B.order));
}

// comparator
bool argroup::operator<(const argroup &B ) const
{
   if (order!=B.order) {
      std::cout << "groups have not the same order" << std::endl;
   }
   return ((imod % order) < (B.imod % B.order));
}

/**
 *
// self-addition
argroup argroup::operator+=(const argroup &B ) const
{
   if (order!=B.getorder()) {
      std::cout << "groups have not the same order" << std::endl;
   }
   else {
      imod = ((imod+B.getimod()) % order);
   }
   return *this;
}

// self-substraction
argroup argroup::operator-=(const argroup &B ) const
{
   if (order!=B.getorder()) {
      std::cout << "groups have not the same order" << std::endl;
   }
   else {
    imod = ((imod-B.getimod()) % order);
    if (imod<0) imod=imod+order; // to have all elements between 0..order-1
   }
   return *this;
}

// self-multiplication
argroup argroup::operator*=(const argroup &B ) const
{
   if (order!=B.getorder()) {
      std::cout << "groups have not the same order" << std::endl;
   }
   else {
    imod = ((imod * B.getimod()) % order);
    if (imod<0) imod=imod+order; // to have all elements between 0..order-1
   }
   return *this;
}
*/

/*** function template ***/
template <class T>
T Min(T x, T y)
{
    return x<y ? x : y;
}

/*** main ***/

int main(int argc, char **argv) {
   int m,d1,d2;
   sscanf(argv[1],"%d", &m);
   sscanf(argv[2],"%d", &d1);
   sscanf(argv[3],"%d", &d2);
   
   argroup A(m,d1), B(m,d2);
   
   argroup C=Min(A,B);
   std::cout << "Resultat du min: " << C.getimod() << std::endl;
   return 0;
}
