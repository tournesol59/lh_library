#include <iostream>
#include <cstring>

/**
 * Class Declaration
 */

class A
{
public:
   A(int d);
   A(const A &ac);
   A &operator=(const A &);
   int getDonnee();
protected:
   int Donnee;
};

// Heritage de la classe A, virtuelle
class B : virtual public A
{
public:
   B(int d);
   B(int d, int db);
   B(const B &bc);
   B &operator=(const B &);

   protected:
      int Valeur_B;  // autre que "Donnee"
};

// A est toujours virtuelle
class C : virtual public A
{
public:
   C(int d);
   C(int d, int dc);
   C(const C &cc);
   C &operator=(const C &);

   protected:
      int Valeur_C;  // autre que "Donnee"
                   // "Donnee" est acquise par heritage
};

class D : public B, public C  // Ici, Donnee n est pas duplique
{
public:
   D(int d);
   D(const D &dc);
   D &operator=(const D &);

   // Definition de la classe D
};

/**
 * Implementation of classes A,B,C,D
 */

A :: A(int d) :
   Donnee(d)
{
   return;
}

A :: A(const A &ac) {
   Donnee= ac.Donnee;
   return;
}

A &A :: operator=(const A &as) 
{
   Donnee = as.Donnee;
   // return;
}

A :: getDonnee() 
{
   return Donnee;
}

/* ----  B ----  */
B :: B(int d) : A(d) 
{
   return;
}

B :: B(int d, int db) : A(d) 
{
   Valeur_B = db;
   return;
}

B :: B(const B &bc) : A(bc)
{
   Donnee = bc.Donnee;
   Valeur_B = bc.Valeur_B;
   return;
}

B &B :: operator=(const B &bs) 
{
   Donnee = bs.Donnee;
   Valeur_B = bs.Valeur_B;
   // return;
}

/* ---- C ---- */
C :: C(int d) : A(d) 
{
   return;
}

C :: C(int d, int dc) : A(d) 
{
   Valeur_C = dc;
   return;
}

C :: C(const C &cc) : A(cc)
{
   Donnee = cc.Donnee;
   Valeur_C = cc.Valeur_C;
   return;
}

C &C :: operator=(const C &cs) 
{
   Donnee = cs.Donnee;
   Valeur_C = cs.Valeur_C;
   //return;
}

/* ---- D ---- */

D :: D(int d) : B(d), C(d), A(d)
{
   return;
}

D :: D(const D &dc) : B(dc), C(dc), A(dc)
{
   return;
}

D &D :: operator=(const D &ds) // : B(ds), C(ds), A(ds)
{
  Donnee = ds.Donnee;
  Valeur_B = ds.Valeur_B;
  Valeur_C = ds.Valeur_C;
 //  return;
}

/**
 * Main
 */

int main(int argc, char **argv) {

   A *pa;
   B *pb,*pba,*pbd;
   C *pc;
   D *pd;

   A ia= A(1);
   B ib= B(2,2);
   C ic= C(3,3);
   D id= D(4);
   std::cout << "Donnees de la classe D: " << id.getDonnee() << std::endl;

   // test pointers
   pa = &ia;
   pb = &ib;
   std::cout << "Donnees de la classe B par ptr pb: " << pb->getDonnee() << std::endl;
   pd = &id;
   //pba = (A*) pb; // we cannot do this, since B inherits virtual A
//   pbd = (B*) pd; // can we?
   pbd = static_cast<D *>(pb);
   std::cout << "Donnees de la classe B par ptr pd: " << pbd->getDonnee() << std::endl;
   std::cout << "We are at the end" << std::endl;
   return 0;
}
