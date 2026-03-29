#include <iostream>

struct Base
{
    int i;
};

struct Derivee : public Base
{
    int i;
    int LitBase(void);
};

int Derivee::LitBase(void)
{
    return Base::i; // Renvoie la valeur i de la classe de base.
}

int main(void)
{
    Derivee D;
    D.i=1;          // Accède à l'entier i de la classe Derivee.
    D.Base::i=2;    // Accède à l'entier i de la classe Base.
    std::cout << "Deriv : " << D.i << ", Base :" << D.Base::i <<  std::endl;
    return 0;
}

