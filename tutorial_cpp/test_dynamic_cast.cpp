#include <iostream>

struct A
{
    virtual void f(void)
    {
        return ;
    }
};

struct B : virtual public A
{
};

struct C : virtual public A, public B
{
};

struct D
{
    virtual void g(void)
    {
        return ;
    }
};

struct E : public B, public C, public D
{
};

int main(void)
{
    E e;        // e contient deux sous-objets de classe B
                // (mais un seul sous-objet de classe A).
                // Les sous-objets de classe C et D sont
                // fr�res.
    A *pA=&e;   // D�rivation l�gale : le sous-objet
                // de classe A est unique.
    C *pC=(C *) pA;// Ill�gal : A est une classe de base
                      // virtuelle (erreur de compilation).
    std::cout << "Normally we raise an error" << std::endl;
    //C *pC=dynamic_cast<C *>(pA);  // L�gal. Transtypage
                                  // dynamique vertical.
    D *pD=dynamic_cast<D *>(pC);  // L�gal. Transtypage
                                  // dynamique horizontal.
    B *pB=dynamic_cast<B *>(pA);  // L�gal, mais �chouera
                                  // � l'ex�cution (ambigu�t�).
    std::cout << "We are at the end" << std::endl;
    return 0 ;
}