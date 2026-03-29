#include <iostream>

class A {
    public:
	    static A & getInstance();
    private:
	    A();
	    A(const &ac);
	    A &operator=(const & as);
};

A::A() {}
A &A::getInstance() {
    static A a;
    return a;
}

int main() {
   // A err;  // doit lever une erreur
   A & s = A::getInstance();
   std::cout << "Objet A cree" << std::endl;
}
