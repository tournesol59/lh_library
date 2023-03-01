#include <iostream>

class A {
    public:
	    static A &getInstance(void);
    private:
	    A();
	    A(const &ac);
	    A &operator=(const & as);
};

A :: A() {}

A &A::getInstance() {
    static bool active=0;
    static A a;
    if (active==0) {
        a=A();
	active=1;
    }
    return a;
}

int main() {
//   A a=A(j);  // doit lever une erreur
   A &s = A::getInstance();
   std::cout << "Objet A cree" << std::endl;
}
