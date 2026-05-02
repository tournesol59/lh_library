#include <cstdio>

class A {
public:
    virtual void Hello() = 0;
};

void A::Hello() {
    printf("A::Hello\n");
}

class B : public A {
public:
    void Hello() {
        printf("B::Hello\n");
        A::Hello();
    }
};

int main() {
    /* Prints:
           B::Hello
           A::Hello
    */
    B b;
    b.Hello();
    return 0;
}

