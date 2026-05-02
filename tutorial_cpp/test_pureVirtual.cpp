#include <iostream>

using namespace std;

class Base
{
    public:
        virtual void Sub1() = 0;
        virtual void Sub2();
        virtual void Sub3();
        void Sub4();
};

class A : public Base
{
    public:
        void Sub2();
        void Sub4();
};

class B : public A
{
    public:
        virtual void Sub1();
        void Sub2();
};

class C : public Base
{
    public:
        virtual void Sub1();
        virtual void Sub4();
};

void Base::Sub2()
{
    cout << "Hello from Base::Sub2()" << endl;
}

void Base::Sub3()
{
    cout << "Hello from Base::Sub3()" << endl;
    Sub1(); // DONT MISS THIS CALL IN YOUR ANSWER
}

void Base::Sub4()
{
    cout << "Hello from Base::Sub4()" << endl;
}

void A::Sub2()
{
    cout << "Hello from A:Sub2()" << endl;
}
void A::Sub4()
{
    cout << "Hello from A:Sub4()" << endl;
}

void B::Sub1()
{
    cout << "Hello from B:Sub1()" << endl;
}
void B::Sub2()
{
    cout << "Hello from B:Sub2()" << endl;
}

void C::Sub1()
{
    cout << "Hello from C:Sub1()" << endl;
}
void C::Sub4()
{
    cout << "Hello from C:Sub4()" << endl; //error used to say from Sub2
}

void Sub(Base& x)
{
    x.Sub1();
    x.Sub2();
    x.Sub3();
}
void AnotherSub(A& a)
{
    a.Sub1();
    a.Sub2();
    a.Sub4();
}

int main()
{
    A a; // wont compile
    B b;
    C c;
    Sub(a);
    Sub(b);
    Sub(c);
    AnotherSub(b);
}
