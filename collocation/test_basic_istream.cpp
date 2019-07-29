#include <iostream>
#include <sstream>
#include <vector>
#include <array>
 
int main()
{
    std::istringstream input("INT:INIEQN:END");
    std::vector<std::array<char, 8>> listchar;
 
    // note: the following loop terminates when std::ios_base::operator bool()
    // on the stream returned from getline() returns false
    for (std::array<char, 8> a; input.getline(&a[0], 8, ':'); ) {
        listchar.push_back(a);
    }
 
    for (auto& a :listchar) {
        std::cout << &a[0] << '\n';
    }
}

