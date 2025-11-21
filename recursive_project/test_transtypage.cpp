#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

int main() {
	string g {"Green is my favorite color"};
	cout<<g.substr(0,5)<<endl;
	cout<<g.substr(6,2)<<endl;
	cout<<g.substr(21,5)<<endl;
	cout<<endl; return 0;
}
