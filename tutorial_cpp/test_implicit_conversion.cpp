#include <iostream>

class chaine
{
    size_t Taille;
    char * s;

public:
    chaine(void);
    // Ce constructeur permet de préciser la taille de la chaîne
	// à sa création :
    explicit chaine(unsigned int);
    ~chaine(void);
};

int main(int argc, char **argv) {
//Avec cette declaration, l expression suivante :

//int j=2;
//chaine s = j;
//n'est plus valide, alors quelle l etait lorsque le constructeur n etait pas declare explicit.

//On prendra garde au fait que le mot cle explicit n'empeche l utilisation du constructeur dans les operations de transtypage que dans les conversions implicites. 
int j=2;
chaine s = (chaine) j;

return 0;
}

