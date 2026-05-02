// Class.cpp : définit le point d'entrée pour l'application console.
//
 
//#include "stdafx.h"
#include <iostream>
 
using namespace std;
 
class moteur	// classe moteur
{
public:
	moteur(int nombrecv);
	~moteur();
 
private:
	int sonmoteur;
};
 
 
moteur::moteur(int nombrecv)	// Constructeur moteur	
{
	sonmoteur = nombrecv;
}
 
moteur::~moteur()						// destructeur moteur
{
 
}

 
class voiture
{
public:
	voiture(int vitesseinitiale, int nombredecv);				// déclaration Constructeur
	~voiture();									// déclaration Destructeur
	int lirevitesse();							// déclaration des méthodes/fonctions de la classe
	void stopper();
	void accelerer();
 
    moteur deuxchevaux;           // c'est ici que se situe l'erreur
 
private:
	int vitesse;
};
 
voiture::voiture(int vitesseinitiale, int nombredecv) : deuxchevaux(nombredecv)	// Constructeur
{
	vitesse = vitesseinitiale;
 
}
voiture::~voiture()								//Destructeur 
{
	cout << "appel destructeur" << endl;
}
 
 
 
int voiture::lirevitesse()
{
	cout << "vitesse : " << vitesse << endl;
	return vitesse;
}
 
void voiture::stopper()
{
	vitesse = 0;
}
 
void voiture::accelerer()
{
	vitesse = 100;
}
 
 
int main()
{
	//	moteur test(2);              fonctionne

	voiture peugeot(0,2);
	peugeot.lirevitesse();
 
	peugeot.accelerer();
	peugeot.lirevitesse();
	peugeot.stopper();
	peugeot.lirevitesse();
	peugeot.~voiture();
	peugeot.lirevitesse();
 
 
	return 0;
}
