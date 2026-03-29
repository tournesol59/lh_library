#include "decl_MonContainer.hpp"

// cette fonction n'est appelable que par la class MonObjet

unsigned long int MonObjet::new_handle(void)
{
   static unsigned long int hc = 0; 
   return hc = hc + 1;
}

// le constructeur de MonObjet doit etre appele par les classes derivees
// plusieurs constructeurs possibles
MonObjet :: MonObjet() 
{ 
   h = new_handle();
   taille = 0;
   return;
}

MonObjet :: MonObjet(int dtaille, char * s)
{
   h = new_handle();  
   taille = dtaille;
   strncpy(name, s, dtaille);

   return;
}

MonObjet::~MonObjet(void) 
{
  return;
}

/********* Methodes de MonContainer ********/
// Le constructeur : a
MonContainer :: MonContainer(void) : MonObjet()
{
// par defaut
  return;
}


MonContainer :: ~MonContainer(void)
{
  MonObjetList *tmp = head;
  while (tmp != NULL) {
     tmp = tmp -> next;
     delete head;
     head = tmp;  //astuce pour n'utiliser qu'une fois delete dans cette boucle
  }

  return;
}

void MonContainer :: print(void) const // Fonction d'affichage du Container
{

   MonObjetList *tmp1 = head;
   char *s1;
   const char *c1 = (const char*) content();
   strncpy(s1, c1, taille);

   std::cout << "Container n°" << handle() << ", has name " << s1 << std::endl; // do not use tmp-> attributes
   std::cout << "\t";
   while (tmp1->next != NULL) {
      tmp1 = tmp1->next;
      tmp1->ptr->print(); // do not use tmp-> attributes
   }
    return;	
}

unsigned long int MonContainer :: handle(void) const // function returning h
                                         // wh. is the object serial Nr
{
   return head->ptr->h;
}

char*  MonContainer :: content(void) const // retuning the name
{
   return head->ptr->name;
}

     
bool MonContainer :: has(unsigned long int h) const  // ne modifie pas le Container
{
   MonObjetList *tmp = head;
   while (tmp != NULL) {
      if (tmp->ptr->handle() == h)   // cherche l'objet
         break;
   }
   return (tmp != NULL);
}

bool MonContainer :: is_empty(void) const   // true si le Container est vide.
{
   return (head == NULL);
}

void MonContainer :: add(MonObjet &source)  // Ajoute un objet en tete de liste
{
   MonObjetList *tmp = new MonObjetList();
   tmp->ptr = &source;
   tmp->next = head;
   head = tmp;
   return;
}

void MonContainer :: remove(MonObjet &source) // Retire un objet correspondant
{
   MonObjetList *tmp1, *tmp2; 
   tmp1 = head;

   if (tmp1->next != NULL) {
      tmp2=tmp1->next;
      while ((tmp2->next != NULL) && (tmp2->ptr->handle() != source.handle())) {
         tmp1=tmp2;  //cherche l'objet
	 tmp2=tmp2->next;
      }
      // l'a t on trouve?
      if (tmp2 != NULL) {
         tmp1->next = tmp2->next; // il fallait garder un pointeur tmp1 sur l'objet avant pour cette operation
	 delete tmp2;
      }
      else {
         // rien a faire
      }
   }
   else {  // cas particulier de l'objet en tete de liste
      head = head-> next;
      delete tmp1;
   }
   return;
}

/********** la classe derivant de la classe abstraite MonObjet: MonObjetDeriv
***/
MonObjetDeriv :: MonObjetDeriv(void) : MonObjet()
{
  return;       
}

MonObjetDeriv :: MonObjetDeriv(int taille, char* s) : MonObjet(taille,s)
{
  return;       
}

MonObjetDeriv :: ~MonObjetDeriv(void) 
{
   return;
}

void MonObjetDeriv::print(void) const 
{
  std::cout<< "Object Nr :" << h << " has name :" << name << std::endl;
}

unsigned long int MonObjetDeriv::handle(void) const {
   return h;
}

char * MonObjetDeriv::content(void) const {
   return name;
}

