#include <iostream>
#include <vector>
#include <string.h>
#include <cstring>
#include <sstream>
#include "../include/lhTypes.hpp"

using namespace lhlib;

// demonstration d'une classe virtuelle


/*********** la classe de Base **********/
class MonObjet {
   unsigned long int new_handle(void);


   public:
     MonObjet(void);
     MonObjet(int dtaille, char *s);
     virtual ~MonObjet(void);   // the virtual destructor
     virtual void print(void) const = 0;  // the virtual method pure
     virtual unsigned long int handle(void) const =0; // function returning h
                                         // wh. is the object serial Nr
     virtual char* content(void) const =0; // function returning string "name"
     
   //protected: // do not work if not this uncommented
     unsigned long int h; 
     int taille;
     char *name;

};

/********** la classe MonContainer *********/
class MonContainer : public MonObjet  //  Container class, wh. inherits
				// of MonObjet, because a Container can 
				// contain another Cpontainer. The container 
				// is implemented with linked list
{
   struct MonObjetList {
	MonObjetList *next;
        MonObjet *ptr;
   };

   MonObjetList *head;

   public:
     MonContainer(void);        // The contructor : calls those of Object.
     ~MonContainer(void);       // The destructor.
     void print(void) const; // display of the linked list content
     unsigned long int handle(void) const; // function returning h
                                         // wh. is the object serial Nr
     char* content(void) const; // retuning the name

     bool has(unsigned long int h) const;  // returns true if Container contains
                              // an Object of serial = h
     bool is_empty(void) const;   // true if Container is empty.
     void add(MonObjet &);          // Add an object of type 'MonObjet'.
     void remove(MonObjet &);       // Remove an object.
};

/********** la classe derivant de la classe abstraite MonObjet: MonObjetDeriv
***/
class MonObjetDeriv : public MonObjet {
   public:
     MonObjetDeriv(void);
     MonObjetDeriv(int dtaille, char *s);
     ~MonObjetDeriv(void);  	   
     void print(void) const;
     unsigned long int handle(void) const;
     char* content(void) const;

};

