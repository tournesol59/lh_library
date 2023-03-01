#include "./decl_MaClasse.hpp"

using namespace lhlib;

/***
 * class MaClassBase methods
***/
MaClasseBase::MaClasseBase(std::vector<double> entrees, int dsize, int dnp) :
	valeurs(std::vector<double> (dsize, 0.0)),
        newvaleurs(std::vector<double> (dsize, 0.0)),	
	size(dsize),
	np(dnp),
        coeffs(std::vector<double> (dnp, 0.0))

{
   for (int i=0; i< dsize; i++) {
      valeurs[i]=entrees[i];
   } 
}
/*        
MaClasseBase::MaClasseBase(const MaClasseBase &source)
{
}

MaClasseBase & MaClasseBase :: operator=(const MaClasseBase &source)
{
}
*/
MaClasseBase::~MaClasseBase() 
{
}

bool MaClasseBase::entre(void) {
    for (int i=0; i<np; i++) {
       std::cout << "Base: Entrez valeur i= " << i << std::endl;
       std::cin >> coeffs[i];
    }
    return 0;
}
    
bool MaClasseBase::recalcule(void) {
    entre();
    std::cout << "Base Re-calculated values:" << std::endl;
    for (int i=0; i <size; i++) {
       newvaleurs[i] = coeffs[0];
       std::cout << newvaleurs[i] << ", ";
    }
    std::cout << std::endl;
    return 0;
}


/***
 * class MaClassDerivee methods
***/

MaClasseDerivee::MaClasseDerivee(std::vector<double> entrees, int dsize, int dnp) : MaClasseBase(entrees, dsize, dnp) 
/*
MaClasseDerivee::MaClasseDerivee(std::vector<double> entrees, int dsize, int dnp) :
	valeurs(std::vector<double> (dsize, 0.0)),
        newvaleurs(std::vector<double> (dsize, 0.0)),	
	size(dsize),
	np(dnp),
        coeffs(std::vector<double> (dnp, 0.0))
	*/
{
}
/*
MaClasseDerivee::MaClasseDerivee(const MaClasseDerivee &source) : MaClasseBase(source)
{
}

MaClasseDerivee &MaClasseDerivee :: operator=(const MaClasseDerivee &source)
{
}
*/
MaClasseDerivee::~MaClasseDerivee(void) 
{
}

bool MaClasseDerivee::entre(void) {
   bool r=MaClasseBase::entre();  // call to mother class, raises memory error at call of recalcule()
//
//    for (int i=0; i<np; i++) {
//       std::cout << "Derived Entrez valeur i= " << i << std::endl;
//       std::cin >> coeffs[i];
//    }

   std::cout << "We are in Derived Class after entre() " << std::endl;
   return r;
}


bool MaClasseDerivee::recalcule(void) 
{
   entre();  // commented
   if (size >= np) {
	   // check variables accessability
      std::cout << "List of params coeffs" << std::endl;
      for (int k=0; k < np; k++) {  // k < np
	      std::cout << coeffs[k] << ", ";  // check
      }
      std::cout << std::endl;
      std::cout << "Derived Recalculated values:" << std::endl;
      for (int i=0; i < np; i++) {
          newvaleurs[i]=valeurs[i];
      }    
      for (int i=np; i <size; i++) {
	  newvaleurs[i]=0.0;
	  for (int k=0; k < np; k++) {
              newvaleurs[i] += coeffs[k]*newvaleurs[i-k-1]; // yes -1
	  }
	  std::cout << newvaleurs[i] << ", ";
      }
      std::cout << std::endl;
      
   }
   return 0;
}

/***
 * class MaDeuxClassDerivee methods
**
MaDeuxClasseDerivee::MaDeuxClasseDerivee(std::vector<double> entrees, std::vector<double> com, int dsize, int dna, int dnb) : MaClasseDerivee(entrees, dsize, dna) 

{
   na=dna;
   nb=dnb;
   np=dna+dnb;
  
   for (int k=0; k < dnb; k++) {
      coeffs.push_back(0.0);   // with this, allocate the last nb entries (na entries have already been allocated by the base constructor
   }
   for (int k=0; k < dna+dnb; k++) {
      states.push_back(0.0);  // with this, allocate the good dimension
   }
   for (int i=0; i< dsize; i++) {
      commande.push_back(com[i]);
      valeurs[i]=entrees[i];
   }
  
}

// copy constructor
MaDeuxClasseDerivee::MaDeuxClasseDerivee(const MaDeuxClasseDerivee &source) : MaClasseDerivee(source)
{

   size=source.size;
  na=source.na;
  nb=source.nb;
   np=source.np;
   for (int k=0; k < source.np; k++) {
      coeffs.push_back(source.coeffs[k]);   // no previous alloc by the base constructor, since we redefine the copy method here
      states.push_back(source.states[k]); 
   }
   for (int i=0; i< source.size; i++) {
      valeurs.push_back(source.valeurs[i]);
      newvaleurs.push_back(source.newvaleurs[i]);
      commande.push_back(source.commande[i]);
   }
  
}


//affect constructor
MaDeuxClasseDerivee &MaDeuxClasseDerivee :: operator=(const MaDeuxClasseDerivee &source)
{

   size=source.size;
  // na=source.na;
   //nb=source.nb;
   np=source.np;
   for (int k=0; k < source.np; k++) {
      coeffs.push_back(source.coeffs[k]);   // no alloc by the base constructor
      states.push_back(source.states[k]); 
   }
   for (int i=0; i< source.size; i++) {
      valeurs.push_back(source.valeurs[i]);
      newvaleurs.push_back(source.newvaleurs[i]);
      commande.push_back(source.commande[i]);
   }
 
}

// destructor
MaDeuxClasseDerivee::~MaDeuxClasseDerivee(void) 
{
}

bool MaDeuxClasseDerivee::entre(void) {
   bool r=MaClasseBase::entre();  // call to mother class, raises memory error at call of recalcule()

   std::cout << "We are in Deux Derived Class after entre() " << std::endl;
   return r;
}

bool MaDeuxClasseDerivee::recalcule(void) 
{
   entre(); 
   if (size >= np) {
	   // check variables accessability
      std::cout << "List of params coeffs" << std::endl;
      for (int k=0; k < np; k++) {  // k < np
	      std::cout << coeffs[k] << ", ";  // check
      }
      std::cout << std::endl;
          // for each sample, calculate predicted value
      std::cout << "Derived Deux Recalculated values:" << std::endl;
      for (int i=0; i < na; i++) {
          newvaleurs[i]=valeurs[i];
      }    
      for (int i=na; i <size; i++) {
	  newvaleurs[i]=0.0;
          for (int k=0; k < na; k++) {  // k < na
	     states[k] = valeurs[i-k-1];
          }
	  // ATTENTION CHANGE na to nb
          for (int k=0; k < nb; k++) {  // k < nb
	     states[na+k] = commande[i-k-1];
          }		  
	  for (int k=0; k < np; k++) { // k < np = na + nb
             newvaleurs[i] += coeffs[k]*states[k]; // yes -1
	  }
	  std::cout << newvaleurs[i] << ", ";
      }
      std::cout << std::endl;
   }
   return 0;
}
*/

/***
 * class MaTroisClassDerivee methods

MaTroisClasseDerivee::MaTroisClasseDerivee(std::vector<double> entrees, std::vector<double> com, int dsize, int dna, int dnb, int idelay) : MaDeuxClasseDerivee(entrees, com, dsize, dna, dnb) 
{
   delay=idelay;
}	

//copy constructor
MaTroisClasseDerivee::MaTroisClasseDerivee(const MaTroisClasseDerivee &source) : MaDeuxClasseDerivee(source)
{
   delay=source.delay;
}

MaTroisClasseDerivee &MaTroisClasseDerivee :: operator=(const MaTroisClasseDerivee &source)
{
   delay=source.delay;
}

// destructor
MaTroisClasseDerivee::~MaTroisClasseDerivee(void) 
{

}

bool MaTroisClasseDerivee::recalcule(void) 
{
   entre();  // call to virtal method, hence of the current class
   
   for (int i=0; i < na; i++) {
       newvaleurs[i]=valeurs[i];
   }    
   for (int i=na; i < size; i++) {
       newvaleurs[i]=0.0;
       for (int k=0; k < na; k++) {
          states[k] = valeurs[i-k-1];
       }
   // ATTENTION CHANGE na to nb
       for (int k=0; k < nb; k++) {  // k < nb
          states[na+k] = commande[i-k-delay];  // important assumption: delay < na
       }		  
       for (int k=0; k < np; k++) { // k < np = na + nb
          newvaleurs[i] += coeffs[k]*states[k]; // yes -1
        }
	std::cout << newvaleurs[i] << ", ";
      }
      std::cout << std::endl;

   return 0;
}
*/

