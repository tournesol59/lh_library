#include <iostream>

#define MAX_NUMBER  20

template <class T>
class Tableau
{
  // Definition de la classe template tableau
   public:
     Tableau(int nc, const char *initdata);
//     Tableau(const &Tab);
//     Tableau &operator=(const &Tab);
     ~Tableau();
     getI(int i);
   private:
     int ncol;
     T m_tableau[MAX_NUMBER];
};

/*** implementation de la classe tableau ***/
template <class T>
Tableau<T>::Tableau(int nc, const char *initdata) :
	ncol(nc)
{
   double item;
   //m_tableau = new T[MAX_NUMBER];
   for (int i; i<nc; i++) {
      sscanf((initdata+i*sizeof(char)),"%f", &item);
      // if item would not be of type double..
      m_tableau[i]=(T) item; // that not very in the sense of template
      			     // but the aim is testing tableau
   }
}

template <class T>
Tableau<T>::~Tableau()
{
   delete m_tableau;
}


/*** main program ***/

int main()  
{
   char initdata[5];  // we need this trick (not very Cpp) to copy data from a priori unknown type
   for (int i=0; i<5; i++) {
      sprintf((initdata+i*sizeof(char)), "%f", (double) i);
   }
   Tableau<double> tab_A(5, (const char*) initdata); 	
   return 0;
}

