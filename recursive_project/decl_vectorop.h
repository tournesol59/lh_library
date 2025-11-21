/* ------- declaration of an own precomplete class vector as basis container ---------*/
class vectorprec {
   public:
     vectorprec(int dsize);
     vectorprec(const vectorprec &);
     vectorprec &operator=(const vectorprec &);
     ~vectorprec(); 
   public:  // oblige de mettre public pour avoir acces en C, is there any other opt?
     int size;
     double *values;
};

/* ------- declaration of an own class vector (just three constructors and at, pushback ---------*/

class vectorop : public vectorprec {
   public:
     vectorop(int dsize);
     vectorop(const vectorop &);
     vectorop &operator=(const vectorop &);
     ~vectorop(); 
     double &operator() (int i);
     double operator() (int i) const;
     double &at(int i);
     double at(int i) const;
     void push_back(double val);
   // TBD encore: scal_prod and add
     double scal_prod(const vectorop & );
     vectorop &operator+=(const vectorop & ); 
     void addvector(const vectorop & );
     vectorop &operator-=(const vectorop & );
     void decvector(const vectorop & );
     void multscal(double scalar);

   public:  // oblige de mettre public pour avoir acces en C, is there any other opt?
     int size;
     double *values;
};

