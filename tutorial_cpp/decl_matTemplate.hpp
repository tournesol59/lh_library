/**
 * Fred: implementer les matrices pleines de taille indiquee dans un structure
 * mais le faire avec completude : produit, trans, access row,col
 */
#include <iostream>

/**
 * classe matrice comme double tableau
 */


template <class T>
class matrix {
  public:
     matrix(int nr,int mc);
     matrix(const & M);
     matrix &operator=(const & M);
     ~matrix();
     matrix operator+=(const matrix &M);
     matrix operator-=(const matrix &M);

     int nr_, mc_;
     std::vector<T> matrix_;  

};

template <class T>
matrix::matrix(int nr, int mc) :
   nr_(nr),
   mc_(mc),
   matrix_(vector<T>(nr*mc), 0.0)
{
   return;
}

template <class T>
matrix::matrix(const & M) :
   nr_(M.nr_),
   mc_(M.mc_),
   matrix_(std::vector<T> (M.nr_*M.nc_), 0.0)
{
//   nr_=M.nr_;
//   mc_=M.mc_;
   for (int i; i<nr_; i++) { // COL_MAJOR
      for (int j; j<mc_; j++) {
      matrix_[j*nr_+i]=M.matrix_[j*nr_+i];
      }
   }
   return;
}
template <class T>
matrix &matrix:: operator=(const & M) :
   nr_(M.nr_),
   mc_(M.mc_),
   matrix_(std::vector<T> (M.nr_*M.nc_), 0.0)
{
   for (int i; i<nr_; i++) { // COL_MAJOR
      for (int j; j<mc_; j++) {
      matrix_[j*nr_+i]=M.matrix_[j*nr_+i];
      }
   }
   return;
}


