/* Global Scope: Parameter Estimation, ML
 * Scope:
 * Solve a second-order linear (for the moment) equation
 * coefficients could be a polynom up to degree two.
 * The collocation method (with Chebyshev interpolation)
 * and Lanczos method is used
 * Ref: Wright (60's), Lanczos (1938), Clenshaw (50's)
 *
 * LIEBHERR TOULOUSE
 *******************************************************/
#include "ident05_coll.hpp"

#ifndef __TEST_COLL_ONLY__
#define __TEST_COLL_ONLY__
#endif

//using namespace std;
using namespace lhlib;

/** attemps to access optimization parameters via macros */

//#define P1_EPS (ident05_nlp::get_p1_eps)
//#define P2_OM (ident05_nlp::get_p2_om)

// constructor
IDENT05_COLL::IDENT05_COLL(Index typode, Index ord, const char * FileNameInp )
{
/*eqn Problem type: */
  type_eqn=typode;
  equ1[0]=16.0; equ1[1]=0.5; equ1[2]=1.0;  //  (y"+0.5*y'+16y=0)
  equ2[0]=0.0;equ2[1]=1.0;equ2[2]=0.0;  equ2[3]=1.0;equ2[4]=0.0;equ2[5]=0.0; equ2[6]=0.0;equ2[7]=1.0;equ2[8]=0.0;  //  (xy"+y'+xy=0) Bessel
 boundry[0]=0.0;
 boundry[1]=1.0;

/* Input File Name */
  lenFileNameInp = 8;  // do not change this value
  strncpy(strFileNameInp, FileNameInp, lenFileNameInp);   // works on MInGW not on Linux
  std::cout << "Have input file name: " << strFileNameInp << "\n";
/*  t0=(Number)malloc(7*sizeof(Number)); */
  t0[0]=1.0;t0[1]=0.0;t0[2]=0.0;t0[3]=0.0;t0[4]=0.0;t0[5]=0.0;t0[6]=0.0;
  t1[0]=0.0;t1[1]=1.0;t1[2]=0.0;t1[3]=0.0;t1[4]=0.0;t1[5]=0.0;t1[6]=0.0;
  t2[0]=-1.0;t2[1]=0.0;t2[2]=2.0;t2[3]=0.0;t2[4]=0.0;t2[5]=0.0;t2[6]=0.0;
  t3[0]=0.0;t3[1]=-3.0;t3[2]=0.0;t3[3]=4.0;t3[4]=0.0;t3[5]=0.0;t3[6]=0.0;
  t4[0]=1.0;t4[1]=0.0;t4[2]=-8.0;t4[3]=0.0;t4[4]=8.0;t4[5]=0.0;t4[6]=0.0;
  t5[0]=0.0;t5[1]=5.0;t5[2]=0.0;t5[3]=-20.0;t5[4]=0.0;t5[5]=16.0;t5[6]=0.0;
  t6[0]=-1.0;t6[1]=0.0;t6[2]=18.0;t6[3]=0.0;t6[4]=-48.0;t6[5]=0.0;t6[6]=32.0;

  t7[0]=0.0; t7[1]=-7.0; t7[2]=0.0; t7[3]=56.0; t7[4]=0.0; t7[5]=-112.0; t7[6]=0.0; t7[7]=64.0;
  order=ord;
  num_ranges=1;
  num_total_coeffs=(order+1)*num_ranges;

}

IDENT05_COLL::~IDENT05_COLL() 
{
}
bool IDENT05_COLL::read_parse_file() 
{
   Index row,col;
 //  char ReadC[39];  // 39 seems to be the largest character's number in a line

// Open the input file, which is a text file with 3 columns
  std::ifstream lh_file;
  lh_file.open("TheFile", std::ifstream::in);  // test w/o variable Name
  
  row=0;
  while (!lh_file.eof()) {
//  while (lh_file.good()) {
      std::string line; 
      std::getline(lh_file, line);    
    
      std::stringstream ss(line);
      col = 0;
      while(ss >> dataarray[row][col]) col++;
#ifdef __COLL_TEST_ONLY__
      std::cout << "tf= "  << dataarray[row][0] << " ";
      std::cout << "utf= " << dataarray[row][1] << " ";
      std::cout << "ytf= " << dataarray[row][2] << "\n";
#endif
      row++;
  }
   lh_file.close();

   return 0;
}

bool IDENT05_COLL::ExpandSeriesLinearSys_ref1()
{
   Index i,j,k,trunc_order;
   Number Hk[7][7];
   Number Fk[7][7];
   Number Ek[7][7];
   Number Rk[7][2];
   Number Kk[7][7];
   Number Mk[7][7];

// step1: create system: H*Di+F*Bi+E*Ai=0
    for (i=0;i<=order;i++) { //i is corresp. to terms X^i
       //corresp. to coeff Taui (order) in Chebyshev Sum of contins term 
       Rk[i][0]=t5[i];
       Rk[i][1]=t6[i];
      
    }
  //terms of matrix Ek corresponds to continuous term
   for (i=0;i<=order;i++) { //i is corresp. to terms X^i
      for (j=0;j<=order;j++) { //j corresp. to coeff Ai (order) in Chebyshev Sum of contins term 
         if (i>j) {
           Ek[i][j]=0.0;
         }
         else {
             switch(j) {
               case 0: 
			Ek[i][j]=t0[i]*equ1[0];
			break;
               case 1: 
			Ek[i][j]=t1[i]*equ1[0];
			break;
               case 2: 
			Ek[i][j]=t2[i]*equ1[0];
			break;
               case 3: 
			Ek[i][j]=t3[i]*equ1[0];
			break;
               case 4: 
			Ek[i][j]=t4[i]*equ1[0];
			break;
               case 5: 
			Ek[i][j]=t5[i]*equ1[0];
			break;
               case 6: 
			Ek[i][j]=t6[i]*equ1[0];
			break;
            }
         }//endif 
      }//end for j
   }//end for i

   //terms of matrix Fk corresponds to continuous term
   for (i=0;i<=order;i++) { //i is corresp. to terms X^i
      for (j=0;j<=order-1;j++) { //j corresp. to coeff Bi (order-1) in Chebyshev Sum of derivative term 
         if (i>j) {
           Fk[i][j]=0.0;
         }
         else {
             switch(j) {
               case 0: 
			Fk[i][j]=t0[i]*equ1[1];
			break;
               case 1: 
			Fk[i][j]=t1[i]*equ1[1];
			break;
               case 2: 
			Fk[i][j]=t2[i]*equ1[1];
			break;
               case 3: 
			Fk[i][j]=t3[i]*equ1[1];
			break;
               case 4: 
			Fk[i][j]=t4[i]*equ1[1];
			break;
               case 5:  
			Fk[i][j]=t5[i]*equ1[1];
			break;
               case 6: 
                        Fk[i][j]=0.0; 
			break;  
           }
         }//endif 
      }//end for j
      Fk[i][order]=0.0;  // fill the unused part of the matrix with zero
   }//end for i

   //terms of matrix Hk corresponds to continuous term
   for (i=0;i<=order;i++) { //i is corresp. to terms X^i
      for (j=0;j<=order-1;j++) { //j corresp. to coeff Di (order-2) in Chebyshev Sum of 2nd derivative term 
         if (i>j) {
           Hk[i][j]=0.0;
         }
         else {
             switch(j) {
               case 0: 
			Hk[i][j]=t0[i]*equ1[2];
			break;
               case 1: 
			Hk[i][j]=t1[i]*equ1[2];
			break;
               case 2: 
			Hk[i][j]=t2[i]*equ1[2];
			break;
               case 3: 
			Hk[i][j]=t3[i]*equ1[2];
			break;
               case 4: 
			Hk[i][j]=t4[i]*equ1[2];
			break;
          //     case 5: Hk[i][j]=t5[i]*equ1[2];
          //     case 6: Hk[i][j]=t6[i]*equ1[2];
             }
         }//endif 
      }//end for j
      Hk[i][order-1]=0.0;
      Hk[i][order]=0.0;
   }//end for i


// step2: replace terms F(H in F) => K
   trunc_order=order-2;
   for (i=0;i<=order;i++) {
      for (j=0;j<=order;j++) {
	   Kk[i][j]=Fk[i][j];
      }
   }
   for (i=0;i<=order;i++) { //i is corresp. to terms X^i
      // 2*j*Bj=Dj-1 - Dj+1
      for (j=0;j<=order-2;j++) {
          if (fabs(Hk[i][j]) > 1.0E-10) {
             if (j==0) {
		   Kk[i][j+1]=Kk[i][j+1]+2*Hk[i][j];
	     }
	     else {
		for (k=j;k<=trunc_order;k=k+2) {
		   Kk[i][k+1]=Kk[i][k+1]+(2*k)*Hk[i][j];

		}
	    }
          }
      }
      
   }//end for i

// step3: replace terms E(K in E) => M
   trunc_order=order-1;
   for (i=0;i<=order;i++) {
      for (j=0;j<=order;j++) {
	   Mk[i][j]=Ek[i][j];
      }
   }
   for (i=0;i<=order;i++) { //i is corresp. to terms X^i
      // 2*j*Bj=Dj-1 - Dj+1
      for (j=0;j<=order-1;j++) {
          if (fabs(Kk[i][j]) > 1.0E-10) {
             if (j==0) {
		   Mk[i][j+1]=Mk[i][j+1]+2*Kk[i][j];
	     }
	     else {
		for (k=j;k<=trunc_order;k=k+2) {
		   Mk[i][k+1]=Mk[i][k+1]+(2*k)*Kk[i][j];
		}
	    }
          }
      }
   }//end for i


//prepare for solving with lapack
   for (i=0; i<=order; i++) {
      for (j=0; j<=order; j++) {
	A_l1[j*(order+2)+i]=Mk[i][j]; //transpose and pass to an array for Clapack
      }
      A_l1[(order+1)*(order+2)+i]=-Rk[0][i];
      A_l1[(order+2)*(order+2)+i]=-Rk[1][i];
   }
   for (j=0; j<=order; j=j+2) {
       A_l1[j*(order+2)+(order+1)]=1.0;
   }
   for (j=1; j<=order-1; j=j+2) {
       A_l1[j*(order+2)+(order+1)]=-1.0;
   }
// last row of boundary by hand:
  A_l1[0*(order+2)+(order+2)]=0.0;
  A_l1[1*(order+2)+(order+2)]=2.0;
  A_l1[2*(order+2)+(order+2)]=-4.0;
  A_l1[3*(order+2)+(order+2)]=6.0;
  A_l1[4*(order+2)+(order+2)]=-16.0;
  A_l1[5*(order+2)+(order+2)]=20.0;
  A_l1[6*(order+2)+(order+2)]=-36.0;

   for (j=0; j<=order; j++) {
       B_l1[j]=0.0;
   }
   B_l1[order+1]=boundry[0];
   B_l1[order+2]=boundry[1];
//end matrices for Clapack

#define __TEST_COLL_ONLY__
#ifdef __TEST_COLL_ONY__
// display matrices for debug:
  std::cout <<"Content of matrix Rk for debug: \n"; 
  for (i=0;i<=order;i++) {
     std::cout << "Ri= "<< i <<" ";
     for (j=0;j<=1;j++) {
        std::cout << " " << Rk[i][j];
     }
    std::cout << "\n";
  }
#else 
  std::cout <<"Content of matrix Ek for debug: \n" ;
  for (i=0;i<=order;i++) {
     std::cout << "Ei= "<< i <<" ";
     for (j=0;j<=order;j++) {
        std::cout << " " << Ek[i][j];
     }
    std::cout << "\n";
  }
  std::cout <<"Content of matrix Fk for debug: \n" ;
  for (i=0;i<=order;i++) {
     std::cout << "Fi= "<< i <<" ";
     for (j=0;j<=order;j++) {
        std::cout << " " << Fk[i][j];
     }
    std::cout << "\n";
  }
 std::cout <<"Content of matrix Hk for debug: \n" ;
 for (i=0;i<=order;i++) {
     std::cout << "Hi= "<< i <<" ";
     for (j=0;j<=order;j++) {
        std::cout << " " << Hk[i][j];
     }
    std::cout << "\n";
  }
  std::cout <<"Content of matrix Kk for debug: \n" ;
  for (i=0;i<=order;i++) {
     std::cout << "Ki= "<< i <<" ";
     for (j=0;j<=order;j++) {
        std::cout << " " << Kk[i][j];
     }
    std::cout << "\n";
  }
  std::cout <<"Content of matrix Mk for debug: \n" ;
  for (i=0;i<=order;i++) {
     std::cout << "Mi= "<< i <<" ";
     for (j=0;j<=order;j++) {
        std::cout << " " << Mk[i][j];
     }
    std::cout << "\n";
  }

#endif //end test code

   return 0;
}

bool IDENT05_COLL::SolveSeriesLinearSys_ref1() 
{
   int n, rhs, lda, ldb, info;
   int ipiv[9];
 
   n=order+2;
   rhs=1;
   lda=n;
   ldb=n;
   info=1;
 
   // call to ../MathFunctions/mySolveLinearLapack.c
   mySolveLinearLapack(n,rhs,A_l1,lda,ipiv,B_l1,ldb,info);
 
   //should try to display the solution in solarray[2000] points:

}

Number IDENT05_COLL::evalChebyshevPolynom(Number t, Index i)
{
  int j;
  Number sum;
  sum=0.0;
  for (j=0;j<=order;j++) {
     switch(j) {
	case 0:
		sum=sum+t0[j]*pow(t,j);
		break;
	case 1:
		sum=sum+t1[j]*pow(t,j);
		break;
	case 2:
		sum=sum+t2[j]*pow(t,j);
		break;
	case 3:
		sum=sum+t3[j]*pow(t,j);
		break;
	case 4:
		sum=sum+t4[j]*pow(t,j);
		break;
	case 5:
		sum=sum+t5[j]*pow(t,j);
		break;
	case 6:
		sum=sum+t6[j]*pow(t,j);
		break;	
     }
  }
}
