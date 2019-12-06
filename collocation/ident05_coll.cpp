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
#include "../include/ident05_coll.hpp"

#ifndef __TEST_COLL_ONLY__
#define __TEST_COLL_ONLY__
#endif

//using namespace std;
using namespace lhlib;

/** attemps to access optimization parameters via macros */

//#define P1_EPS (ident05_nlp::get_p1_eps)
//#define P2_OM (ident05_nlp::get_p2_om)

// constructor
IDENT05_COLL::IDENT05_COLL(Index ord, const char * FileNameInp, const char * CodeNameInp )
{
/*eqn Problem type: */
//  type_eqn=typode;
  equ1[0]=1.0; equ1[1]=0.0;equ2[2]=1.0;
  equ2[3]=1.0;equ2[4]=0.0;equ2[5]=0.0; equ2[6]=0.0;equ2[7]=1.0;equ2[8]=0.0;  //  (xy"+y'+xy=0) Bessel
 boundry[0]=0.0;
 boundry[1]=1.0;

/* Input File Name */
  lenFileNameInp = 8;  // do not change this value
  strncpy(strFileNameInp, FileNameInp, lenFileNameInp);   // works on MInGW not on Linux
  std::cout << "Have input(1) file name: " << strFileNameInp << "\n";
 // strncpy(strCodeNameInp, CodeNameInp, lenFileNameInp);   // works on MInGW not on Linux
 // std::cout << "Have input(2) file name: " << strCodeNameInp << "\n";

/*  t0=(Number)malloc(7*sizeof(Number)); */
  t0[0]=1.0;t0[1]=0.0;t0[2]=0.0;t0[3]=0.0;t0[4]=0.0;t0[5]=0.0;t0[6]=0.0;
  t1[0]=0.0;t1[1]=1.0;t1[2]=0.0;t1[3]=0.0;t1[4]=0.0;t1[5]=0.0;t1[6]=0.0;
  t2[0]=-1.0;t2[1]=0.0;t2[2]=2.0;t2[3]=0.0;t2[4]=0.0;t2[5]=0.0;t2[6]=0.0;
  t3[0]=0.0;t3[1]=-3.0;t3[2]=0.0;t3[3]=4.0;t3[4]=0.0;t3[5]=0.0;t3[6]=0.0;
  t4[0]=1.0;t4[1]=0.0;t4[2]=-8.0;t4[3]=0.0;t4[4]=8.0;t4[5]=0.0;t4[6]=0.0;
  t5[0]=0.0;t5[1]=5.0;t5[2]=0.0;t5[3]=-20.0;t5[4]=0.0;t5[5]=16.0;t5[6]=0.0;
  t6[0]=-1.0;t6[1]=0.0;t6[2]=18.0;t6[3]=0.0;t6[4]=-48.0;t6[5]=0.0;t6[6]=32.0;

  t7[0]=0.0; t7[1]=-7.0; t7[2]=0.0; t7[3]=56.0; t7[4]=0.0; t7[5]=-112.0; t7[6]=0.0; t7[7]=64.0;
  order=ord; //6
  tinit=0.0;
  tend=2.0;
  num_ranges=2;
  num_points=20;
  num_rows=num_ranges*num_points;
  num_total_coeffs=(order+1)*num_ranges;
  predictparams[0]=1.0;
  predictparams[1]=2.0;
}

IDENT05_COLL::~IDENT05_COLL() 
{
}

bool IDENT05_COLL::read_parse_file() 
{
   Index row,col;

// Open the input file, which is a text file with 6 columns
  std::ifstream lh_file;
  lh_file.open("TheFile", std::ifstream::in);  // test w/o variable Name
  
  row=0;
  while ((!lh_file.eof()) && (row<41)) {
//  while (lh_file.good()) {
      std::string line; 
      std::getline(lh_file, line);    
    
      std::stringstream ss(line);
      col = 0;
      while(ss >> dataarray[row][col]) col++;
#ifdef __TEST_COLL_ONLY__
      std::cout << "tf= "  << dataarray[row][0] << " ";
      std::cout << "utf= " << dataarray[row][1] << " ";
      std::cout << "ytf= " << dataarray[row][2] << " ";
      std::cout << "c2f= " << dataarray[row][3] << " ";
      std::cout << "c1f= " << dataarray[row][4] << " ";
      std::cout << "c0f= " << dataarray[row][5] << "\n";
#endif
      row++;
  }
   num_rows=row;
   lh_file.close();

   return 0;
}


bool IDENT05_COLL::read_parse_code() {
   int maxrow, row, col;
   char tuple_t0[8]="REL";
   char tuple_t1[12]="INT:INIPAR";
//   const char* tuple_t;
   Number doublearray[1][2];
   Index intarray[1][2];
/* *
 * Open "TheCode" which is a parameter file with 6 rows of doublets:
 * 1: int:bvp if (1)or ivp, type of diff eqn
 * 2: prediction(1) or given boundry[1] if (0), repeat boundry or use once
 * 4: double:tinit,tend
 * 3: double:boundry[0], boundry[1] as bvp or ivp
 * 5: int: num_ranges for calc, num_points for output
 * 6: double:prediction[0] and prediction[1] (amplitude and period of sinusoide)
 * */
   std::ifstream lh_code;
   lh_code.open("TheCode", std::ifstream::in);
   std::string line;
   std::getline(lh_code, line);
   std::stringstream ss(line);
   ss >> maxrow;
   std::cout << "have " << maxrow << " Code Parameters\n";
   row=1;
   std::div_t dv{}; dv.quot = row;
   while ( (std::getline(lh_code, line)) && (row<(4*maxrow)) ) {
      std::stringstream ss(line);
      col=0;
      dv.quot=row;
      dv = std::div(dv.quot, 2);
      std::vector<std::array<char, 8>> listchar;   //stackoverflow

      if (dv.rem) {     // odd rows: typ and name of params
//        while (ss >> tuple_t[0][col]) col++;  // do not work for const char*
          //std::istringstream input("INT:INIEQN:END");

          // note: the following loop terminates when std::ios_base::operator bool()
          // on the stream returned from getline() returns false

         // for (std::array<char, 8> ar; ss.getline(&ar[0], 8, ':'); ) {//stackoverflow$
         //       listchar.push_back(ar);  // name
         // }

	  int count=1;
           for (count=1;count<3;count++) {
	       const std::string tmp = ss.str();
          //for (auto& ar :listchar) {
           //    std::cout << &ar[0] << '\n';   // Test the format (type,name)
               if (count==1) {
		 strncpy( tuple_t0, (const char*) tmp.c_str(), 3);
		 std::cout << "type= " << tuple_t0 << "\n";       
	       }

	       else {
	//	 tuple_t=(const char*) tmp.c_str();
		 strncpy( tuple_t1, (const char*) tmp.c_str(), 11);
                 std::cout << "name= " << tuple_t1 << "\n"; 
	       }
       	  }
	  
      }//end odd rows
      else { 
	      // odd rows: params values
       if (!(strcmp(tuple_t0, "INT"))) {    // if not == strcmp  MATCH INT
         while (ss >> intarray[0][col]) col++;
         if (!strcmp(tuple_t1, "INT:IBVEQN:")) {
           type_ovp = intarray[0][0];
	   type_eqn = intarray[0][1];
	 }	 
         if (!strcmp(tuple_t1, "INT:PRELOG:")) {
           type_predict = intarray[0][0];
	   repeat_predict = intarray[0][1];
	 }	
         if (!strcmp(tuple_t1, "INT:NUMPTS:")) {
           num_ranges = intarray[0][0];
	   num_points = intarray[0][1];
	 }
         std::cout << "valuei1= " << intarray[0][0] << " and valuei2= " << intarray[0][1] << "\n";	 
       }
       else if (!(strcmp(tuple_t0, "DOU"))) {     // if not == strcmp MATCH DOU
         while (ss >> doublearray[0][col]) col++;
         if (!strcmp(tuple_t1, "DOU:INIEND:")) {
	    tinit = doublearray[0][0];
            tend  = doublearray[0][1];
	 }
         if (!strcmp(tuple_t1, "DOU:BVALUE:")) {
		 //both initial and boundry are filled with the same values
		 // but only one of them shall be used whether type_ovp=1 (boundry)
		 // or 0 (initial)
            initial[0] = doublearray[0][0];
	    initial[1] = doublearray[0][1];
            boundry[0] = doublearray[0][0];
            boundry[1]  = doublearray[0][1];

	 }
         if (!strcmp(tuple_t1, "DOU:PREVAL:")) {
	    predictparams[0] = doublearray[0][0];
            predictparams[1] = doublearray[0][1];
	 }
         std::cout << "valued1= " << doublearray[0][0] << " and valued2= " << doublearray[0][1] << "\n"; 	 
       }//end if INT:DOU

      }//end even rows, (if dv.rem)
      row=row+1;
   }//end while rows

#ifdef __TEST_COLL_ONLY__
      std::cout << "number of code params= " << maxrow << "\n";
      std::cout << "type ovp problem= " << type_ovp << " , order problem= " << order << "\n";      
      std::cout << "use prediction function= " << type_predict << " , repeat prediction all ranges= " << repeat_predict << "\n";       
      std::cout << "initial time= " << tinit << " , tend= " << tend << "\n";
      std::cout << "numranges= " << num_ranges << " , num points= " << num_points << "\n";
      std::cout << "sinus prediction amplitude= " << predictparams[0] << "period= " << predictparams[1] << "\n";

#endif
   lh_code.close();
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
  std::cout << "Expanding equation: " << equ1[2] << "*Y'' +" << equ1[1] << "*Y' +" << equ1[0] << "*Y = 0 \n";
// step1: create system: H*Di+F*Bi+E*Ai=R
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

   //terms of matrix Fk corresponds to first derivative term
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

   //terms of matrix Hk corresponds to second derivative term
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

// step3: replace terms E(K in E) => M now have the system: H*Di+F*Bi+E*Ai=0 into M*Ai=R

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
	A_l1[j*(order+3)+i]=Mk[i][j]; //transpose and pass to an array for Clapack
      }
      A_l1[(order+1)*(order+3)+i]=-Rk[i][0];
      A_l1[(order+2)*(order+3)+i]=-Rk[i][1];
   }
/***********
   for (j=0; j<=order; j=j+2) {
       A_l1[j*(order+3)+(order+1)]=1.0;
   }
   for (j=1; j<=order-1; j=j+2) {
       A_l1[j*(order+3)+(order+1)]=-1.0;
   }
*********/
// last row of boundary by hand:
  A_l1[0*(order+3)+(order+1)]=1;
  A_l1[1*(order+3)+(order+1)]=-1;
  A_l1[2*(order+3)+(order+1)]=1;
  A_l1[3*(order+3)+(order+1)]=-1;
  A_l1[4*(order+3)+(order+1)]=1;
  A_l1[5*(order+3)+(order+1)]=-1;
  A_l1[6*(order+3)+(order+1)]=1;


// last row of boundary by hand:
  if (type_ovp==0) {
  A_l1[0*(order+3)+(order+2)]=0.0;
  A_l1[1*(order+3)+(order+2)]=1.0;//2.0//2.0
  A_l1[2*(order+3)+(order+2)]=-4.0;//-4.0//0.0
  A_l1[3*(order+3)+(order+2)]=9.0;//6.0//0.0
  A_l1[4*(order+3)+(order+2)]=-16.0;//-16.0//0.0
  A_l1[5*(order+3)+(order+2)]=25.0;//20.0//10.0
  A_l1[6*(order+3)+(order+2)]=-36.0;//-36.0//0.0
  for (j=0; j<=order; j++) {
       B_l1[j]=0.0;
   }
   B_l1[order+1]=initial[0];
   B_l1[order+2]=initial[1];  
  }
  else 
  {
  A_l1[0*(order+3)+(order+2)]=1.0;
  A_l1[1*(order+3)+(order+2)]=1.0;
  A_l1[2*(order+3)+(order+2)]=1.0;
  A_l1[3*(order+3)+(order+2)]=1.0;
  A_l1[4*(order+3)+(order+2)]=1.0;
  A_l1[5*(order+3)+(order+2)]=1.0;
  A_l1[6*(order+3)+(order+2)]=1.0;	  
  for (j=0; j<=order; j++) {
       B_l1[j]=0.0;
   }
   B_l1[order+1]=boundry[0];
   B_l1[order+2]=boundry[1];
  } 
//check the four last col-row are equal to zero:
  A_l1[7*(order+3)+(order+1)]=0.0; 
  A_l1[7*(order+3)+(order+2)]=0.0; 
  A_l1[8*(order+3)+(order+1)]=0.0; 
  A_l1[8*(order+3)+(order+2)]=0.0; 


//end matrices for Clapack

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

  std::cout <<"Content of matrix A_l1 for debug: \n" ;
  for (j=0;j<=order+2;j++) {
     std::cout << "Aj= "<< j <<" ";
     for (i=0;i<=order+2;i++) {
        std::cout << " " << A_l1[j+i*(order+3)];
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
 
   n=order+3;
   rhs=1;
   lda=n;
   ldb=n;
   info=1;
 
   // call to ../MathFunctions/mySolveLinearLapack.c
   mySolveLinearLapack(n,rhs,A_l1,lda,ipiv,B_l1,ldb,info);
return 0;
}

bool IDENT05_COLL::NumRangesCalcBoundary(Number* boundary_all, Number* times_end, Index k) {
  Number t;
        
	/* REPEAT BOUNDARY FOR NEXT RANGE: this function should be called at init k=0 and after SolveSeriesLinearSys based on boundary_all[2*k,2*k+1] */
	 /* Hence this function should be called with arg3=k+1 */

   if (type_ovp == COL_TYP_INITIAL) {  // Initial Value Problem
   // in the case of Initial Value Problem, boundry[0] is the state and boundry[1] is the derivative
      if (k==0) {
         boundary_all[0]=initial[0];
         boundary_all[1]=initial[1];
      } 
      else { // temptative of jointure of piecewise polynoms but will this not diverge? 
	 boundary_all[2*(k)]=solarray[num_points*(k-1)+num_points-1][1];
	 boundary_all[2*(k)+1]=( solarray[num_points*(k-1)+num_points-1][1] 
                         - solarray[num_points*(k-1)+num_points-2][1] ) 
                        / ( solarray[num_points*(k-1)+num_points-1][0] 
                         - solarray[num_points*(k-1)+num_points-2][0] );
      }
       // for an initial value problem the second boundary condition is an initial condition and shall be scaled to the interval [-1,1]
      boundary_all[2*(k)+1] = boundary_all[2*(k)+1] / (2/(times_end[k+1]-times_end[k]));
   } // if type ovp:
   else {  // Boundary Value Problem
     if (type_predict==COL_TYP_BVP_PRED_TRICKED) {
       if (k==0) {
         boundary_all[0]=boundry[0];
         boundary_all[1]=boundry[1];
       } 
       else {
         t=times_end[k];
         boundary_all[2*(k)]=predictparams[0]*sin(2*3.14156*t/predictparams[1]); // does not work if the period is not captured
         boundary_all[2*(k)+1]=predictparams[0]*2*3.14156/predictparams[1]*cos(2*3.14156*t/predictparams[1]); // defaut but does not work if the period is not captured
       }
     }
     else { //type_predict=1 (COL_TYP_BVP_PRED_NORM , use a function)
        if (repeat_predict==COL_TYP_BVP_NREPEAT) {  // use the predictor function as many times as ranges
          t=times_end[k];
	  boundary_all[2*(k)]=predictparams[0]*sin(2*3.14156*t/predictparams[1]);

	  t=times_end[k+1];
	  boundary_all[2*(k)+1]=predictparams[0]*sin(2*3.14156*t/predictparams[1]);

         } 
         else {  // use the predictor function the same manner as the first range, although this is not the philosophy of collocation
         boundary_all[2*(k)] = 0.0;
         boundary_all[2*(k)+1] = predictparams[0]*2*3.14156; // /predictparams[1]
         }
     } // endif type_predict
  }//endif type_ovp
  return 0;
}


bool IDENT05_COLL::SolveNumRangesSys_ref1()
{
   Index k, j;
   Number t, x, times_end[5];
   Number boundary_all[10]; 
// 10 should be set to a variable index in the future as it corresponds to a number of points in a simulation intervall (10) and is time-sampling dependant
   Number equ_all[12];

//   equ_all[0]=equ1[0];
//   equ_all[1]=equ1[1];
//   equ_all[2]=equ1[2];
   for (k=0; k<num_ranges; k++) {
       // iterate on intervals
       times_end[k]=tinit + k*(tend-tinit)/num_ranges;
   }   
   times_end[num_ranges]=tend;
   k=0;
   if (NumRangesCalcBoundary(boundary_all, times_end, k)) {
      std::cout << "Error Calc Init Boundary\n"; 
   }
   for (k=0; k<num_ranges; k++) {
	// import equation coefficients, which are set in the Data TheFile.txt
       equ_all[3*k+2]=dataarray[k*num_points+0][3];
       equ_all[3*k+1]=dataarray[k*num_points+1][4];
       equ_all[3*k]=dataarray[k*num_points+2][5];
	// correction of equation coefficient, due to reduction of intervall size to -1..1:
       equ_all[3*k]=equ_all[3*k]*(1); 
       equ_all[3*k+1]=equ_all[3*k+1]*(2/(times_end[k+1]-times_end[k]));
       equ_all[3*k+2]=equ_all[3*k+2]*pow(  (2/(times_end[k+1]-times_end[k])), 2);
#ifdef __TEST_COLL_ONLY__       
	std::cout << "    check boundary conditions: y0 " <<  boundary_all[2*k] << " and dy0 " <<  boundary_all[2*k+1] << " \n";
#endif
       boundry[0]=boundary_all[2*k];
       boundry[1]=boundary_all[2*k+1];

       equ1[0]=equ_all[3*k];
       equ1[1]=equ_all[3*k+1];
       equ1[2]=equ_all[3*k+2];
      
        // call to form the system
       ExpandSeriesLinearSys_ref1();
	// now the system is set in matrices A_l1, B_l1
       SolveSeriesLinearSys_ref1();
	// now the system is solved in matrices B_l1
       for (j=0;j<(order+1);j++) {   // 7 is equal to order plus 1
           //copy the result from B_l1
           coeffarray[(order+1)*k+j]=B_l1[j];
       }

#ifdef __TEST_COLL_ONLY__
        std::cout << " ****** Solution of Linear System: " << k << " c0=" << equ1[0] << " c1=" << equ1[1] << " c2=" << equ1[2] << ": ******** \n";
	std::cout << "    with boundary conditions: y0 " <<  boundry[0] << " and dy0 " <<  boundry[1] << " \n";
        for (j=0;j<(order+1);j++) {   // 7 is equal to order plus 1
          // display the result
           std::cout << coeffarray[(order+1)*k+j] << '\n';
       }
#endif           

	//evaluate the result of the solution in solarray vector
       for (j=0; j<=num_points; j++) {
           x=-1+2.0/float(num_points)*j;
           t=(times_end[k+1]-times_end[k])/2.0*x + (times_end[k+1]+times_end[k])/2.0;  // is there an error?
	   solarray[num_points*k+j][0]=t;
	   solarray[num_points*k+j][1]=evalCollocation(x, B_l1); // x, not t
       }

      // anticipate next iteration and calc correspondant boundaries            
      if (NumRangesCalcBoundary(boundary_all, times_end, k+1)) {
          std::cout << "Error iteration of Calc Init Boundary Condition\n";
      }

   } //end for k
   
   // SAVE ON DISK the solution: coefficients and time values
   // SHALL BE REPLACED BY ANOTHER DEDICATED CLASS METHOD
  int row;
  std::ofstream lh_spec;
  lh_spec.open("TheSpec", std::ofstream::out);  // test w/o variable Name
  row=0;
  lh_spec << num_ranges << "\n";
  while (k<num_ranges) {
        for (j=0;j<(order+1);j++) {   // 7 is equal to order plus 1
          // display the result
           lh_spec << coeffarray[(order+1)*k+j] << "\n";
         }
   }
  lh_spec.close();

  std::ofstream lh_test;
  lh_test.open("TheTest", std::ofstream::out);  // test w/o variable Name
  row=0;
  while (row<40) {
//  while (lh_test.good()) {
      lh_test << solarray[row][0] << "\t";
      lh_test << solarray[row][1] << "\n";
      row++;
  }
  lh_test.close();

#ifdef __TEST_COLL_ONLY__  
  std::cout << " *************** This is the solution solved by collocation ****************** \n"; 
  row=0;
  while (row<40) {
      std::cout << "tf= " << solarray[row][0] << " ";
      std::cout << "ytf= " << solarray[row][1] << "\n";
      row++;
  }
#endif

  return 0;
}

Number IDENT05_COLL::evalCollocation(Number t, Number * coeff) 
{
  int j;
  Number sum;
  sum=0.0;
  for (j=0;j<=order;j++) {
//#ifdef __TEST_COLL_ONLY__ 
//	std::cout << "coeff[" << j << "]= " << coeff[j] << "\n";
//#endif	  
     sum=sum+coeff[j]*evalChebyshevPolynom(t,j);
  }
  return sum;
}

Number IDENT05_COLL::evalChebyshevPolynom(Number t, Index i)
{
  int j;
  Number sum;
  sum=0.0;
  for (j=0;j<=order;j++) {
     switch(i) {
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
  return sum;
}

bool IDENT05_COLL::pass_dataarray_col(Index dim, li_doubles &exports) {

   Index row;
   Number data[2];

   for (row=0; row< dim; row++) {
 // REVERSE OF exportClassInst.arrayToExport[row][sig]=ptr_sig[row];
     data[0]=solarray[row][0];
     data[1]=solarray[row][1];
     exports.push_back(data[0]);
     exports.push_back(data[1]);     
   }
  return 1;
}

