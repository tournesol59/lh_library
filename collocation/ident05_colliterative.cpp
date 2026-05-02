/*
 * Solve a second-order linear (for the moment) equation
 * coefficients could be a polynom up to degree two.
 * The collocation method (with Chebyshev interpolation)
 * This part deals with the analysis of the previous step solution
 * to update the boundary values
 *
 * This file contains the procedures:
 * UpdateBoundaryIterativeSlope()
 * UpdateBoundaryIterativePhase()
 *
 * UpdateBoundaryIterative() as a whole is still present but deprecated
 *
 * LIEBHERR TOULOUSE
 *******************************************************/

#include "../include/ident05_coll.hpp"

#ifndef __TEST_COLL_ONLY__
#define __TEST_COLL_ONLY__
#endif

#ifndef __TEST_COLL_SIMPLEDPHI__
#define __TEST_COLL_SIMPLEDPHI__
#endif

bool IDENT05_COLL::UpdateBoundaryIterativeSlope()
{
  Number LAMBDA1 = 0.7/predictparams[1];
  Number LAMBDA2 = 0.6/predictparams[1];
  Number COLL_MIN_DERIVATIVE = predictparams[0]*2*3.14156/predictparams[1]*0.16;
  Number COLL_MAX_MISMATCH = 1.0*COLL_MIN_DERIVATIVE;

  Index k; // subinterval no and various indexes: for search (h) and location of extremum of y (h0) and of y' (h1)
  //Number kcoeffs[8];
  Number x_a, x_b, y_b, dy_b, dyb_next, pred_dy_a, pred_dy_b, 
	 y_a, 
	 diff_dy_b, new_y_a, new_y_b; // 
  std::cout << "DEBUG, predictparams " << predictparams[0] << " " << predictparams[1] << " " << predictparams[2] << " " << predictparams[3] << "\n";

  // LOGIC is present here!
 Miter++; // increase of global index of iteration
 int enable_new_dyabs=0;
 if ( enable_new_dyabs == 1 ) {

  for (k=0; k < num_ranges-1; k++) {  // we only evaluate k-1
     // normally import TheCode[k] here, but for a two points config we can test w/o it
     y_a = boundary_all[2*k];
     y_b = boundary_all[2*k+1];

     //dy_a = boundary_derivall[2*k];
     dy_b = boundary_derivall[2*k+1];  //REMARK: this convention for interv k
     dyb_next = boundary_derivall[2*k+2];

     x_a = times_end[k];
     x_b = times_end[k+1];
//     x_c = times_end[k+1];    
//     we neglect the exponential decay effect for the moment 
     pred_dy_a = 2*3.14156/predictparams[1]*predictparams[0]*sin(2*3.14156*x_a/predictparams[1]+predictparams[3]);
     pred_dy_b = 2*3.14156/predictparams[1]*predictparams[0]*sin(2*3.14156*x_b/predictparams[1]+predictparams[3]);
     //std::cout << "DEBUG predicted param3 " << predictparams[3] << " ypa " << pred_dy_a << " ypb " << pred_dy_b << "\n";
     //diff_dy_a = 0; // unused
     diff_dy_b = - dy_b + dyb_next; // same REMARK
     if (k > 0) {
        new_y_a = boundary_all[2*(k-1)+1];  // corresp to new_y_b[k-1], and do nothing if if k==0
     } 
     else {
        new_y_a = boundary_all[0]; // ok to recopy inverted after
     }
     if (fabs(pred_dy_b) > COLL_MIN_DERIVATIVE ) {
	     // fred: here add a modification by authorizing the modification
	     // of y_c ! , always based on slopes difference at x_b and only for
	     // a derivative threshold of dy_b !
	     // fred (15 oct 25): I have changed the sign from minus to plus
        new_y_b = (fabs(diff_dy_b) > COLL_MAX_MISMATCH)? y_b + LAMBDA1*diff_dy_b : y_b;
	//new_y_c = (fabs(diff_dy_b) > COLL_MAX_MISMATCH)? y_c + LAMBDA2*diff_dy_b : y_c;
     }
     else {
       // preview another algorithm
       new_y_b = boundary_all[2*k+1];
     }
     if (Miter > 0) {
	std::cout << "Boundary update, k=" << k << " x_a=" << x_a << ", x_b=" << x_b << "\n";
	std::cout << "Read calc value: ya " << y_a << " yb " << y_b << "\n";
        std::cout << "Predict derivative: dypa " << pred_dy_a << ", dypb " << pred_dy_b << ", min threshd: " << COLL_MIN_DERIVATIVE << "\n";
	std::cout << "Read derivative: dyb " << dy_b << ", dybnext " << boundary_derivall[2*(k+1)] << "\n";
        std::cout << "Deviation derivative: diffdyb " << diff_dy_b << ", threshd: " << COLL_MAX_MISMATCH << ", factors: lambdab=" << LAMBDA1 << ", lambdac=" << LAMBDA2 << "\n";
        std::cout << "new bound value: nya " << new_y_a << " nyb " << new_y_b << "\n";

     }     
     // write the result in array, but this must be written also on file
     boundary_all[2*k]=new_y_a;
     boundary_all[2*k+1]=new_y_b;
  }
 }
  return 0;
}

// subroutine used by UpdateBoundaryIterativePhase()
bool IDENT05_COLL::exploreForExtremumFirstOdd(int k, int &found1, int &found2, double &yabs, double &dyabs, double &dphi)
{
   Number COLL_MIN_DERIVATIVE = predictparams[0]*2*3.14156/predictparams[1]*0.35;
  Index h, h_s, h0, h1, hend; // subinterval no and various indexes: for search (h) and location of extremum of y (h0) and of y' (h1)
  Number yh, dyh, yh_s, dyh_s, yh1, dyh1, yhend, yM, dyM; 
  //Number y_a, y_b; 
  Number dphi1; //dphi0;

  // TBD insert 2 while loop
        h=k*num_points;
	h_s = h;
	yh_s=solarray[h_s][1];
	dyh_s=solarray[h_s][2];
	yM=solarray[h][1]; // value
	yh=solarray[h+1][1];
	while ((h < (k+1)*num_points) && (fabs(yM) < fabs(yh))) { // look up for a maximum of y
          h++;
	  yM=yh;
          yh=solarray[h+1][1];
	}  
	h0=h;
        // reinit
	h=k*num_points;
	dyM=solarray[h][2]; // derivative
	dyh=solarray[h+1][2];
	while ((h < (k+1)*num_points) && ((fabs(dyM) < fabs(dyh)) || (fabs(dyM) < 2.0*COLL_MIN_DERIVATIVE) )) { // look up for a minimum of y'
          h++;
	  dyM=dyh;
	  dyh=solarray[h+1][2];
	  yh=solarray[h+1][1];
	}
	h1=h;
	if (h1==h0) { // security against NaN values
	     h1=h0+1;
       	}
	yh1 = yh;
	dyh1 = dyh;
	hend=(k+1)*num_points-1;
	yhend=solarray[hend][1];

	if (((yh1 > 0) && (yM > 0)) || ((yh1 < 0) && (yM < 0)))  { // eg. phi=k*Pi...(k+1/2)*Pi-eps
            if ((h0-h_s) <=1) {
	// formula, correct version is line one below: (with scale of phase)
               dphi1 = atan( fabs(dyh1)/fabs(yh1)/(6.283/predictparams[1])); // checked
	       if (h1 < hend) {
                   dphi1 = dphi1 + atan(fabs(yhend-yh1)*6.283/predictparams[1]/fabs(dyh1)); //TBC
	       }	       	
	       std::cout << "DEBUG JUST AFTER phi" << k+1 << " " << dphi1 << " calc: yh1=" << yh1 << " dyh1=" << dyh1  << " t1=" << solarray[h0][0] << " t2=" << solarray[h1][0] << "\n";
            }
	    else if ((h0-h_s) > 1) {
	       dphi1 = 3.1415 - atan( fabs(yh_s)*6.283/predictparams[1]/fabs(dyh_s));
	       if (h1 < hend) {
                   dphi1 = dphi1 + atan(fabs(yhend-yh1)*6.283/predictparams[1]/fabs(dyh1)); //TBC
	       }
	       	std::cout << "DEBUG JUST AFTER phi" << k+1 << " " << dphi1 << " calc: yh_s=" << yh_s << " dyh_s=" << dyh_s  << " t1=" << solarray[h0][0] << " t2=" << solarray[h1][0] << "\n";
	    }
	}
	else if (((yh1 < 0) && (yM > 0)) || ((yh1 > 0) && (yM < 0))) { // eg. phi=k*Pi...(k+1/2)*Pi+eps
            if ((h0-h_s) <=1) {
               dphi1 = 1.5707+atan( fabs(yh1)*(6.283/predictparams[1])/fabs(dyh1)); //TBC
	       if (h1 < hend) {
                   dphi1 = dphi1 + atan(fabs(yhend-yh1)*6.283/predictparams[1]/fabs(dyh1)); //TBC
	       }		
	       std::cout << "DEBUG JUST AFTER phi" << k+1 << " " << dphi1 << " calc: yh1=" << yh1 << " dyh1=" << dyh1  << " t1=" << solarray[h0][0] << " t2=" << solarray[h1][0] << "\n";
	    }
	    else if ((h0-h_s) > 1) {
	       dphi1 = 3.1415 - atan( fabs(yh_s)*6.283/predictparams[1]/fabs(dyh_s));
	       if (h1 < hend) {
                   dphi1 = dphi1 + atan(fabs(yhend-yh1)*6.283/predictparams[1]/fabs(dyh1)); //TBC
	       }	
	       std::cout << "DEBUG JUST AFTER phi" << k+1 << " " << dphi1 << " calc: yh_s=" << yh_s << " dyh_s=" << dyh_s  << " t1=" << solarray[h0][0] << " t2=" << solarray[h1][0] << "\n";

	    }
	}

  // result
  found1 = h0;
  found2 = h1;
  yabs = yM;
  dyabs = dyM;
  dphi = dphi1;
  return 0;
}

/*** second subroutine called by UpdateBoundaryIterativePhase ****/
bool IDENT05_COLL::exploreForExtremumSecEven(int k, int &found1, int &found2, double &yabs, double &dyabs, double &dphi)
{   
  Number COLL_MIN_DERIVATIVE = predictparams[0]*2*3.14156/predictparams[1]*0.35;
  Index h, h1, h2; // subinterval no and various indexes: for search (h) and location of extremum of y (h0) and of y' (h1)
  Number yh, dyh, yh1, dyh1, yh2, dyh2, yM, dyM; 
  //Number y_a, y_b; 
  Number dphi2; //dphi1

  // insert while loop
        h=num_points*k;
	dyM=solarray[h][2];
	dyh=solarray[h+1][2];
        while ((h < (k+1)*num_points) && ((fabs(dyM) < fabs(dyh)) || (fabs(dyh) < 2.0*COLL_MIN_DERIVATIVE)) ) { // look up y'(max)
           h++;
	   dyM=dyh;
	   dyh=solarray[h+1][2];
	}
	h1=h;
	dyh1=solarray[h1][2];
	yh1=solarray[h1][1];
	//extremum_alldphi[2*k+1]=asin(fabs(dyM)/yM/(6.283/4./(solarray[h1][0]-solarray[h0][0])) ); // fred: leave here as comment since another formula is used below
        h=num_points*k;
	yM=solarray[h][1];
	yh=solarray[h+1][1];
	while ((h< (k+1)*num_points) && ( ((fabs(yM) < fabs(yh)) && (yh<0)) || ((fabs(yM) > fabs(yh)) && (yh>0)) )) { // look up for a minimum of yi, FRED: OK!
          h++;
	  yM=yh;
          yh=solarray[h+1][1];
	}
	h2=h;
	dyh2=solarray[h2][2];
	yh2=solarray[h2][1];
	
	if ((dyh2 > 0) && (yh2 < 0)) {  // corrected here: sign(dyh2) test
	   // FRED: les formule ci-dessouss
	    dphi2=(3.1415-atan(fabs(yh2)*(6.28/predictparams[1])/fabs(dyh2)));
	   // FRED : essayer en prenant en compte les abcisses pour 'mesurer' la dilatation 
	 // dphi2 = (3.1415-asin( fabs(yh2)*(6.28/predictparams[1])/( fabs(dyM)) ));
	  //extremum_alldphi[2*k+1]=dphi2*(predictparams[1]/4./(solarray[h2][0]-solarray[h1][0]));
	   // FRED : this might be OK but prefer the following 
	  //extremum_alldphi[2*k+1]=3.1415-asin( fabs(yh2)*(6.28/predictparams[1])/( fabs(dyM)) );
	} 
	else if ((dyh2 < 0) && (yh2 < 0)) {
	   // FRED: la formule ci-dessous qui est la meilleure est la premiere TBD
	   // ERREUR dphi2 = atan(fabs(dyh2)/fabs(yh2)/(6.28/predictparams[1]));
	   dphi2 = atan(fabs(yh2)*(6.28/predictparams[1])/fabs(dyh2));
           //extremum_alldphi[2*k+1]=atan(fabs(dyh2)/fabs(yh2)/(6.28/4./(solarray[h2][0]-solarray[h1][0])));
	  // RETURN Nan dphi2 = asin( fabs(yh2)*(6.28/predictparams[1])/( fabs(dyM)) );
	  //extremum_alldphi[2*k+1]=dphi2*(predictparams[1]/4./(solarray[h2][0]-solarray[h1][0]));
	   //extremum_alldphi[2*k+1]=asin( 0.9*fabs(yh2)*(6.28/predictparams[1])/(fabs(dyM)) );

	}
        else if ((dyh2 > 0) && (yh2 >0)) {
	    //std::cout << "logique PAS ENCORE Implementee pour PHI" << "\n";
	    dphi2=(3.1415-atan(fabs(yh2)*(6.28/predictparams[1])/fabs(dyh2)));
	}
        else {
	    //std::cout << "logique PAS ENCORE Implementee pour PHI" << "\n";
	    dphi2 = atan(fabs(yh2)*(6.28/predictparams[1])/fabs(dyh2));
	}	//endif 1
	   std::cout << "DEBUG JUST AFTER phi" << k+1 << " " << dphi2 << " calc: yh2=" << yh2 << " dyh2=" << dyh2  << " t1=" << solarray[h1][0] << " t2=" << solarray[h2][0] << "\n";

  // result
  found1 = h1;
  found2 = h2;
  yabs = yM;
  dyabs = dyM;
  dphi = dphi2;
  return 0;
}


bool IDENT05_COLL::UpdateBoundaryIterativePhase() {
  Number LAMBDA_PHI = 0.50/1.57;
  Number COLL_MIN_DPHI = 0.1;
  Number COLL_DPHI_EPS = 0.025;
  Index found1, found2;
  Index k; // subinterval no and various indexes: for search (h) and location of extremum of y (h0) and of y' (h1)
  Number yabs, dyabs; 
  Number y_a, y_b, new_y_a, new_y_b, diff_dphi; 
  Number dphi1, dphi2, dphi3;

  // now look for the phase increase
// first thing is to populate the extremum_alllocus, _allvalues tables
  for (k=0; k < num_ranges; k++) {
	  
//     if (((k%4 == 0) || (k%4 == 1)) && (boundary_all[2*k] > 0.5*predictparams[0])) { // part decreasing
     if (k%4 == 0) {       
	// cosinus-like part, decreasing
	
        exploreForExtremumFirstOdd(k, found1, found2, yabs, dyabs, dphi1);
	
	// store result
	if ((found1>=0) && (found2>=0)) {
	   extremum_alllocus[2*k]=found1;
	   extremum_allval[2*k]=yabs;
	   extremum_alldphi[2*k]=0.;
	   extremum_alllocus[2*k+1]=found2;
	   extremum_allval[2*k+1]=dyabs;
#ifdef __TEST_COLL_SIMPLEDPHI__
	   extremum_alldphi[2*k+1]=dphi1;
           std::cout << "DEBUG Phase increase, ph"<< 2*k << "=" << extremum_alldphi[2*k] << ", ph" << 2*k+1 << "=" << extremum_alldphi[2*k+1] << "\n" ;

#else
	   extremum_alldphi[2*k+1]=dphi1 * predictparams[1]/4./(solarray[found2][0]-solarray[found1][0]);
#endif
	}
	else {
           std::cout << "No extremum found in exploreForExtremum with k=" << k << "\n";
	}

     }
//    else if (((k%4 == 0) || (k%4 ==1)) && (fabs(boundary_all[2*k]) < 0.2*fabs(predictparams[0]) )) { // sinus increasing or cosinus decreasing part
     else if (k%4 == 1) {  // sinus-like part
        exploreForExtremumSecEven(k, found1, found2, yabs, dyabs, dphi2);
        if ((found1>=0) && (found2>=0)) {
          extremum_alllocus[2*k]=found1;
	  extremum_allval[2*k]=dyabs;
	  extremum_alldphi[2*k]=0.;

	  extremum_alllocus[2*k+1]=found2;
	  extremum_allval[2*k+1]=yabs;
#ifdef __TEST_COLL_SIMPLEDPHI__
	  extremum_alldphi[2*k+1]=dphi2;
#else
	  extremum_alldphi[2*k+1]=dphi2 * predictparams[1]/4./(solarray[found2][0]-solarray[found1][0]);
#endif
	}
	else {
		std::cout << "No extremum found in exploreForExtremum with k=" << k << "\n";
	}
     } // endif 2   
   //    if (((k%4 ==2) && (boundary_all[2*k] < 0.5*predictparams[0])) || ((k%4==0) && (boundary_all[2*k] > 0.5*predictparams[0])) ) { // cosinus increase or sinus increasing
     if (k%4 == 2) {  // cosinus-like part, increasing
	// cosinus-like part, decreasing
	
        exploreForExtremumFirstOdd(k, found1, found2, yabs, dyabs, dphi3);
	
	// store result
	if ((found1>=0) && (found2>=0)) {
	   extremum_alllocus[2*k]=found1;
	   extremum_allval[2*k]=yabs;
	   extremum_alldphi[2*k]=0.;
	   extremum_alllocus[2*k+1]=found2;
	   extremum_allval[2*k+1]=dyabs;
#ifdef __TEST_COLL_SIMPLEDPHI__
	   extremum_alldphi[2*k+1]=dphi3;
#else
 	   extremum_alldphi[2*k+1]=dphi3 * predictparams[1]/4./(solarray[found2][0]-solarray[found1][0]);
#endif
	}
	else {
           std::cout << "No extremum found in exploreForExtremum with k=" << k << "\n";
	}
     } // endif1  
  }// end for k

  // now re iterate loop over k intervals, this time knowing the max locus
  
  for (k=0; k < num_ranges-1; k++) {
     std::cout << "Phase increase, ph"<< 2*k << "=" << extremum_alldphi[2*k] << ", ph" << 2*k+1 << "=" << extremum_alldphi[2*k+1] << "\n" ;
  }

  // now the update
  int enable_new_dphi=1;
  int ident_case_dphi=0;
  if ( enable_new_dphi == 1 ) {
	// loop must start at k=1 !
   for (k=0; k < num_ranges; k++) {  // for (k=1: < num_ranges-1; k++) {
    if (k%4 == 0) {
     y_a=boundary_all[2*k];
     y_b=boundary_all[2*k+1];
     std::cout << "recheck boundary condition before update for k=" << k << ", y_a=" << y_a << " y_b=" << y_b << "\n";
       if (extremum_alldphi[2*k+1] > 1.5707+COLL_MIN_DPHI) {
           //strategy1: give more amplitude at right-end /
	   ident_case_dphi = 0;
	   diff_dphi = (extremum_alldphi[2*k+1] - 1.5707);
	   if (y_a > 0) {
	      new_y_a = y_a;
	      new_y_b = y_b + LAMBDA_PHI*diff_dphi;
	   } else {
	      new_y_a = y_a;
	      new_y_b = y_b - LAMBDA_PHI*diff_dphi;
	   }
       }
       else { // dphi < 1.5707
	   ident_case_dphi = 1;
	   diff_dphi = (extremum_alldphi[2*k+1] - 1.5707);
	   if (y_a > 0) {
	      new_y_a = y_a;
	      new_y_b = y_b + LAMBDA_PHI*diff_dphi;
	   } else {
	      new_y_a = y_a;
	      new_y_b = y_b - LAMBDA_PHI*diff_dphi;
	   }
      
       }
    }
    else if ((k%4 ==1 ) || (k%4 == 2)) {
     y_a=boundary_all[2*k];
     y_b=boundary_all[2*k+1];
     std::cout << "recheck boundary condition before update for k=" << k << ", y_a=" << y_a << " y_b=" << y_b << "\n";
      // OLD VERSION:
      // compare phase increase from previous interval
      // diff_dphi = (extremum_alldphi[2*k+1]-extremum_alldphi[2*k-1]);
      // - (extremum_alldphi[2*k-1]-extremum_alldphi[2*k-2]); 
      // if (fabs(diff_dphi) > COLL_MIN_DPHI ) {
      // new_y_a = y_a;
      //   new_y_b = y_b - LAMBDA_PHI*diff_dphi; // neg. extremum
      //
      //   NEW VERSION:
       if ((extremum_alldphi[2*k-1] > 1.5707+COLL_MIN_DPHI) && (extremum_alldphi[2*k+1] < 1.5707+COLL_DPHI_EPS)) {
	ident_case_dphi = 0;
        //strategy1: give more amplitude at beginning // NO less positive amplitude
	//  diff_dphi = (extremum_alldphi[2*k-1] - 1.5707);
	//  if (y_a > 0) {
        //    new_y_a = y_a - LAMBDA_PHI*diff_dphi;
	//    new_y_b = y_b;
 	  new_y_a = y_a;
          new_y_b = y_b;
       }
       else if ((extremum_alldphi[2*k-1] > 1.5707+COLL_MIN_DPHI) && (extremum_alldphi[2*k+1] > 1.5707+COLL_MIN_DPHI)) {
	   ident_case_dphi = 1;
	 // to think about: a good idea would be to override predictparams[1] and redefine sub-intervals abcissas
	   diff_dphi = (extremum_alldphi[2*k+1] - 1.5707);
	   new_y_a = y_a;
	   new_y_b = y_b + LAMBDA_PHI*diff_dphi;
       }
       else if ((extremum_alldphi[2*k-1] < 1.5707+COLL_DPHI_EPS) && (extremum_alldphi[2*k+1] > 1.5707+COLL_MIN_DPHI)) {
	 ident_case_dphi = 2;
          // strategy1: give more amplitude at end
         diff_dphi = (extremum_alldphi[2*k+1] - 1.5707);
	 if (y_b < 0) {
	  // new_y_a = y_a + LAMBDA_PHI*diff_dphi; // we could do this but creates mismatch with prv interval
	   new_y_a = y_a;
	   new_y_b = y_b + LAMBDA_PHI*diff_dphi;
	 }
        }
        else if ((extremum_alldphi[2*k-1] < 1.5707+COLL_DPHI_EPS) && (extremum_alldphi[2*k+1] < 1.5707+COLL_DPHI_EPS)) {
	  ident_case_dphi=3;
	// ne rien faire/do nothing
 	  new_y_a = y_a;
          new_y_b = y_b;
        }// endif
    }//endif k%4==1
     
    //for debug:
     std::cout << "Deviation : diff_dphi: " << diff_dphi << ", case: " << ident_case_dphi << ", new_y_a: " << new_y_a << ", new_y_b: " << new_y_b << ", factor: " << LAMBDA_PHI << "\n";

     // write the result in array, but this must be written also on file
     boundary_all[2*k]=new_y_a;
     boundary_all[2*k+1]=new_y_b;
     
     //erase lefft-boundary condition with right end from precedent interval:
     Number new_y_a_mid;
     if (k>0) {
       //new_y_a_mid = (boundary_all[2*k]+boundary_all[2*k-1])/2.0;
       new_y_a_mid = boundary_all[2*k-1];

       boundary_all[2*k-1] = new_y_a_mid;
       boundary_all[2*k] = new_y_a_mid;
    }
   }//endfor
  }//endif enable

return 0;
}


/********** deprecated **********/
/*
bool IDENT05_COLL::UpdateBoundaryIterative()
{
  Number LAMBDA1 = 0.7/predictparams[1];
  Number LAMBDA2 = 0.6/predictparams[1];
  Number LAMBDA_PHI = 0.15/1.57;
  Number COLL_MIN_DPHI = 0.25;
  Number COLL_MIN_DERIVATIVE = predictparams[0]*2*3.14156/predictparams[1]*0.35;
  Number COLL_MAX_MISMATCH = 1.5*COLL_MIN_DERIVATIVE;
  Index k, h, h0, h1, h2, h3, h4; // subinterval no and various indexes: for search (h) and location of extremum of y (h0) and of y' (h1)
  //Number kcoeffs[8];
  Number x_a, x_b, y_b, dy_b, dyb_next, pred_dy_a, pred_dy_b, 
	 y_a, dy_a, x_c, y_c, dy_c, dyc_next, diff_dy_a, 
	 diff_dy_b, new_y_a, new_y_b, new_y_c, diff_dphi; // 
  Number yh, dyh, yh0, dyh0, yh1, dyh1, yh2, dyh2, yh3, dyh3, yh4, dyh4, yM, dyM; 
  Number dphi0, dphi1, dphi2;
  // 
  std::cout << "DEBUG, predictparams " << predictparams[0] << " " << predictparams[1] << " " << predictparams[2] << " " << predictparams[3] << "\n";
  //
  Miter++; // increase of global index of iteration
  for (k=0; k < num_ranges-1; k++) {  // we only evaluate k-1
     // normally import TheCode[k] here, but for a two points config we can test w/o it
     y_a = boundary_all[2*k];
     y_b = boundary_all[2*k+1];
//      y_c = boundary_all[2*k+2];
     dy_a = boundary_derivall[2*k];
     dy_b = boundary_derivall[2*k+1];  //REMARK: this convention for interv k
     dyb_next = boundary_derivall[2*k+2];
//     dy_c = boundary_derivall[2*k+2];
     x_a = times_end[k];
     x_b = times_end[k+1];
//     x_c = times_end[k+1];    
//     we neglect the exponential decay effect for the moment 
     pred_dy_a = 2*3.14156/predictparams[1]*predictparams[0]*sin(2*3.14156*x_a/predictparams[1]+predictparams[3]);
     pred_dy_b = 2*3.14156/predictparams[1]*predictparams[0]*sin(2*3.14156*x_b/predictparams[1]+predictparams[3]);
     //std::cout << "DEBUG predicted param3 " << predictparams[3] << " ypa " << pred_dy_a << " ypb " << pred_dy_b << "\n";
     diff_dy_a = 0; // unused
     diff_dy_b = - dy_b + dyb_next; // same REMARK
     if (k > 0) {
        new_y_a = boundary_all[2*(k-1)+1];  // corresp to new_y_b[k-1], and do nothing if if k==0
     } 
     else {
        new_y_a = boundary_all[0]; // ok to recopy inverted after
     }
     if (fabs(pred_dy_b) > COLL_MIN_DERIVATIVE ) {
	     // fred: here add a modification by authorizing the modification
	     // of y_c ! , always based on slopes difference at x_b and only for
	     // a derivative threshold of dy_b !
	     // fred (15 oct 25): I have changed the sign from minus to plus
        new_y_b = (fabs(diff_dy_b) > COLL_MAX_MISMATCH)? y_b + LAMBDA1*diff_dy_b : y_b;
	//new_y_c = (fabs(diff_dy_b) > COLL_MAX_MISMATCH)? y_c + LAMBDA2*diff_dy_b : y_c;
     }
     else {
       // preview another algorithm
       new_y_b = boundary_all[2*k+1];
     }
     if (Miter > 0) {
	std::cout << "Boundary update, k=" << k << " x_a=" << x_a << ", x_b=" << x_b << "\n";
	std::cout << "Read calc value: ya " << y_a << " yb " << y_b << "\n";
        std::cout << "Predict derivative: dypa " << pred_dy_a << ", dypb " << pred_dy_b << ", min threshd: " << COLL_MIN_DERIVATIVE << "\n";
	std::cout << "Read derivative: dyb " << dy_b << ", dybnext " << boundary_derivall[2*(k+1)] << "\n";
        std::cout << "Deviation derivative: diffdyb " << diff_dy_b << ", threshd: " << COLL_MAX_MISMATCH << ", factors: lambdab=" << LAMBDA1 << ", lambdac=" << LAMBDA2 << "\n";
        std::cout << "new bound value: nya " << new_y_a << " nyb " << new_y_b << "\n";

     }     
     // write the result in array, but this must be written also on file
     boundary_all[2*k]=new_y_a;
     boundary_all[2*k+1]=new_y_b;
  }

// now look for the phase increase
// first thing is to populate the extremum_alllocus, _allvalues tables
  for (k=0; k < num_ranges-1; k++) {
//     if (((k%4 == 0) || (k%4 == 1)) && (boundary_all[2*k] > 0.5*predictparams[0])) { // part decreasing
     if (k%4 == 0) {  // cosinus-like part, decreasing
        h=k*num_points;
	yM=solarray[h][1]; // value
	yh=solarray[h+1][1];
	while ((h < (k+1)*num_points) && (fabs(yM) < fabs(yh))) { // look up for a maximum of y
          h++;
	  yM=yh;
          yh=solarray[h+1][1];
	}  
	h0=h;
	extremum_alllocus[2*k]=h;
	extremum_allval[2*k]=fabs(yM);
	extremum_alldphi[2*k]=0.;
        // reinit
	h=k*num_points;
	dyM=solarray[h][2]; // derivative
	dyh=solarray[h+1][2];
	while ((h < (k+1)*num_points) && ((fabs(dyM) < fabs(dyh)) || (fabs(dyM) < 2.0*COLL_MIN_DERIVATIVE) )) { // look up for a minimum of y'
          h++;
	  dyM=dyh;
	  dyh=solarray[h+1][2];
	  yh=solarray[h+1][1];
	}
	h1=h;
	if (h1==h0) { // security against NaN values
	     h1=h0+1;
       	}
	
	extremum_alllocus[2*k+1]=h;
	extremum_allval[2*k+1]=dyM;
	if (yh > 0) {
	// formula, correct version is line one below: (with scale of phase)
           dphi1 = atan( fabs(dyh)/fabs(yh)/(6.283/predictparams[1]));
	   extremum_alldphi[2*k+1]=dphi1 * predictparams[1]/4./(solarray[h1][0]-solarray[h0][0]);
	//dphi1 = asin(fabs(dyM)/yM/(6.283/predictparams[1]));
	//extremum_alldphi[2*k+1]=dphi1 * predictparams[1]/4./(solarray[h1][0]-solarray[h0][0]);
	//extremum_alldphi[2*k+1]=(asin(fabs(dyM)/yM/(6.283/predictparams[1])))*predictparams[1]/4./(solarray[h1][0]-solarray[h0][0]);
        //extremum_alldphi[2*k+1]=asin(fabs(dyM)/(fabs(yM)*6.28/predictparams[1])) ;

	}
	else { // yh<0
	// formula, correct version is line one below: (with scale of phase)
           dphi1 = 1.5707+atan( fabs(yh)*(6.283/predictparams[1])/fabs(dyh));
	   extremum_alldphi[2*k+1]=dphi1 * predictparams[1]/4./(solarray[h1][0]-solarray[h0][0]);

	}
	std::cout << "DEBUG JUST AFTER phi" << 2*k+1 << " " << dphi1 << " calc: yh1=" << yh << " dyh" << dyh  << " t1=" << solarray[h0][0] << " t2=" << solarray[h1][0] << "\n";

     }
//    else if (((k%4 == 0) || (k%4 ==1)) && (fabs(boundary_all[2*k]) < 0.2*fabs(predictparams[0]) )) { // sinus increasing or cosinus decreasing part
     else if (k%4 == 1) {  // sinus-like part
        h=num_points*k;
	dyM=solarray[h][2];
	dyh=solarray[h+1][2];
        while ((h < (k+1)*num_points) && ((fabs(dyM) < fabs(dyh)) || (fabs(dyh) < 2.0*COLL_MIN_DERIVATIVE)) ) { // look up y'(max)
           h++;
	   dyM=dyh;
	   dyh=solarray[h+1][2];
	}
	h1=h;
	dyh1=solarray[h1][2];
	yh1=solarray[h1][1];
	extremum_alllocus[2*k]=h;
	extremum_allval[2*k]=dyM;
	//extremum_alldphi[2*k+1]=asin(fabs(dyM)/yM/(6.283/4./(solarray[h1][0]-solarray[h0][0])) ); // fred: leave here as comment since another formula is used below
        h=num_points*k;
	yM=solarray[h][1];
	yh=solarray[h+1][1];
	while ((h< (k+1)*num_points) && ( ((fabs(yM) < fabs(yh)) && (yh<0)) || ((fabs(yM) > fabs(yh)) && (yh>0)) )) { // look up for a minimum of yi, FRED: OK!
          h++;
	  yM=yh;
          yh=solarray[h+1][1];
	}
	h2=h;
	dyh2=solarray[h2][2];
	yh2=solarray[h2][1];
	
	extremum_alllocus[2*k+1]=h2;
	extremum_allval[2*k+1]=fabs(yM);
	if (dyh2 > 0) {  // corrected here: sign(dyh2) test
	   // FRED: la formule ci-dessous qui est la meilleure est la premiere TBD
	   extremum_alldphi[2*k+1]=(3.1415-atan(fabs(dyh2)/fabs(yh2)/(6.28/predictparams[1])))*(6.28/4./(solarray[h2][0]-solarray[h1][0]));
	   // FRED : essayer en prenant en compte les abcisses pour 'mesurer' la dilatation 
	  //dphi2 = (3.1415-asin( fabs(yh2)*(6.28/predictparams[1])/( fabs(dyM)) ));
	  //extremum_alldphi[2*k+1]=dphi2*(predictparams[1]/4./(solarray[h2][0]-solarray[h1][0]));
	   // FRED : this might be OK but prefer the following 
	  //extremum_alldphi[2*k+1]=3.1415-asin( fabs(yh2)*(6.28/predictparams[1])/( fabs(dyM)) );
	} 
	else {
	   // FRED: la formule ci-dessous qui est la meilleure est la premiere TBD
	   dphi2 = atan(fabs(dyh2)/fabs(yh2)/(6.28/predictparams[1]));
	   extremum_alldphi[2*k+1]=dphi2*(6.28/4./(solarray[h2][0]-solarray[h1][0]));
           //extremum_alldphi[2*k+1]=atan(fabs(dyh2)/fabs(yh2)/(6.28/4./(solarray[h2][0]-solarray[h1][0])));
	  //dphi2 = asin( fabs(yh2)*(6.28/predictparams[1])/( fabs(dyM)) );
	  //extremum_alldphi[2*k+1]=dphi2*(predictparams[1]/4./(solarray[h2][0]-solarray[h1][0]));
	   //extremum_alldphi[2*k+1]=asin( 0.9*fabs(yh2)*(6.28/predictparams[1])/(fabs(dyM)) );
	   std::cout << "DEBUG JUST AFTER phi" << 2*k+1 << " " << dphi2 << " calc: yh2=" << yh2 << " dyM" << dyh1  << " t1=" << solarray[h1][0] << " t2=" << solarray[h2][0] << "\n";

	} //endif 1
     } // endif 2   
   //    if (((k%4 ==2) && (boundary_all[2*k] < 0.5*predictparams[0])) || ((k%4==0) && (boundary_all[2*k] > 0.5*predictparams[0])) ) { // cosinus increase or sinus increasing
     if (k%4 == 2) {  // cosinus-like part, increasing
        h=k*num_points;
	yM=solarray[h][1]; // value
	yh=solarray[h+1][1];
	while ((h < (k+1)*num_points) && (fabs(yM) < fabs(yh))) { // look up for a minimum of y (max(fabs))
          h++;
	  yM=yh;
          yh=solarray[h+1][1];
	}  
	h3=h;
	extremum_alllocus[2*k]=h3;
	extremum_allval[2*k]=fabs(yM);
	extremum_alldphi[2*k]=0.;
        // reinit
	h=k*num_points;
	dyM=solarray[h][2]; // derivative
	dyh=solarray[h+1][2];
	while ((h < (k+1)*num_points) && (fabs(dyM) < fabs(dyh))) { // look up for a maximum of y'
          h++;
	  dyM=dyh;
	  dyh=solarray[h+1][2];
	}
	h4=h;
	if (h4==h3) { // security against NaN values
	     h4=h3+1;
       	}
	extremum_alllocus[2*k+1]=h4;
	extremum_allval[2*k+1]=dyM;
	// formula, correct version is line two below: (difference of phase)
	// extremum_alldphi[2*k+1]=asin(fabs(dyM)/yM/(6.283/4./(solarray[h4][0]-solarray[h3][0])) );
	//extremum_alldphi[2*k+1]=asin(fabs(dyM)/(fabs(yM)*6.28/predictparams[1])) *predictparams[1]/4./(solarray[h4][0] - solarray[h3][0]);

        extremum_alldphi[2*k+1]=asin(fabs(dyM)/(fabs(yM)*6.28/predictparams[1])) ;

     } // endif1  
  }// end for k

  // now re iterate loop over k intervals, this time knowing the max locus
  
  for (k=0; k < num_ranges-1; k++) {
     std::cout << "Phase increase, ph"<< 2*k << "=" << extremum_alldphi[2*k] << ", ph" << 2*k+1 << "=" << extremum_alldphi[2*k+1] << "\n" ;
  }

  // now the update
  int enable_new_dphi=0;
if ( enable_new_dphi == 1 ) {
  for (k=1; k < num_ranges-1; k++) {
    if (k%4==1) {
     y_a=boundary_all[2*k];
     y_b=boundary_all[2*k];
      // compare phase increase from previous interval
       diff_dphi = (extremum_alldphi[2*k+1]-extremum_alldphi[2*k-1]);
      // - (extremum_alldphi[2*k-1]-extremum_alldphi[2*k-2]); 

       if (fabs(diff_dphi) > COLL_MIN_DPHI ) {
	 new_y_a = y_a;
         new_y_b = y_b - LAMBDA_PHI*diff_dphi; // neg. extremum
     }

     //for debug:
     std::cout << "Deviation : diff_dphi: " << diff_dphi << ", new_y_b: " << new_y_b << ", factor: " << LAMBDA_PHI << "\n";

     // write the result in array, but this must be written also on file
     boundary_all[2*k]=new_y_a;
     boundary_all[2*k+1]=new_y_b;
    }
  }
}
  return 0;
}
*/
