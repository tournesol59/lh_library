#another version with smarter functions that test equality of points and
#proven algorithm
import numpy as np
import math

NCOMP=2 #number of first order state or d
M=[2,2] #vector number of degree max derivative of each state
KP=3 # number of collocation points in each interval
MSTAR=4 # sum,i m(i)

#adapt new mesh from previous values, fixpoints or user defined
# to do: preview compute and save in VALS the values of previous solution
ALEFT=0.1
ARIGHT=0.8
#Ncomp=2
NPOLD=8
NP=15
NFIXP=2
FIXPNT=[0.15, 0.5]
XIOLD=[0.1, 0.15, 0.2, 0.4, 0.5, 0.6, 0.75, 0.8]
XI=[]
PRECIS=1e-9

# anticipating solving: structures that contain the solution
ZVALS=np.zeros((NCOMP,NP*KP))  # values of first orders states
DMZVALS=np.zeros((MSTAR,NP*KP))  # values of derivatives of states i: l=0,m(i)
# WGTMSH(l) are precomputed CONSTS2(k+mj,j) divided by the relative TOL_l, where
# j=JTOL(l) is the index in the state vector z_j(u)=(u_1, u^(1)_1 ...
# k=KP is the collocation order and l the appropriate tolerance index
WGTMSH=np.zeros((NCOMP*MSTAR)) # should be of same dimension as state vector
CONSTS2=np.zeros((28))  # values taken from Asher
CONSTS2=[0.125, 2.604E-3, 8.019E-3, 2.17E-5, 
          7.453E-5,  5.208E-4,  9.689E-8,  3.689E-7,  3.100E-6,
          2.451E-5,  2.691E-10, 1.120E-9,  1.076E-8,  9.405E-8,
          1.033E-6,  5.097E-13, 2.290E-12, 2.446E-11, 2.331E-10,
          2.936E-9,  3.593E-8,  7.001E-16, 3.363E-15, 3.921E-14,
          4.028E-13, 5.646E-12, 7.531E-11, 1.129E-9  ]

#COEF = []
  # run the colnew program to output RK coeff one time

#****************************************************
# mesh halving:
def testpointequal(x1,x2):
#    if (((x1-x2)<PRECIS) and ((x2-x1)>-PRECIS)) or (((x2-x1)<PRECIS) and ((x1-x2)>-PRECIS)) :
    if (math.fabs(x1-x2)<PRECIS):
        return 1
    else:
        return 0


def newmesh(xi,xiold,fixpnt,Np,Npold,Nfixpnt,aleft,aright):
   #remark: this program must be tested because it will be much used 
#  for test:
#   Npold=NPOLD #
#   Nfixpnt=NFIXP #
#   Np=2*Npold-1
#   aleft=ALEFT
#   aright=ARIGHT
   for i in range(0,Np): #init
      xi.append(0.0)
   i=0 # index xiold
   k=0 # index xi
   for j in range(0,Nfixpnt): #esperons pas de debordement
      print("#.j= "+str(j)) #debug
      if (j==0):
         xleft=aleft  # searched pnts in interval [aleft, Nfixpnt]
         xright=fixpnt[j]
      elif (j<=(Nfixpnt-1)): # simpler 
         xleft=fixpnt[j-1]
         xright=fixpnt[j]
      else: # (j==(Nfixpnt)):
         xleft=fixpnt[Nfixpnt-1]
         xright=aright
 
      while ((i<Npold-2) and (xiold[i+1] < xright-PRECIS)): 
      #it will go along the table xiold except two config with last 3 nodes
         dx=xiold[i+1]-xiold[i] # which must be halved and whever fix position
         print("#...A :"+str(i)+" "+str(k)) #debug
         print("  ... i+1 "+str(i+1)+" xright "+str(xright)+" xiold "+str(xiold[i+1]))
         xi[k]=xiold[i]
         k=k+1
         xi[k]=xiold[i]+dx/2.
         k=k+1
         i=i+1

      if (i<Npold-2) and (xiold[i+1]>=xright):
      #we are always betwenn two fixed points but next step, nextfixed interval
         if (testpointequal(xiold[i+1],xright)==1):
            xi[k]=xiold[i]
            dx=xiold[i+2]-xiold[i]
            if ((xiold[i]+dx/2)<xright-PRECIS):
               print("#...B :"+str(i)+" "+str(k)) #debug
               k=k+1
               xi[k]=xiold[i]+dx/2.
               k=k+1
               xi[k]=xright
               k=k+1
               i=i+2
            elif ((xiold[i]+dx/2)>xright+PRECIS):
               print("#...C :"+str(i)+" "+str(k)) #debug
               k=k+1
               xi[k]=xright
               k=k+1
               xi[k]=xiold[i]+dx/2.
               k=k+1
               i=i+2
            else: # coincidence of xiold[i]+dx/2 and xright
                print("#...D :"+str(i)+" "+str(k)) #debug
                k=k+1
                xi[k]=xright
                k=k+1
                i=i+2 #whereas only one pnt inserted in new, next: xiold[i+2]
        #endif xiold[i+1]==xright but seems always good
      elif (i>=Npold-2): # caution indent
         if (fixpnt[j]<xiold[i]) and not (testpointequal(fixpnt[j],xiold[i])):
             #WAS: fixpnt[j]<xiold[i-1]
    # tip: always give an initial mesh that has no fixed pnt in this pos
            print("#...Q :"+str(i)+" "+str(k)) #debug
            break # quit j loop

   #after the j loop: two configurations
   while (i<=Npold-2) and (fixpnt[Nfixpnt-1]<xiold[i]-PRECIS): # second part only necessary in first loop iteration
      print("#...copy A :"+str(i)+" "+str(k)+" "+str(xiold[i])) #debug
      # recommend: always give an initial mesh that has fixed pnt in this pos
      dx=xiold[i+1]-xiold[i]
      xi[k]=xiold[i]
      k=k+1
      xi[k]=xiold[i]+dx/2.
      k=k+1
      i=i+1
   xi[k]=xiold[i] #terminaison
   if (i==Npold-3) and (testpointequal(fixpnt[Nfixpnt-1],xiold[i+1])==1):
      dx=aright-xiold[Npold-3]
      if (xiold[i]+dx/2 < fixpnt[Nfixpnt-1]-PRECIS):
         print("#...copy B :"+str(i)+" "+str(k)+" "+str(xiold[i])+" "+str(xiold[i+1])) #debug
         xi[k]=xiold[i]+dx/2.
         k=k+1
         xi[k]=fixpnt[Nfixpnt-1]
         k=k+1
         i=i+2
         xi[k]=xiold[i] #terminaison
      if (xiold[i]+dx/2 > fixpnt[Nfixpnt-1]+PRECIS):
         print("#...copy C :"+str(i)+" "+str(k)+" "+str(xiold[i])+" "+str(xiold[i+1])) #debug
         #xi[k]=fixpnt[Nfixpnt-1] # NO!
         #k=k+1
         xi[k]=xiold[i]
         k=k+1
         xi[k]=xiold[i]+dx/2.
         k=k+1
         i=i+2
         xi[k]=xiold[i] #terminaison
      else: # coincidence of xiold[i]+dx/2 and xright
         print("#...copy D :"+str(i)+" "+str(k)) #debug
         xi[k]=xright
         k=k+1
         i=i+2 
         xi[k]=xiold[i] #terminaison
   print("NP is "+str(k+1)+"\n")
   #end, completed

#**********************************************
#aux subroutines
LTOL=[1,2,1,2]
JTOL=[1,1,2,2]

def meshhorder(i,uhigh,hi,dmz,ncomp,k):
# determine highest order (piecewise constant) derivatives
# of the current collocation solution
   dn=1.0/h**(k-1)
   uhigh=np.zeros((NCOMP))
   kin=1
   idmz=(i-1)*k*ncomp+1
   for j in range(0,k-1):
      fact=dn*coeff[kin]
      for id in range(0,ncomp-1):
         uhigh[id]=uhigh[id]+fact*dmz[idmz]
         idmz=idmz+1
      kin=kin+k

def meshapprox(j,TOLj,zvalues,xi,xiold,accum): # think about other args
# compute the rooted accumulated value for the old mesh and then the halved mesh
   accum=np.zeros((1,Npold))
#xleft=aleft
#xright=aright

#init weight table
#kp=3
   koff=kp*(kp+1)/2
   iz=0
   for i in range(0,NCOMP):
      for l in range(0,KP):
         WGTMSH[i*KP+l]=CONSTS2[koff-m[i]+iz] / TOLj[i*KP+l] # weight
         ROOT[i*KP+l]=1.0/(KP+m[i]-LTOL[i*KP+l])  # reciprocal exponent

# init 
   accum[0]=0.0
   for i in range(0,Npold): #loop over the mesh
      # init
      accleft=0.0
      # select weight and root
      wgt=WGTMSH[j]
      root=ROOT[j]
    
