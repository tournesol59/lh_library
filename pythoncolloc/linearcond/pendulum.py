Pythoncolloc\pendulumlinear.py


#calculer, dans le cas d'un pendule double avec un amortissement et une force du 
#ressort, ou bien d'une équation de Mathieu au choix  la trajectoire temporelle pour différentes conditions initiales de X(0) (R4, ou bien R2) et 
#comparer différents cas de parametres (m1,m2,k1,mu1) ou bien (a,q) avec linearisation autour de 
#Theta10#0 et Theta20=0

import numpy as num
import packcollocpts as ppts 

#----------------------------------------------------------
#initialize some global structures every prog component should have access to: 
class ApproxRes(object)
   def __init__(self,nc,ms,N,M):
       # N is no of intervals
       # M is collocation order: M+1 pts
       mstar=0
       for i in range(0,nc):
           mstar=mstar+ms[i]
    #Z is a structure with the solutions at the legendre points hence size=N*(M+1)
       self.__Z=num.zeros(N*(M+1))
#DMZ is a structure with all the derivative up to ms[] of the solution at the legendre
#points hence size=mstar*N*(M+1)
       self.__DMZ=num.zeros(N*(M+1)*mstar)
#Approx is a structure with the special derivatives y^{M+ms} at the legendre points
# hence size=N*(M+1)
       self.__Approx=num.zeros(N*(M+1))

   def approx(self,i,xi,x):
       # perform a Lagrange interpolation with legendre pts values in Z
       xext=num.zeros((M+2)) #M+1+1 because the pqrticular value of x is attached
       for l in range(0,M+1):
           xext[l]=xi[l]
       xext[M+1]=x
       aMatLagr=MatLagrKernel(xext,M)
       KM=aMatLagr.__KM[:,:]
       yval=0.0 # sum of evaluation of Lagr. polynoms
       for l in range(0,M+2):
           prodx=self.__Z[l]
           for j in range(0,M+1):
               if (l~=j):
                   prodx=prod*KM(l,j)
            yval=yval+prodx
        return yval

# end of class ApproxRes

#--------------------------------------------
class LinearEqnSolver(object)
    def __init__(self,nc,ms,N,M,X0,trange,params):
         #param is a tuple with (m1,m2,l1,l2,mu1,mu2,k)
         #nc is the number of subsystem
         self.__nc=nc 
         #ms is an array of size nc, ms[i] = order of subsystem i thus
         mstar=0
         for i in range(0,nc):
            mstar=mstar+ms[i]
         self.__mstar=mstar
         self.__X0= num.zeros((4))
         for i in range(0,4):
             self.__X0[i]=X0[i]
         self.__ni=num.zeros((N+1))
         self.__M=M # order of collocation
         self.__N=N # no intervals thus N+1 node points

         for i in range(0,N+1):
             self.__ni[i]=i*trange/(N)
         m1=params[0]
         m2=params[1]
         l1=params[2]
         l2=params[3]
         mu1=params[4]
         mu2=params[5]
         k=params[6]
         self.__m1=params[0]
         self.__m2=params[1]
         self.__l1=params[2]
         self.__l2=params[3]
         self.__mu1=params[4]
         self.__mu2=params[5]
         self.__k=params[6]

         self.__J1=1/3.0*m1*l1**2+m2*l1**2+m2*l1*l2+1/3.0*m2*l2**2
         self.__J2=1/3.0*m2*l2**2+1/2*m2*l1*l2
         self.__J12=1/3*m2*l2**2+1/2*m2*l1*l2
         self.__J21=1/6*m2*l2**2

# get the quadrature/collocation points in interval [-1,1] and scale them to [ti,ti+1]

   def setallpoints(self):
        self.__ti=num.zeros((self.__N*self.__M+1))
        for i in range(0,self.__N):
            for l in range(0,self.__M):
                nu=ppts.diff_nodes[self.__M,l]
                self.__ti[i*self.__M+l]=(self.__ni[i]+self.__ni[i+1])/2+(self.__ni[i+1]-self.__ni[i])/2*nu

   def setmatrixW(self,i):
      #skeleton of this function which init the spectral derivative matrix
      #and uses it to compute the linear system to solve the diff equations
      #get access to mstar and M
      mstar=self.__mstar
      M=self.__M
      params=(M,0.0,0.5,3) # time interval [ti,ti+1] shall be decomp by 3
      params(1)=self.__ti[i*self__M]
      params(2)=self.__ti[(i+1)*self.__M]
      xi=num.zeros(M-1)
      for l in range(0,M-1):
          xi[l]=self.__ti[i*self__M+l+1]  #these are interior nodes
      #init a class MatDiffPoints(param,x) with x=mesh and size M
      aMatDiff=ppts.MatDiffPoints(param,xi)
      DM1=aMatDiff.__DM[1,:,:]

     #    A=[ 0,   1,   0,   0]
     #-1/J1 [K1, mu1, K12, mu2s]
     #      [0,   0,   0,   1]
     #-1/J2 [K21,mu1s,K2,mu2]
      
      # V=(DM_-A^) where
      # DM_ is DM1 copied each/two lines: one line for deriv st1, one for 2nd deriv st1, ...
      # A^ is A block copied as many times as colloc points


   def setmatrixGam(self,i):
     # assembly matrix Gam=
     # [1 0 0 0     ] [x1] = (y0)  wh. y0 is the boundary value
     # [0W1 0 0 0 ..  [x1] = (b1-a1) wh. b1=rhs, a1=A[:,end]*y1 wh. y1 is a guess value
     # [0 0 W2 ...  ] [x2] = (b2-a2)
     # ...
     # [0 0 0 ... WN] [xN] = (bN-aN)
     # this is verified by methods of Weidemann is matlab

