#calculer, dans le cas d'un pendule double avec un amortissement et une force du 
#ressort, ou bien d'une équation de Mathieu au choix  la trajectoire temporelle pour différentes conditions initiales de X(0) (R4, ou bien R2) et 
#comparer différents cas de parametres (m1,m2,k1,mu1) ou bien (a,q) avec linearisation autour de 
#Theta10#0 et Theta20=0

import numpy as num
import packcollocpts as ppts 

class LinearEqnSolver(object)
# contains the definition of the (linear) system to solve
# the names of the function "fun" and "jac" must be generic
# as they are used in other part in the program

    def __init__(self,N,M,X0,range,params):
         #param is a tuple with (m1,m2,l1,l2,mu1,mu2,k)
         self.__X0= num.zeros((4))
         for i in range(0,4):
             self.__X0[i]=X0[i]
         self.__ni=num.zeros((N+1))
         self.__M=M # order
         self.__N=N # no interval
         for i in range(0,N+1):
             self.__ni[i]=i*range/(N)
         self.__m1=m1
         self.__m2=m2
         self.__l1=l1
         self.__l2=l2
         self.__mu1=mu1
         self.__mu2=mu2
         self.__k=k
         self.__J11=1/3.0*m1*l1**2+m2*l1**2+m2*l1*l2+m2*l1*l2+1/3.0*m2*l2**2
         self.__J22=1/3.0*m2*l2**2+1/2.0*m2*l1*l2
         self._J12=1/3.0*m2*l2**2+1/2.0*m2*l1*l2
         self.__J21=1/6.0*m2*l2**2
# to complete
   def setallpoints(self):
        self.__ti=num.zeros((self.__N*self.__M+1))

        for i in range(0,self.__N):
            for l in range(0,self.__M):
# get the quadrature/collocation points in interval [-1,1] and scale them to [ti,ti+1]
                nu=ppts.diff_nodes[self.__M,l]
                self.__ti[i*self.__M+l]=(self.__ni[i]+self.__ni[i+1])/2+(self.__ni[i+1]-self.__ni[i])/2*nu

   def setmatrixV(self):
#Tbc plus comments

      #    A=[ 0,   1,   0,   0]
      #-1/J1[K1, mu1, K12, mu2s]
     #         [0,   0,   0,   1]
     #-1/J2[K21,mu1s,K2,mu2]

   def setmatrixW(self):
