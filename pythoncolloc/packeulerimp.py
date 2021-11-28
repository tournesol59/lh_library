# oser un Shema implicite pour avoir le 
# point milieu, puis le point t+1 pour avoir 
# la tangente comme CI de collocation

import numpy as num
#import packcollocpts as cpts 
import math as math 

def f(y):
   return 1+(math.tan(y))**2

def df(y):
   return 2*(math.tan(y))*(1+y**2)

class pkeulerImp(object):
    #takes a vector with aquidistant pts
    def __init__(self,xi,yi,nx):
        self.__nx=nx
        self.__dx=xi[1]-xi[0]
        self.__xi=num.zeros((nx))
        self.__yg=num.zeros((nx))
        self.__y=num.zeros((nx))
        self.__y1=num.zeros((nx))
        self.__diag=num.zeros((nx))
        self.__lowd=num.zeros((nx))
        self.__rhs=num.zeros((nx))
        for i in range(0,nx):
           self.__xi[i]=xi[i]
           self.__yg[i]=yi[i]
           self.__y[i]=yi[i]

     def __setimpl__(self):
           nx=self.__nx
           self.__diag[0]=1-self.__dx*df(self.__y[0])
          self.__rhs[0]=self.__yg[0]-self.__y[0]+self.__dx*f(self.__y[0])
           for i in range(1,nx-1):
              self.__diag[i]=1-self.__dx*df(self.__y[i])
              self.__lowd[i]=-1
              self.__rhs[i]=(-self.__y[i]+self.__y[i-1])+self.__dx*f(self.__y[i])

     def __lsor__(self):
         omega = 0.5
#lsor algorithm for solving sparse 2 diags system 
         for i in range(0,i):
            sum =0 
            for j in range(0,i-1):
          sum=sum+self.__lowd[j]*self.__y[j]*(1-om)-om*self.__lowd[j]*self_y1[j]
         self._y1[i]=1/self.__diag[i]*(self.__rhs[i]+sum)

    def __solve__(self):
       k=0
       res=0 
       while (k<11) and (res>0.01):
           self.__lsor__()
           res=0.0
           for i in range(0,self.__nx): 
                 res=res+math.abs(self.__y1[i]-self.__y[i])
                 self.__y[i]=self.__y1[i]

    def __str__(self):
       nxe =self.__nx
       for i in range(0,nxe):
           print(str(self.__y1[i])

#main
Nx=10
Xi=num.linspace(0.0, 1.0, Nx)
Yi=0.8*math.tan(Xi)
Exeuler=pkeulerImp(Xi,Yi,Nx)
Exeuler.__setimpl__()
Exeuler.__solve__()
Exeuler.__str__()
