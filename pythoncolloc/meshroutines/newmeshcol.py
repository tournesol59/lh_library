	#adaptation of newmesh from colnew.f
import packcollocpts as ppts
import numpy as num
import meshfunc as meshfc

#plan: 1) define methods meshhalve and meshuniform and which variable they share

CONSTS2=num.zeros((28))  # values taken from Asher
CONSTS2=[0.125, 2.604E-3, 8.019E-3, 2.17E-5, 
          7.453E-5,  5.208E-4,  9.689E-8,  3.689E-7,  3.100E-6,
          2.451E-5,  2.691E-10, 1.120E-9,  1.076E-8,  9.405E-8,
          1.033E-6,  5.097E-13, 2.290E-12, 2.446E-11, 2.331E-10,
          2.936E-9,  3.593E-8,  7.001E-16, 3.363E-15, 3.921E-14,
          4.028E-13, 5.646E-12, 7.531E-11, 1.129E-9  ]


PRECIS=1e-6
class newmesh(object):
    def __init__(self, mstar,n,nold,nfix,xi,xiold,Z,DMZ,accum,fixpnt):
        self.__aleft=xiold[0]
        self.__aright=xiold[nold-1]
        self.__mstar=mstar
        self.__n=n
        self.__nold=nold 
        self.__nfix=nfix
        self.__xi=num.copy(xi)
        self.__xiold=num.copy(xiold)
        self.__Z=num.copy(Z)
        self.__DMZ=num.copy(DMZ)
        self.__accum=num.copy(accum)
        self.__fixpnt=num.copy(fixpnt)

#end of class newmesh data definition 

    def meshhalve(self): 
    # intern data manip:
    # xi, xiold, fixpnt, slope, aleft, aright 
    # extern data manip: 
    # k, order of sheme 
    # valstr copy of approx at some subinterval end points
    #remark: this program is far more complicated as needed for halving
    # but some parts will be used to program the other part (mesh select)
        Nold=self.__nold
        Nfixpnt=self.__nfix
        Np=2*Nold
        aleft=self.__aleft
        aright=self.__aright

        #init
        self.__xi=[] #.reshape((Np)) # or free + append
        for i in range(Np):
            self.__xi.append(0.0)
        # call to extern function newmesh
        meshfc.newmesh(self.__xi,
                  self.__xiold,
                  self.__fixpnt,
                  self.__n,
                  self.__nold,
                  self.__nfix,
                  self.__aleft,
                  self.__aright)

        if (Np==self.__n):
           print("newmesh: same number"+str(Np))
        else:
           print("newmesh: mismatch possibly due to fixpnt"+str(Np))
        self.__n=Np
 
#end of class method meshhalve

    def meshuniform(self, Nintval):
    # caution: unlike meshhalve, this should reinitialize the mesh as uniform
    # intern data manip:
    # xi, xiold, fixpnt, slope, aleft, aright
        
        Nfixpnt=self.__nfix
        N=2*Nold-1
        aleft=self.__aleft
        aright=self.__aright     # denote N=2+Nfixpnt+order
        self.__n=2+Nfixpnt+Nintval-1
   # loop over fixpoint k:
        self.__xi[0]=aleft
        i=0  #to order
        k=0
        step=(aright-aleft)/Nintval
        self.__xi[i]=aleft
        xp1=self.__xi[i]+step
        i=i+1
        while (k<=Nintval):
            if (xp1 < self.__fixpnt[k]):
                self.__xi[i]=xp1
                i=i+1
                xp1=self.__xi[i]+step
            elif (k<Nfixpnt):
                 self.__xi[i]=self.__fixpnt[k]
                 k=k+1
                 i=i+1
            else:
                print("meshuniform: error fixpnt overcount aright bound")

        if (self.__n != i):
            print("meshuniform: error about number of points "+str(i))
#end of class method meshnonuniform

    def meshprint():
       print("A new mesh has "+str(self.__n)+" points")
       print(str(self.__xi))

    def calcslope(pts,Z,DMZ,order,mstar,i,M,j,ltol,mn,slope):
   # algorithm for giving a criterion of error-reduced to hi 
   # a tolerance tol(j) is applied: if ltol(j)=l then tol is indicated to
   # bound the error related to the l-th component of |z(u)-z(v)| wh. v
   # is the approximation of the solution and operator z() gives its states
   #
   # |zj(u)-zj(v)|     < tol * fabs(z(u))
   #              ltol                  ltol

   # mstar is the number of states in solution vector
   # i is the intervals index (up Np-1 in the denomination of mesh)
   # j is the index of tolerance
   # n is the index of the state, whose m(n) derivative is subject to tol
   # (from theory:
   # |zj(u)-zj(v)|     < (coeff(order,order-vu))*|z^(k+mn)(u)|*hi^(k+mn-l)
   #             ltol
   # where mn is the current state max order of derivative in the system
   # and vu=mn-l
   # IMPORTANT ! since the derivative of a k-solution of order k+mn is unknown
   # and should be zero, a solving of order 2*k > k+mn will be performed
   # in the solver prgm (see pendulum.py)

   # store current interval (i) result
        deriv=num.zeros((order,M+1))
        maxderiv=math.fabs(DMZ[0])
        for r in range(0,order):
            for l in range(0,M+1):
                deriv[r,l]=DMZ[r*(M+1)+l]
                if (math.fabs(deriv[r,l]) > maxderiv):
                    maxderiv=math.fabs(deriv[r,l])
       # get the coeff correspondg. to mn,j
        coeff=CONSTS2[order*4+(order-mn+ltol)]
       # compute with formula above
        rootj=1/(order+mn-ltol)
        slope=(coeff/1e-6*maxderiv)**(rootj)
        return slope

#end of class newmesh definition


