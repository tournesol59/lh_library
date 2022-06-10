	#adaptation of newmesh from colnew.f
import packcollocpts as pkcol
import numpy as np
#aleft=0.1
#aright=0.8
#Ncomp=1
#NPOLD=8
#NP=15
#NFIXP=3
#KWARGS={'Np':NP, 'Npold':NPOLD, 'Nfixpnt':NFIXP}
#FIXPNT=[0.15, 0.5, 0.75]
#XIOLD=[0.1, 0.15, 0.2, 0.4, 0.5, 0.6, 0.75, 0.8]
#XI=[]
#PRECIS=1e-9

#plan: 1) define methods meshhalve and meshuniform and which variable they share
PRECIS=1e-6
class newmesh(object):
    def __init__(self, mstar,n,nold,nfix,xi,xiold,Z,DMZ,accum,fixpnt):
        self.__aleft=xiold[0]
        self.__aright=xiold[nold-1]
        self.__mstar=mstar
        self.__n=n
        self.__nold=nold 
        self.__nfix=nfix
        self.__xi=np.copy(xi)
        self.__xiold=np.copy(xiold)
        self.__Z=np.copy(Z)
        self.__DMZ=np.copy(DMZ)
        self.__accum=np.copy(accum)
        self.__fixpnt=np.copy(fixpnt)

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
        N=2*Nold
        self.__n=N
        aleft=self.__aleft
        aright=self.__aright
        #init
        self.__xi.reshape((N)) # or free + append 
        i=0 # xiold mesh index
        k=0 # xi mesh index
        # now start follow the xiold and fixpnt:
        for j in range(0,Nfixpnt):
        #    xleft=fixpnt[0]
        #    xright=fixpnt[1]
        #before first fix point:
            if (j==0):
                xleft=aleft  # searched pnts in interval [aleft, Nfixpnt]
                xright=self.__fixpnt[j]
            if (j==(Nfixpnt-1)):
                xleft=self.__fixpnt[Nfixpnt-1]
                xright=aright
            else: 
                xleft=fixpnt[j]
                xright=fixpnt[j+1]
        #tip : fix points are part of the mesh halving
            while (i < Nold-1) and (self.__xiold[i] < xright-PRECIS): 
                self.__xi[k]=self.__xiold[i]
                k=k+1
                if (i < Nold-2):  # perhaps this case could be 
                          # part of the while
                    dx=self.__xiold[i+1]-self.__xiold[i]
                    self.__xi[k]=self.__xi[k-1]+dx/2
                    k=k+1
                    i=i+1  #for treating xiold(i+1) next i
                else: 
                    break
        #end while current fix point
        #preview end
            if (j>=Nfixpnt-1) and (i < Nold-2): #test last fix point
                self.__xi[k]=xleft
                k=k+1
                dx=self.__xiold[i+1]-self.__xiold[i]
                self.__xi[k]=self.__xi[k-1]+dx/2
                #k=k+1 #not needed 
                i=i+1
                break
            elif (i==Nold-2): #test last point if it is not a fix poitn
                self.__xi[k]=self.__xiold[i]
                k=k+1
                dx=self.__xiold[i+1]-self.__xiold[i]
                self.__xi[k]=self.__xi[k-1]+dx/2
                k=k+1 
                self.__xi[k]=aright
                break
            elif (i==Nold-1):
                break #nothing else
 #end for
 
        if (N==k-1):
            print("newmesh: same number"+str(N))
        else:
            print("newmesh: mismatch"+str(N))
    #args['Np']=Np

#end of class method meshhalve

    def meshnonuniform(self, order):
    # caution: unlike meshhalve, this should take Legendre points instead
    # intern data manip:
    # xi, xiold, fixpnt, slope, aleft, aright
        Nold=self.__nold
        Nfixpnt=self.__nfix
        N=2*Nold-1
        aleft=self.__aleft
        aright=self.__aright     # denote N=2+Nfixpnt+order
        self.__n=2+Nfixpnt+order
   # loop over fixpoint j:
        self.__xi[0]=aleft
        i=0  #to order
        k=1
        xleft=aleft + (pkcol.diff_nodes[order,i]+1)/2.*(aright-aleft)
   # xleft shall be the far most advanced node <= fixpnt[j]
        for j in range(0, Nfixpnt):
            while (i<=order) and (xleft<=self.__fixpnt[j]-PRECIS):
                self.__xi[k]=xleft
                k=k+1
                i=i+1
                if (i>order):
                    print("too much points")
                    break
                else:
                    xleft=aleft + (pkcol.diff_nodes[order,i]+1)/2.*(aright-aleft)
            self.__xi[k]=self.__fixpnt[j]
            k=k+1
        self.__xi[k]=aright
        k=k+1
        if (self.__n != k):
            print("error about number of points "+str(k))
#end of class method meshnonuniform


#end of class newmesh definition
MSTAR=2
N=12
NOLD=6
NFIX=1
FIXPNT=[5.5]
XIOLD=[0.,1.5,3.,6.,9.,10.5,12.]
XI=[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
Z=[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
DMZ=[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
ACCUM=[0.]

instNewmesh=newmesh(MSTAR,N,NOLD,NFIX,XI,XIOLD,Z,DMZ,ACCUM,FIXPNT)

#test:
for i in range(0,NFIX):
    print(FIXPNT[i])
print(" ")
for i in range(0,NOLD):
    print(XIOLD[i])

instNewmesh.meshhalve()
print("Np= "+str(N))
for i in range(0,N):
    print(XI[i])


