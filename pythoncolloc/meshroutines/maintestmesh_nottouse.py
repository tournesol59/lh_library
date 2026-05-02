## test program, to comment if package is vaidated:
import math as math 
import numpy as num
import packcollocpts as ppts

#test with third-order
#data
N=3
xtest=num.zeros((3))
ytest=num.zeros((3))
for i in range(0,3):
    xtest[i]=ppts.diff_nodes[3,i]
    ytest[i]=math.sin((xtest[i]+1)/2.0*math.pi/2)
#list: order, left, right, no intervals:
params = [3,-1.0,1.0,1] 
# first derivative
instMatDiff =ppts.MatDiffPoints(params,xtest)
instMatDiff.setallDiffMat()
#display
instMatDiff.__str__()
#node interpolation
instMatLagr = ppts.MatLagrPol(xtest,ytest,N)
instMatLagr.__str__()
print("verif for xv=1..4/8.0 ")
xv=num.zeros((4))
yv=num.zeros((4))
for i in range(0,4):
   xv[i]=-0.75+(i)/2.0
yv=instMatLagr.interp((xv-num.ones((4)) )*2.0)
print(str(xv))
print(str(yv))

#instanciate a linear solver class:
import pendulum as linear
import newmeshcol as newmesh
nc=2
ms=[2,2]
N=2
M=3
X0=num.array([0.157,0.0,0.314,0.0])
trange=0.25
params=[1.2,0.7,0.5,0.3,0.02,0.02,2.5]
aApprox=linear.ApproxRes(nc,ms,N,M)
aLinearEqn=linear.LinearEqnSolver(nc,ms,N,M,X0,trange,params)
mstar=4
n=10
nold=6
nfix=2
xi=num.zeros((10))
xiold=num.linspace(0,0.25,6)
fixpnt=num.array([0.041667,0.20333])
Z=num.zeros((10))
DMZ=num.zeros((4,10))
accum=num.zeros((4))
aNewMesh=newmesh.newmesh(mstar,n,nold,nfix,xi,xiold,Z,DMZ,accum,fixpnt)
 



