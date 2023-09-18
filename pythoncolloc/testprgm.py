import numpy as num
import packcollocpts as cpts

#test third-order of a linear ode1
B_1=0.9  #shrink time to avoid inf singularity of 1/(1-x**2) at -1,1
eps_2=0.1
B_2=0.9 # same signficance as B_1 but applies to poslnfunc
y_init=1.0
N=3

#define scalar functions a(x) of the variable x vs time for a(x)y'+y=b(x)*q(x)
# diff eqn w/o transform: y'(x)*(lnfunc(x))+y(x)=b(x)  over [-1,1]
def lnfunc(x):  
   return (1-x**2)
   
# diff eqn with transform: t<-c(1-x)/(1-Bx) 
#     gives y'(x)*(poslnfunc(x))+y(x)=b*q(x) with
#     poslnfunc(x)=...
#     posqfunc(x)=...
def poslnfunc(x):
   cmin=(1-B_2)/2-eps_2
   return 1/cmin*((1-B_2*x)**2-cmin**2*(1+x)**2)/((1+B_2)**2)
def posqfunc(x):  #q(x)
   cmin=(1-B_2)/2-eps_2
   return 1/2.0*((1-B_2*x)**2-cmin**2*(1+x)**2)/((1-B_2*x)**2)

# option 1: take package pre-computed Legendre nodes   
xtest=num.zeros((N))
for i in range(0,N):
   xtest[i]=cpts.diff_nodes[N,i]
#list
params = [N,-1.0,1.0,'norm',1] 
# first derivative
instMatDiffP = cpts.MatDiffPoints(params,xtest)
instMatDiffP.calcKernelFac()
instMatDiffP.setallDiffMat()
instMatDiffP.__str__()
# recatch the data structures from instMatDiff:
xtest_ext=instMatDiffP.get_xext()
factors=instMatDiffP.get_factors()
print(str(factors))
Amat=num.identity(N+2);
# perform derivative matrix computation
D=instMatDiffP.getDM1()
b=num.zeros((N+2))
for i in range(0,N+2):
   D[i,:]=D[i,:]*poslnfunc(xtest_ext[i])
   b[i]=b[i]*posqfunc(xtest_ext[i]) #if not zeros
Amat=Amat+D
a1=Amat[:,1]  # first column as opposed to Matlab version (last column)
# SOLVING THE SYSTEM
yvec=num.matmul(num.linalg.inv(Amat), b-y_init*a1)

print("And the solved values should be :")
print(str(yvec))
M=4
# compute polynoms expressions
def polyint(yvector,M,N):
    Np=(N+1)*(M)+1
    xint=num.zeros((Np))
    yint=num.zeros((Np))
    for i in range(0,N+1):
        for j in range(0,M):
            xij=(xtest_ext[i+1]-xtest_ext[i])/(M-1)*j
            xint[i*M+j]=xij
            yij=0.0
            table=instMatDiffP.get_LX(xij)  #instMatDiffP.instKernel.fillVarVector(xij)
            for l in range(0,N+2):
# computes P_l(xij)
                prod_l=1.0
                for k in range(0,l-1):
                    prod_l=prod_l*table[k]
                for k in range(l+1,N+2):
                    prod_l=prod_l*table[k]
                prod_l=prod_l/factors[l]  # divide by times,k (Xl-Xk) k!=l
                yij=yij+yvec[l]*prod_l
                print(str(yij))
    result=num.zeros((Np,2))
    result[:,0]=xint
    result[:,1]=yint
    return result

print("And the interpolated values should be :")
yplot=polyint(yvec,M,N)
#print(str(yplot[1,:]))
