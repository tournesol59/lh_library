import numpy as num
import volteralegend.packcollocpts as cpts
#os.chdir("C:\\MinGW\\msys\\1.0\home\\lh_library\\pythoncolloc")
#import matplotlib.pylab as plt

poles=['-1','1']
#parameters associated to different changes of vars, thus associated to different poles
changevar=[{'B_2': 0.9, 'cmin': 0.05}, {'B_2': 0.95, 'cmin': 1.0}, {'B_2': 0.95, 'cmin': 1.0}]
#test third-order of a linear ode1 over subinterval
#cmin=0.1
#B_2=0.9 # shrink time to avoid inf singularity of the poles functions, applies to poslnfunc
y_init=1.0
N=3

#define scalar functions a(x) of the variable x vs time for a(x)y'+y=b(x)*q(x)
# diff eqn w/o transform: y'(x)*(lnfunc(x))+y(x)=b(x)  over [-1,1]
def lnfunc(x):
    if (int(poles[0])==-1) and (int(poles[1])==1):
        return (1-x**2)
    elif (int(poles[0])==-1) and (int(poles[1])==0):
        return x*(1+x)
    elif (int(poles[0])==0) and (int(poles[1])==1):
        return x*(1-x)
 
# diff eqn with transform: t<-c(1-x)/(1-Bx) 
#     gives y'(x)*(poslnfunc(x))+y(x)=b*q(x) with
#     poslnfunc(x)=...
#     posqfunc(x)=...
def poslnfunc(x):
    if (int(poles[0])==-1) and (int(poles[1])==1):
        B_2=(changevar[0])['B_2']
        cmin=(changevar[0])['cmin']
        return 1/cmin/2.*((1-B_2*x)**2-cmin**2*(1+x)**2)/((1+B_2))
    elif (int(poles[0])==-1) and (int(poles[1])==0):
   #     B_2=(changevar[1])['B_2']
   #     cmin=(changevar[1])['cmin']
        return (1+x)
    elif (int(poles[0])==0) and (int(poles[1])==1):
   #     B_2=(changevar[2])['B_2']
   #     cmin=(changevar[2])['cmin']
        return (1-x)

def poschfunc(x):
    # change of variable
    if (int(poles[0])==-1) and (int(poles[1])==1):
        B_2=(changevar[0])['B_2']
        cmin=(changevar[0])['cmin']
        return (-x)/(1-B_2+B_2*x)        
      #  return cmin*(1+x)/(1-B_2*x)
    elif (int(poles[0])==-1) and (int(poles[1])==0):
        B_2=(changevar[1])['B_2']
        cmin=(changevar[1])['cmin']
        return (1-B_2)/(1+B_2*x)
    elif (int(poles[0])==0) and (int(poles[1])==1):
        B_2=(changevar[2])['B_2']
        cmin=(changevar[2])['cmin']
        return (1-B_2)/(1-B_2*x)

def posqfunc(x):  #q(x)
    if (int(poles[0])==-1) and (int(poles[1])==1):
        B_2=(changevar[0])['B_2']
        cmin=(changevar[0])['cmin']
        return 1/2.*((1-B_2*x)**2-cmin**2*(1+x)**2)/((1-B_2*x)**2)
    elif (int(poles[0])==-1) and (int(poles[1])==0):
        B_2=(changevar[1])['B_2']
        cmin=(changevar[1])['cmin']
        return (1-B_2)*B_2*(1-x)/(1-B_2*x)**2    # tbc
    elif (int(poles[0])==0) and (int(poles[1])==1):
        B_2=(changevar[2])['B_2']
        cmin=(changevar[2])['cmin']
        return (1-B_2)*B_2*(1+x)/(1+B_2*x)**2 # tbc

def exactsol(x):
    # three exact solutions depending on the poles defining the diff equn
    if (int(poles[0])==-1) and (int(poles[1])==1):
        # solution on interval [-1,1]
        return y_init*(1-x)/(1+x)
    elif (int(poles[0])==-1) and (int(poles[1])==0):
        # solution on interval [-1,0]
        return y_init*(1+x)/x
    elif (int(poles[0])==0) and (int(poles[1])==1):
        # solution on interval [0,1]
        return y_init*(1-x)/x

def exactchsol(z):
    # three exact solutions depending on pole, this time with change of var on [-1,1]
    if (int(poles[0])==-1) and (int(poles[1])==1):
        # solution on interval [-1,1]
        B_2=(changevar[0])['B_2']
        cmin=(changevar[0])['cmin']
        x=cmin*(1+z)/(1-B_2*z)
        return y_init*(1-x)/(1+x)
    elif (int(poles[0])==-1) and (int(poles[1])==0):
        # solution on interval [-1,0]
        B_2=(changevar[1])['B_2']
        cmin=(changevar[1])['cmin']
        x=(1-B_2)/(1+B_2*z)
        return y_init*(1+x)/x
    elif (int(poles[0])==0) and (int(poles[1])==1):
        # solution on interval [0,1]
        B_2=(changevar[2])['B_2']
        cmin=(changevar[2])['cmin']
        x=(1-B_2)/(1-B_2*z)
        return y_init*(1-x)/x

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
# recatch the data structures from instMatDiff:
xtest_ext=instMatDiffP.get_xext()
Kernel=instMatDiffP.get_KM()
factors=instMatDiffP.get_factors()
print("points:");
print(str(xtest_ext))
#print("Kernel")
#print(str(Kernel))
print("factors:")
print(str(factors))
print("diff matrix:")
instMatDiffP.__str__()

Amat=num.identity(N+2);
# perform derivative matrix computation
D=instMatDiffP.getDM1()
b=num.zeros((N+2))
for i in range(0,N+2):
   D[i,:]=D[i,:]*poslnfunc(xtest_ext[i])
   b[i]=b[i]*posqfunc(xtest_ext[i]) #if not zeros
print(str(D))
Amat=Amat+D
a1=Amat[1:N+2,0]  # first column as opposed to Matlab version (last column)
Amat=Amat[1:N+2,1:N+2]
b = b[1:N+2] - y_init*a1
print("with identity sum the system matrix should be")
print(str(Amat))
print("and the right hand side")
print(str(b))
# SOLVING THE SYSTEM
yvec=num.matmul(num.linalg.inv(Amat), b)

print("And the solved values should be :")
print(str(yvec))
M=4
# compute polynoms expressions
def polyint(yvector,M,N):
    Np=(N+1)*(M) #+1
    xint=num.zeros((Np))
    yint=num.zeros((Np))
    for i in range(0,N+1):
        for j in range(0,M):
            xij=(xtest_ext[i+1]-xtest_ext[i])/(M-1)*j
            xint[i*M+j]=xtest_ext[i]+xij
            yij=0.0
            table=instMatDiffP.get_LX(xij)  #instMatDiffP.instKernel.fillVarVector(xij)
            for l in range(0,N+2):  # skip l=0 for yinit
# computes P_l(xij)
                prod_l=1.0
                for k in range(0,l):  # corrected was: l-1, is: l
                    prod_l=prod_l*table[k]
                for k in range(l+1,N+2):
                    prod_l=prod_l*table[k]
                prod_l=prod_l/factors[l]  # divide by times,k (Xl-Xk) k!=l
                if (l==0):      # skip l==0 for y_init
                    yij=yij+y_init*prod_l
                else:
                    yij=yij+yvec[l-1]*prod_l
            yint[i*M+j]=yij
           # print(str(yij))  #debug
    result=num.zeros((Np,2))
    result[:,0]=xint
    result[:,1]=yint
    return result

print("And the interpolated values should be :")
yplot=polyint(yvec,M,N)
print(str(yplot[:,1]))


#plot the functions

xtest=num.zeros(51)
for i in range(0,51):
    xtest[i]=0.+i*(1./50.)
ytest=exactsol(xtest)
#plt.plot(xtest, ytest)

ztest=num.zeros(51)
for i in range(0,51):
    ztest[i]=-1.+i*(2./50.)
ychtest=exactchsol(ztest)
#plt.plot(ztest, ychtest)
#plt.plot(yplot[:,0], yplot[:,1])

#plt.show()
