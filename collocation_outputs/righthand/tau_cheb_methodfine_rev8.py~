# REV 8 try order 6 instead of 4 (hence size of matrix Acpy=8)
#study of the condition number of the system matrix
#deriving from the tau-collocation method applied to the system
# y''+lamb*y(x) = f
import math as math
#import numpy as num

Mmax=9 #unused now
lamb=4
y_limit=-0.35
### FIRST PART: SUBROUTINES
#derivative in terms of truncated series recursion
def der_order_two(M):
   A=[[0. for j in range(M+2)] for i in range(M+2)]  # WAS M,M
   B2=[[0. for j in range(M)] for i in range(M)]  # WAS M,M
   # building first order
   for i in range(0,M+2): # WAS M-1
      for j in range(0,M+2): # WAS M-1
         if ((j>=i) and ((i+j)%2==0)):
            A[i][j]=2*(j+1) # we begin from b1,b2.. vector
   print("A")
   print(A)
   # multiplying by itself towards right just before truncation:
   for i in range(0,M): # WAS M-2
      for j in range(0,M): # WAS M-1
         term_ij=0.
         for k in range(0,M): #WAS M-1
            term_ij=term_ij+A[i][k]*A[k+1][j]
         B2[i][j]=term_ij
   return B2

# inverse a upper triangular matrix
def inverse_U_mat(A1,M):
   S=[[0 for j in range(M)] for i in range(M)]
   for c in range(0,M):
      E1=[0. for i in range(M)]
      E1[c]=1.    # c-th vector of the canonical basis
      S1=[0. for i in range(M)]
      for k in range(0,M):
         term=E1[M-1-k]  # descending (index from M to one since upper trig)
         for l in range(0,k):
            term=term - A1[M-1-k][M-1-l]*S1[M-1-l]
         if ((math.fabs(term) < 0.001) and (math.fabs(A1[M-1-k][M-1-k]) < 0.001)):
            S1[M-1-k]=0.  # 0/0 : priority over zero than infinity
         else:
            S1[M-1-k]=term/A1[M-1-k][M-1-k]
      #recopy
         S[M-1-k][c]=S1[M-1-k]
   return S
#end inv U

# fred: 06/08/modified version for finding eigenvectors to A (M=A-eigv*Id)
def inverse_LU_pivot(B1,M):
   print("LU procedure")
   L=[[0 for j in range(M)] for i in range(M)]
   U=[[0 for j in range(M)] for i in range(M)]   
   for i in range(0,M):
      L[i][i]=1.0
   for i in range(0,M):
      for j in range(0,M):
         U[i][j]=B1[i][j]     # fred this copy operation was necessary, without ot we would get an erasing of original data
   L1=L
   #print(str(L)) #debug
   #print(str(U))  #debug
   for j in range(0,M-2):  # but cols=1.. M-2 and M-1 pivoting are necessary
      for i in range(M-2,M):  #(M) last two rows have off low diagonal element
         L1=L
         print(str(U[i][j])) #debug
         if (math.fabs(U[i][j])< 0.001):
            print("pivot not necessary: ("+str(i)+","+str(j)+")")
            #print(str(L))
            #print(str(U)) #debug
            #U[i][j]=0.001
            continue
         pivot=U[j][j]/U[i][j]        # pivot j
         if (math.fabs(pivot) < 0.001):
            print("pivot ("+str(j)+") is zero")
            break
         for c in range(0,M): 
            U[i][c]=U[i][c]-1./pivot*U[j][c]
            L[c][j]=L1[c][j]+1./pivot*L1[c][i] # this is already the inverse, right multiplication by inverse of pivot
   #last line, last (-1) column:
   L1=L
   if (math.fabs(U[M-1][M-2])< 0.001):
      print("end pivot not necessary")
   else:
      pivot=U[M-2][M-2]/U[M-1][M-2]
      print("end pivot: "+str(pivot))
      #print(str(U)) #debug
      for c in range(0,M):
         U[M-1][c]=U[M-1][c]-1./pivot*U[M-2][c]
         L[c][M-2]=L1[c][M-2]+1./pivot*L1[c][M-1]
    # since two matrices have to be returned a dictionnary structure is best adapted
   output={}
   output["L"]=L
   output["U"]=U
   print("end LU procedure")
   return output
#end LU

def solve_veigenvector(A1,val,M):
   Mat=[[0.0 for j in range(M)] for i in range(M)]
   for i in range(0,M):
      for j in range(0,M):
         Mat[i][j]=A1[i][j]
 
   Y=[0. for i in range(M)]
   Xv=[0. for i in range(M)]
   for i in range(M):
      Mat[i][i]=Mat[i][i]-val  # correspond to a singular system for eigenvalue val
   print("input veigenvector")
   print(str(Mat))
   decLU=inverse_LU_pivot(Mat,M)
   Lv=decLU["L"]
   Uv=decLU["U"]
   print("Lv:\n"+str(decLU["L"]))
   print("Uv:\n"+str(decLU["U"]))
   # if the multiplicity is only one then U has all zeros on its last line
   # Y, solution of LY=0 is assumed to be zero, we modify it by imposing its last
   # coefficient being 1.0, assuming the eigenspace dimension is one
   # prepare table of index for permutation

#   decLUp=inverse_LU_pivot(Msub,M)
   Lp=Lv
   Up=Uv
#   print(str(Lp))
#   print(str(Up))
   #solving Lp.Y=0 :
   if (math.fabs(Lp[M-1][M-1])<0.001):
      Y[M-1]=1.0
   #solving Up.X=Y :
   for i in range(M-1):
      term_res=Y[M-i-1]  # index descending
      for j in range(M-i-1,M):
         term_res = term_res - Uv[M-i-1][j]*Xv[j]
      Xv[M-i-1]=term_res/Uv[M-i-1][M-i-1]  # 
   print("Xv \n");
   print(str(Xv))
   return Xv
# end solve_veigenvector

def norm_abs(A1,M):
   term=0.
   for i in range(0,M):
      for j in range(0,M):
         term=term+math.fabs(A1[i][j])
   return term

##### PROGRAM PRINCIPAL 
# firstly test inversion of upper triangular
def main_subroutine_LU():
   print("test of LU inverse procedure:")
   #Mtest=[[2.,-0.5,3.,0.],[0,2.,0.5,1.],[1,-1,0.5,-1],[1.,1.,1.,1.]]
   Mtest=[[10.,0.,8.,0.],[0,10.,0.,24.],[1.,-1.,1.,-1.],[1.,1.,1.,1.]] 
   outputLU=inverse_LU_pivot(Mtest,4)
   Ltest=outputLU["L"]
   Utest=outputLU["U"]
   print("L\n")
   print(Ltest)
   print("U\n")
   print(Utest)
   print("M\n")
   print(Mtest)
   print("\n")

# DATA MUST BE COPIED FROM MATLAB HERE:
# FRED: there is a numeric problem for my LU with singular matrice:
# THE ONLY SOLUTION IS TO GO ON FULLY WITH MATLAB FOR EIGVEC SEARCH
   print("import eigenvals from Matlab: ") 
   eigval=[11.5030 + 1.1151*1j, 11.5030 - 1.1151*1j, -1.3229, 0.335];# 0.3169]
   Xv3=solve_veigenvector(Mtest,eigval[3], 4)
   print("found an eigenvector to val: "+str(eigval[3]))
   print(str(Xv3))


# !!!! PRGM fix an order 5 for the collocation tau for diff eqn y"+lamb*y=-lamb*cos^2(pi*x)-2*pi^2*cos(2*pi*x)
M=8

def main_subroutine_der(M):
   print(M)

   Mat=[[0. for j in range(M)] for i in range(M)]

   Mat=der_order_two(M+1) #derivative spectral coeff in function of scalar ones
                          # M+1 because by the way there will just M nnz rows
   print(Mat)
   print("\n")
   Matcpy=[[0. for j in range(M+2)] for i in range(M+2)]
#assembly
   for i in range(0,M):
      for j in range(0,M+1):  # WAS: M+1
         Matcpy[i][j+1]=Mat[i][j]  # c0,c1 for Ti determines the rows, but 1st co is zero, hence j+1
      Matcpy[i][i]=Mat[i][i]+lamb/4. # diagona elmt are for scalar operator
   # the two last lines are interpolation traduction of boundary conditio(-1,1)
   for j in range(0,M+2):
      if (j%2 == 0):
         Matcpy[M][j]=1.
      else:
         Matcpy[M][j]=-1.
      Matcpy[M+1][j]=1.
   print("\n")
   print(Matcpy)
   print(norm_abs(Matcpy,M))
   print("inverse\n")
   return Matcpy
#Minv=inverse_U_mat(Matcpy,M+2)




#Bilan: ce n'est pas possible d'utiliser un solveur direct car Acpy n est pas triangulaire a cause des conditions aux limite (seule l'integration pure le permet)
# il faut passer par Matlab ou numpy pour resoudre le systeme .. et s'apercevoir qu il est tres mal conditionne (cond(Acpy) >  132 !)
#import de Matab:
# solution:
#eigenvalues of A=Matcpy and A^(-1):
#  17.0641 + 5.6007i
#  17.0641 - 5.6007i
#  -8.6524 + 7.5626i
#  -8.6524 - 7.5626i
#  -8.0997 + 0.0000i
#  -2.7236 + 0.0000i
# THIRD PART: FOURIER DECOMPOSITION OF THE RHS TERM: output to Matlab
#analysis of the exact solution with the source=-4*cos^2(pi*x) (always combined 
#with Matlab): a short way to do this instead of by means of integration is to
#compute the Fourier series of the source assuming it can be periodicized ofterm
#cos(2*pi*x), sin(2*pi*x) (these approx will be good at least around zero)
# for the decumposition, the coeff of these two terms are computed by Gauss Legendre integration
# L3=5/3*X^3 -X, valid to compute integral of polynoms up to order 5
print("Third part:")
def calc_Fourier_cospi(K):
   if (K==4): # order of decumposition as multiple from cos(3x)
      xiLG5=[-0.90618, -0.538469, 0, 0.538469, 0.90618]
      wiLG5=[0. for i in range(5)]

      for i in range(5):  #re-compute the Gauss Legendre weights wi
         xi=xiLG5[i]
         LPi1=1/16.*(231.*xi**6-315.*xi**4+105*xi**2-5.)  # 4-th Legendre polynom
         dLPi=1/8.*(5*63.*xi**4-3*70*xi**2+15.)               # derivative of 3-th Legendre polynom
         wiLG5[i]=-2./(5.+1.)/LPi1/dLPi  #formule utile
  #### xiLG5[i]=0.5+1/2.*xiLG5[i]  #scale for integration over [0,1] NO PREFER formal derivation cos(2pi*x) => -cos(pi*z)
      print("Legendre coeff");
      print(str(xiLG5))
      print(str(wiLG5))
      print("\n")
      LGcst=[0. for i in range(5)]
      LGcos=[0. for i in range(5)]
      LGsin=[0. for i in range(5)]
      LGcos32=[0. for i in range(5)]
      LGsin32=[0. for i in range(5)]
      LGcos3=[0. for i in range(5)]
      LGsin3=[0. for i in range(5)]
      LGcos6=[0. for i in range(5)]
      LGsin6=[0. for i in range(5)]
      pival=math.pi
      for i in range(5):  # careful, the interval of integration is [0,1] not [-1,1] !
         xi=math.pi/pival*xiLG5[i] #scale
         LGcst[i]=math.cos(pival*xi)*wiLG5[i]
         LGcos[i]=math.cos(pival*xi)*math.cos(2*xi)*wiLG5[i] # try although not correct
         LGsin[i]=math.cos(pival*xi)*math.sin(2*xi)*wiLG5[i]   
         LGcos32[i]=math.cos(pival*xi)*math.cos(3*xi)*wiLG5[i]# the fundamental
         LGsin32[i]=math.cos(pival*xi)*math.sin(3*xi)*wiLG5[i]
         LGcos3[i]=math.cos(pival*xi)*math.cos(6*xi)*wiLG5[i] # twice the fundam
         LGsin3[i]=math.cos(pival*xi)*math.sin(6*xi)*wiLG5[i]   
         LGcos6[i]=math.cos(pival*xi)*math.cos(12*xi)*wiLG5[i] #fourth harmonic
         LGsin6[i]=math.cos(pival*xi)*math.sin(12*xi)*wiLG5[i]

      decump=[0. for i in range(9)]
      decump[0]=1.*sum(LGcst)
      decump[1]=1.*sum(LGcos)
      decump[2]=1.*sum(LGsin)
      decump[3]=1.*sum(LGcos32)
      decump[4]=1.*sum(LGsin32)
      decump[5]=1.*sum(LGcos3)
      decump[6]=1.*sum(LGsin3)
      decump[7]=1.*sum(LGcos6)
      decump[8]=1.*sum(LGsin6)
   return decump
   # end of function calc_Fourier_cospi

decump=calc_Fourier_cospi(4)  # bit strange but ok: cos2,cos3,cos6,cos12
c0=decump[0]
c1=decump[1]
d1=decump[2]
c32=decump[3]
d32=decump[4]
c3=decump[5]
d3=decump[6]
c6=decump[7]
d6=decump[8]
# (c3/2) because c3=ps(cos(pi*x), cos3x) must be divided by length of interval (-1,1)
print("Fourier coeffs: "+str(c0)+" : "+str(c1)+" : "+str(d1)+" : "+str(c32)+" : "+str(d32)+" : "+str(c3)+" : "+str(d3)+" : "+str(c6)+" : "+str(d6))
taun=0.09

# inserted before the chebyshev decomposition: automatic-decumposition from Taylor expansion of cos(6x):
def expand_cos6x():
  # in waiting to complete, take constant Ucos3 and Usin3
    return 0

# the following comes from the inversion x^i <-> Ti(x) of cos(3x), sin(3x) in Matlab
Ucos5=[   -0.3008,  0,  -1.0371, 0,  0.2320, 0]
Usin5=[   0, 0.8906, 0, -0.4922,  0, 0.1266]

#below: sign (-1) because the change of variable gives cos(2*pi*x)~= cos(pi*z+pi)= -cos(pi*z)
Uvar1 = [0.0 for i in range(M)]
Uvar1[0] =  (-(lamb/8.))*(-c32)*Ucos5[0]  # %(-(lamb/8.+pi^2/40))*(c3/2)*Ucos6(1:2:5) ;  /8 =1/(2^2)/2
Uvar1[0] = Uvar1[0] -(lamb/8)*(-c0) -lamb/8 #%Uvar1(1) -(lamb/8+pi^2/40)*(c0/2) -lamb/8 
Uvar1[1] = (-(lamb/8.))*(-d32)*Usin5[1] 
Uvar1[2] =  (-(lamb/8.))*(-c32)*Ucos5[2]
Uvar1[3] = (-(lamb/8.))*(-c32)*Usin5[3]
Uvar1[4] = taun
Uvar1[5] = 0.
Uvar1[6] = y_limit # % see true solution by Fourier series below to see that u(1) ~= -0.25
#export in screen
print("Uvar1")
print(str(Uvar1))
#then a secondary, perturbed tight hand side: Uvar2=Uvar1+delta OR just the perturbation
CONTROL_sel_pb="disturb"
Uvar2 = [0.0 for i in range(M)]
if (CONTROL_sel_pb=="disturb"):
    Uvar2[0] = (math.pi**2/40)*(c32)*Ucos5[0]
    Uvar2[1] = 0.
    Uvar2[2] = (math.pi**2/40)*(c32)*Ucos5[2]
    Uvar2[3] = 0.
    Uvar2[4] = taun
    Uvar2[5] = 0.
    Uvar2[6] = y_limit
elif (CONTROL_sel_pb=="sumcstdis"):
    Uvar2[0] = Uvar1[0] + (math.pi**2/40)*(c32)*Ucos5[0]
    Uvar2[1] = Uvar1[1]
    Uvar2[2] = Uvar1[2] + (math.pi**2/40)*(c32)*Ucos5[2]
    Uvar2[3] = Uvar1[3]
    Uvar2[4] = Uvar1[4]
    Uvar2[5] = Uvar1[5]
    Uvar2[6] = Uvar1[6]
#export in screen
print("Uvar2")
print(str(Uvar2))


# ************** then SOLVING ************* by inversion of numpy

#def main_chebsolver(M, Uvar):
#    print("main cheb solver")
#    #transform matrix MatCpy in numpy array 2x2
#    Anum = num.array(main_subroutine_der(M-2))
#    Anum.reshape(M,M)
    #Anum=Anum[0:M,0:M]  #truncate see the remark in c0,c1 but first column is zero in main_subroutine_der
#    print(str(Anum))
#    print("condition number Anum: "+str(num.linalg.cond(Anum,2)))
    #idem for the rhs, which is still a global variable
#    Uvarnum = num.array(Uvar)
        #and complete with the boundary condition
#    Uvarnum[M-2] = 0.
#    Uvarnum[M-1]  = y_limit
#    Muk1 = num.linalg.solve(Anum, Uvarnum)
#    print("solved coeffs")
#    print(str(Muk1))
#    return Muk1

# CALL SYSTEM 1
#Mukvar1 = main_chebsolver(M, Uvar1)

# THIRD PART BIS: CHEBYSHEV DECOMPOSITION OF THE RHS TERM: output to Matlab
#analysis of the exact solution with the source=-4*cos^2(pi*x) 
# it has been inquired that instead of doing the Fourier decumposition, the
# projection on the Tn(x) by chebyshev weight scalar product leads directly to
# the RHS term
print("Third part bis")
def calc_chebyshev_cospi(Tz,za,nz):
   xiCG5=[-math.cos(math.pi/10), -math.cos(3*math.pi/10), 0., math.cos(3*math.pi/10), math.cos(math.pi/10)]
   wiCG5=[math.pi/5 for i in range(5)] # for Chebyshev this is much simpler!

   print("\n Chebyshev Gauss coeff");
   print(str(xiCG5))
   print(str(wiCG5))
   print("\n")
   CGcst=[0. for i in range(5)]
   CGcos=[0. for i in range(5)]
   CGcos2=[0. for i in range(5)]
   CGcos3=[0. for i in range(5)]
   CGcos4=[0. for i in range(5)]
   CGcos5=[0. for i in range(5)]
   CGcos6=[0. for i in range(5)]
   CGcos7=[0. for i in range(5)]
   pival=math.pi
   for i in range(5):  # Chebyshev weighted scalar product !
      #yi=0.5+1/2.*xiCG5[i]  #scale for integration over [0,1] NO
      xi=xiCG5[i] #scale YES but now we need to approximate
                 # .. over subintervals of [0,T], T being the period of homogeneous true solution or at least a guess of it, that is [n*za*T, (n+1)*za*T]
      vi=(xi+1)*((za*Tz)/(2*pival))+nz*za*Tz/pival
      CGcst[i]=math.cos(pival*vi)*1*wiCG5[i]
      CGcos[i]=math.cos(pival*vi)*(xi)*wiCG5[i]
      CGcos2[i]=math.cos(pival*vi)*(2*xi**2-1)*wiCG5[i]
      CGcos3[i]=math.cos(pival*vi)*(4*xi**3-3*xi)*wiCG5[i]
      CGcos4[i]=math.cos(pival*vi)*(8*xi**4-8*xi**2+1)*wiCG5[i]  
      CGcos5[i]=math.cos(pival*vi)*(16*xi**5-20*xi**3+5*xi)*wiCG5[i]   
      CGcos6[i]=math.cos(pival*vi)*(32*xi**6-48*xi**4+18*xi**2-1)*wiCG5[i]
      CGcos7[i]=math.cos(pival*vi)*(64*xi**7-112*xi**5+56*xi**3-7*xi)*wiCG5[i]
   #au dela du degre 6 l integration de Gauss nest plus exacte:
   decump=[0. for i in range(8)]

   decump[0]=1/math.pi*sum(CGcst) # needed divide by 2 since scalar product of T0=1 is pi
   decump[1]=2/(math.pi)*sum(CGcos)
   decump[2]=2/(math.pi)*sum(CGcos2)
   decump[3]=2/(math.pi)*sum(CGcos3)
   decump[4]=2/(math.pi)*sum(CGcos4)
   decump[5]=2/(math.pi)*sum(CGcos5)
   decump[6]=2/(math.pi)*sum(CGcos6)
   decump[7]=2/(math.pi)*sum(CGcos7)

   return decump 
# end calc_chebyshev_cospi
# then write a write_chebvisu subroutine
def write_chebvisu(decump, Np):
   table=[0. for i in range(Np)]
   for i in range(Np):
      ti = -1. + 2./(Np-1)*i
      si=0.
      si=si+decump[0]*1
      si=si+decump[1]*ti;
      si=si+decump[2]*(2*ti**2-1)
      si=si+decump[3]*(4*ti**3-3*ti)
      si=si+decump[4]*(8*ti**4-8*ti**2+1)
      si=si+decump[5]*(16*ti**5-20*ti**3+5*ti)
      si=si+decump[6]*(32*ti**6-48*ti**4+18*ti**2-1)
      table[i]=si
   return table

# now use this function calc_chebyshev_cospi firstly on the whole interval
Tg=math.pi
za=1.0
ng=0
decump=calc_chebyshev_cospi(Tg,za,ng)
a0=decump[0]
a1=decump[1]
a2=decump[2]
a3=decump[3]
a4=decump[4]
a5=decump[5]
a6=decump[6]
a7=decump[7]      
print("Chebyshev coeffs: "+str(a0)+" , "+str(a1)+" , "+str(a2)+" , "+str(a3)+" , "+str(a4)+" , "+str(a5)+" , "+str(a6)+ ", "+str(a7))

# FOURTH PART
#right hand side: computed manually on the poynom basis of chebyshev
Uvard1=[0. for i in range(M)]
taun=0.09;

Uvard1[0]=(-(lamb/8.))*(a0)-lamb/8 ; #*(0.8644) 2.44; %a0, 1/4 = 1/2^2 for scaling and 1/2 for cos^2
Uvard1[1]=0.0 # a1=0.
Uvard1[2]=(-(lamb/8.))*(a2) # a2=-0.5
Uvard1[3]=(-(lamb/8.))*(a3) # 
Uvard1[4]=taun # a3=0.25 overriden by taun
Uvard1[5]=0. # a4=0 overriden by boundary 1
Uvard1[6]=y_limit # a5=0 overriden by boundary 2
#export to screen
print(str(Uvard1))

# CALL SYSTEM d1
#Mukvard1 = main_chebsolver(M, Uvard1)

# FIFTH PART: perform the approximation over three sub-intervals
#  [0, Tg/4], [Tg/4, Tg/2] and [Tg/2, Tg*3/4]
# first interval:
Tg=math.pi
za=.5
ng=0
decump=calc_chebyshev_cospi(Tg,za,ng)
a0=decump[0]
a1=decump[1]
a2=decump[2]
a3=decump[3]
a4=decump[4]
a5=decump[5]
a6=decump[6]
a7=decump[7]
      
print("first interval Chebyshev coeffs: "+str(a0)+" , "+str(a1)+" , "+str(a2)+" , "+str(a3)+" , "+str(a4)+" , "+str(a5)+" , "+str(a6)+ ", "+str(a7))
# the building and solving of the linear coll problem should be done in a dedicated C++ program that enables intervals-solving and before, reading of rhs coeff:

table1=write_chebvisu(decump, 10)

# second interval:
Tg=math.pi
za=.5
ng=1
decump=calc_chebyshev_cospi(Tg,za,ng)
a0=decump[0]
a1=decump[1]
a2=decump[2]
a3=decump[3]
a4=decump[4]
a5=decump[5]
a6=decump[6]
a7=decump[7]
      
print("second interval Chebyshev coeffs: "+str(a0)+" , "+str(a1)+" , "+str(a2)+" , "+str(a3)+" , "+str(a4)+" , "+str(a5)+" , "+str(a6)+ ", "+str(a7))

table2=write_chebvisu(decump, 10)

# third interval:
Tg=math.pi
za=.5
ng=2
decump=calc_chebyshev_cospi(Tg,za,ng)
a0=decump[0]
a1=decump[1]
a2=decump[2]
a3=decump[3]
a4=decump[4]
a5=decump[5]
a6=decump[6]
a7=decump[7]
      
print("third interval Chebyshev coeffs: "+str(a0)+" , "+str(a1)+" , "+str(a2)+" , "+str(a3)+" , "+str(a4)+" , "+str(a5)+" , "+str(a6)+ ", "+str(a7))

table3=write_chebvisu(decump, 10)

tabfile=open("TheRhs", "w")
for i in range(10):
   tabfile.write(str(table1[i])+"\n")
for i in range(10):
   tabfile.write(str(table2[i])+"\n")
for i in range(10):
   tabfile.write(str(table3[i])+"\n")
tabfile.close()



