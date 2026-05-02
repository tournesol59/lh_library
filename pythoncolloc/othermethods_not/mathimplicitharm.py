import numpy as num
import math as math

def freq_func(xvec):
    yvec=num.zeros((length(xvec))
    for i in range(0,length(xvec)):
       yvec[i]=math.exp(-xvec[i]**2)

def guess_func(xvec):
    yvec=num.zeros((length(xvec))
    for i in range(0,length(vec)):
       yvec[i]=cos(xvec[i])
    return yvec

class implicit:

    def __init__(self, N, M, deltah):
        self.__N=N
        self.__M=M
        self.__varx=num.zeros((N))
        for i in range(0,N):
            self.__varx[i]=i*deltah
        self.__varxij=num.zeros(((N-1)*M+1))
        self.__var_prev=num.zeros(((N-1)*M+1))
        self.__var_fwd=num.zeros(((N-1)*M+1))
        self.__deltah=deltah

    def advance():
        M=self.__M
        N=self.__N
        deltahij=self.__deltah/(M-1)

        Amat=num.zeros(((N-2)*(M-1)+2*M,(N-2)*(M-1)+2*M))
        bvec=num.zeros(((N-2)*(M-1)+2*M))
        qxi=freq_func(self.__xvarx)
        # first block of Amat
        Amat[0,0]=2/deltahij**2
        Amat[0,1]=-1/deltahij**2
        bvec[0]=(-1/deltahij**2-qwi[0]) * (self.__var_prev[0])
        for l in range(1,M-1):
           Amat[l,l-1]=-1/deltahij**2
           Amat[l,l]=2/deltahij**2
           Amat[l,l+1]=-1/deltahij**2
           bvec[l]=(-qwi[0]) * (self.__var_prev[l])
        Amat[M-1,M-2]=-1/deltahij**2
        Amat[M-1,M-1]=2/deltahij**2
        bvec[M-1]=(-1/deltahij**2-qwi[0]) * (self.__var_prev[M-1])

        #middle blocks of size (M-1)
        for i in range(0,N-2): # loop over the no of intervals
        # do not forget: the points at the extremities of intervals belongs to rhs as they are guessed
            Amat[i*(M-1)+M,i*(M-1)+M]=2/deltahij**2
            Amat[i*(M-1)+M,i*(M-1)+M+1]=-1/deltahij**2
            bvec[i*(M-1)+M]=(-1/deltahij**2-qwi[i+1]) * (self.__var_prev[i*(M-1)+M])      
            for l in range(1,M-1):
               Amat[i*(M-1)+M+l,i*(M-1)+M+l-1]=-1/deltahij**2
               Amat[i*(M-1)+M+l,i*(M-1)+M+l]=2/deltahij**2
               Amat[i*(M-1)+M+l,i*(M-1)+M+l+1]=2/deltahij**2
               bvec[i*(M-1)+M+l]=(-qwi[i+1]) * (self.__var_prev[i*(M-1)+M+l)
            Amat[(i+1)*(M-1)+M,(i+1)*(M-1)+M]=-1/deltahij**2







            
            Amat[i*M+M-1,i*M+M-1]=2/deltahij**2
            # values at xvec[i*M]
            # vector b is filled with terms in order zero of previous interval
            bvec[i*M]=-qxi[i]*self.__var_prev[i*M]
        #last rows for A:
        Amat[(N-1)*M+0,(N-1)*M+0]=2/deltahij**2
        Amat[(N-1)*M+0,(N-1)*M+1]=-1/deltahij**2        
        for l in range(1,M):
           Amat[(N-1)*M+l,(N-1)*M+l-1]=-1/deltahij**2
           Amat[i*M+l,i*M+l]=2/deltahij**2

        bvec[M*(N-1)+1]=-qwi[N]*self.__var_prec[i*M]
        
        self.__var_fwd=num.matmul
    # end class implicit
