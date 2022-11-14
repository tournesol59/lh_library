#package of mesh-dependent derivatives of (Legendre) polynoms
# cp Frederic Peugny, LTS
import numpy as num
#from matplotlib import pyplot as plt

diff_nodes=num.zeros((8,7))
#high precision values
diff_nodes[0,:]=[65535] #code infinite value
diff_nodes[1,:]=[0.,0.,0.,0.,0.,0.,0.]
diff_nodes[2,:]=[-0.57735,0.57735,0.,0.,0.,0.,0.]
diff_nodes[3,:]=[-0.7746, 0., 0.7746, 0.,0.,0.,0.]
diff_nodes[4,:]=[-0.86114,-0.340,0.340,0.86114, 0.,0.,0.]
diff_nodes[5,:]=[-0.90618,-0.53857,0.,0.53857, 0.90618, 0.,0.]
#low precision values
diff_nodes[6,:]=[-0.92,-0.65,-0.24,0.24,0.65,0.92, 0.]
diff_nodes[7,:]=[-0.945,-0.74,-0.39,0.,0.39,0.74,0.945]

manual_Btable2=num.zeros((4,4)) #for diff_nodes[2,:]
manual_Btable2[0,:]=[0.11706,0.6856,-1.04645,-0.21328]
manual_Btable2[0,:]=[-0.075392,-1.0464,0.6856,0.17162]
manual_Btable2[0,:]=[-0.0625,0.17524,-0.82476,0.0625]
manual_Btable2[0,:]=[-0.0625,0.17524,-0.82476,0.0625]


class MatLagrKernel(object):
    def __init__(self,xext,nx):
            self.__xext=num.zeros((nx))
            self.__KM=num.zeros((nx,nx))
            self.__LX=num.zeros((nx))	
            self.__factors=num.zeros((nx))				
            self.__nx=nx
            for i in range(0,self.__nx):
                self.__xext[i]=xext[i]

    def fillKernel(self):
        for i in range(0,self.__nx):
            for j in range(0,i):
                self.__KM[i,j]=self.__xext[i]-self.__xext[j]
            self.__KM[i,i]=1e-6 #but care should be taken not to use it
            for j in range(i+1,self.__nx):
                self.__KM[i,j]=self.__xext[i]-self.__xext[j]

    def fillKernelAdapted(self):
    # will be used by MatPseudoPoints to integrate over [-1,r], r <1
    # and only by that class !
        k=int ((self.__nx)/2+1)
        # with rescaled evaluation points diff_nodes
        for l in range(0,self.__nx):
            scal=(self.__xext[l]+1.0)
            for j in range(0,self.__nx):
            # correction: original line: self.__KM[l,j]=(-1.+(diff_nodes[k,l]+1.)*scal)-self.__xext[j]
                self.__KM[l,j]=(-1.+(diff_nodes[k,l]+1.)*scal)-self.__xext[j]

    def calcKernelFac(self): 
	# must be called AFTER a Kernel has been filled
        for l in range(0, self.__nx):
            self.__factors[l]=1.0
            for i in range(0, l-1):
                self.__factors[l]=self.__factors[l]*self.__KM[l,i]
            for i in range(l+1, self.__nx):
                self.__factors[l]=self.__factors[l]*self.__KM[l,i]

    def fillVarVector(self, xv):
	# return a table (xv-xext(i))
        for i in range(0,self.__nx): 
             self.__LX[i]=(xv-self.__xext[i])    
        return self.__LX  # test this!

    def get_KM(self):
        return self.__KM  # test this!
		
    def get_factors(self):
        return self.__factors  # test this!
		
    def __str__(self):
        for i in range(0,self.__nx):
            print(str(i)+"\n")
            for j in range(0,self.__nx):
                print(self.__KM[i,j]);
            print("\n")

#end class MatLagrKernel

class MatDiffPoints(object):

    def __init__(self,param,x):	
    # x: numpy array of points (the mesh)
    # param(0): (list) no. of collocation points per subinterval w/o extrems
    # param(1): left value
    # param(2): right value
    # param(3): choice of nodes: 'norm'=in [-1,1] 'scal'=[p1,p2]
    # param(4) no. of subintervals in the initial mesh.
        self.__xext=num.zeros((param[0]+2))
        self.__xext[0]=param[1] # -1.0 default left value
        self.__xext[param[0]+1]=param[2] # 1.0 default right value
        nx=param[0]+2
        self.__nx=nx
        if (param[3]=='norm'):
            for i in range(1,nx-1):	
                self.__xext[i]=diff_nodes[nx,i]
        elif (param[3]=='scal'):
            for i in range(1,nx-1):
            # formula to test: scale of nodes in interval of span 2.0 to another interval
               self.__xext[i]=self.__xext[0]+(diff_nodes[nx,j]+1.0)/2.0*(self.__xext[nx-1]-self.__xext[0])

        self.__DM = num.zeros((2,self.__nx, self.__nx)) # nx Lagr Points leads to nx polynoms and nx eval pts
        self.__instKernel = MatLagrKernel(self.__xext,nx)

    def calcKernelFac(self):
        self.__instKernel.calcKernelFac()

    def get_factors(self):
        self.__instKernel.calcKernelFac()
        return self.__instKernel.get_factors()

    def setallDiffMat(self):    
    # here implement the lagrange difference and products
    # to have DM(i,j)=deriv(Lagr{j}, dt)(ti)
    #first initialize a Kernel of ALL possible differences combination of the x
        nxe=self.__nx
        # create a class of Lagrange points difference of size nxe=total no.pt
        #instKernel = MatLagrKernel(self.__xext,nxe)
        self.__instKernel.fillKernel()        
        Kernel = self.__instKernel.get_KM() # unclear of whether it is a copy

   #second run a double loop for DM(0,i,j)
        for i in range(0,nxe):
         # set a common factor for D(1)(i,j):
            alphac_i = 1.0
            for l in range(0,nxe):
                if (l != i): # not that alpha_c means completed by (ti-tj)
                    alphac_i = alphac_i*Kernel[i,l]
            for j in range(0,nxe):
                if (j != i):
         # now we have to compute alpha_j to complete DM^(1)(i,j)=al_j/alc_i
                    alpha_j = 1.0
                    for l in range(0,nxe):
                        if (l != j) and (l != i): 
         # this test is necessary to have not simplification by (ti-tj)
                            alpha_j = alpha_j*Kernel[j,l]
                    self.__DM[0,i,j] = alpha_j / alphac_i
                else:
                    #if (j==i):
         # next we have to complete DM^(1)(i,i)=sum 1/(ti-tl)
                    term_ii = 0.0
                    for l in range(0,nxe):
                        if (l != i):
         # this test is necessary to not divide by zero but Kernel(i,l) could be safe
                            term_ii = term_ii + 1.0/Kernel[i,l]
                    self.__DM[0,i,i] = term_ii

# now we have to compute alpha_jp to complete DM"(2)(i,j)=sum,j al_jp/alc_i
  
            for j in range(0,nxe):
               if (j != i):
                  term2_ij = 0.0
                  for p in range(0,nxe):
                      alpha_jp = 1.0
                      for l in range(0,nxe):
                          if (l != j) and (l != i) and (l !=p):
    # this test is necessary to have proper derivation and also not simplification
                             alpha_jp = alpha_jp*Kernel[p,l]
                      term2_ij = term2_ij + alpha_jp
                      self.__DM[1,i,j] = term2_ij/alphac_i
 #end method setallDM() of class MatDiffPoints.

    def __str__(self):
        nxe=self.__nx
        for i in range(0,nxe):
            print(str(i)+"    ")
            for j in range(0,nxe):
                print(self.__DM[0,i,j]);
            print("\n")

    def get_xext(self):
        return self.__xext

    def getDM1(self):
        return self.__DM[0,:,:]

    def get_LX(self,xv):
        return self.__instKernel.fillVarVector(xv)  # test this!		
# end class MatDiffPoints


class MatPseudoPoints(object):
    def __init__(self,param,x):
    # x: numpy array of points (the mesh)
    # param(0): (list) no. of collocation points per subinterval w/o extrems
    # param(1): left value
    # param(2): right value
    # param(3): no. of subintervals in the initial mesh
    # param(...) think about other params..
        self.__xext=num.zeros((param[0]+2))
        self.__xext[0]=-1.0 # default left value
        for i in range(1,param[0]+1):
            self.__xext[i]=x[i-1]
        self.__xext[param[0]+1]=1.0 # default right value
        self.__nx=param[0]+2
        self.__Btable = num.zeros((self.__nx, self.__nx)) # nx Lagr Points

    def setButcherMat(self):    
        nxe=self.__nx
        # create a class of Lagrange points difference of size nxe=total no.pt
        instKernel = MatLagrKernel(self.__xext,nxe)
        instKernel.fillKernel()
        Kernel = instKernel.get_KM()
        instKernel.fillKernelAdapted()        
        KernelAd = instKernel.get_KM() # unclear of whether it is a copy
  #second run a double loop for integrating over [0,xi]
        for i in range(0,nxe-1):
        # set a common factor for D(1)(i,j):
            factor_il = 1.0
            for k in range(0,nxe):
                for l in range(0,nxe):
                   if (l != k): 
                      factor_il = factor_il*KernelAd[i,l]/Kernel[k,l]
                self.__Btable[i,k]=factor_il
        for k in range(0,nxe): #copy
            self.__Btable[nxe-1,k]=self.__Btable[nxe-2,k]


    def checkBtable2():
        nxe=self.__nxe
        # compare with pre calculated values if nxe=4
        if (nxe==4):
            for i in range(0,nxe-1):
                for k in range(0,nxe):
                    if (fabs(self.__Btable[i,k]-manual_Btable2[i,k])>1e-3):
                        print(str(i)+" "+str(k)+" "+str(self.__Btable[i,k]))
                        break

    def __str__(self):
        nxe=self.__nx
        for i in range(0,nxe):
            print(str(i)+"    ")
            for j in range(0,nxe):
                print(self.__Btable[i,j]);
            print("\n")

# end class MatPseudoPoints
