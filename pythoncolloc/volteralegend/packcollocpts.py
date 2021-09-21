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

class MatLagrKernel(object):
    def __init__(self,xext,nx):
            self.__xext=num.zeros((nx))
            self.__KM=num.zeros((nx,nx))
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

    def get_KM(self):
        return self.__KM  # test this!

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
    # param(3): no. of subintervals in the initial mesh
    # param(...) think about other params..
        self.__xext=num.zeros((param[0]+2))
        self.__xext[0]=-1.0 # default left value
        for i in range(1,param[0]+1):
            self.__xext[i]=x[i-1]
        self.__xext[param[0]+1]=1.0 # default right value
        self.__nx=param[0]+2
        self.__DM = num.zeros((self.__nx, self.__nx)) # nx Lagr Points leads to nx polynoms and nx eval pts
         
    def setallDiffMat(self):    
    # here implement the lagrange difference and products
    # to have DM(i,j)=deriv(Lagr{j}, dt)(ti)
    #first initialize a Kernel of ALL possible differences combination of the x
        nxe=self.__nx
        # create a class of Lagrange points difference of size nxe=total no.pt
        instKernel = MatLagrKernel(self.__xext,nxe)
        instKernel.fillKernel()        
        Kernel = instKernel.get_KM() # unclear of whether it is a copy

   #second run a double loop for DM(i,j)
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
                    self.__DM[i,j] = alpha_j / alphac_i
                else:
                    #if (j==i):
         # next we have to complete DM^(1)(i,i)=sum 1/(ti-tl)
                    term_ii = 0.0
                    for l in range(0,nxe):
                        if (l != i):
         # this test is necessary to not divide by zero but Kernel(i,l) could be safe
                            term_ii = term_ii + 1.0/Kernel[i,l]
                    self.__DM[i,i] = term_ii
 #end method setallDM() of class MatDiffPoints.

    def __str__(self):
        nxe=self.__nx
        for i in range(0,nxe):
            print(str(i)+"    ")
            for j in range(0,nxe):
                print(self.__DM[i,j]);
            print("\n")

# end class MatDiffPoints



