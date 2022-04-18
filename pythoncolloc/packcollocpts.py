#package of mesh-dependent derivatives of (Legendre) polynoms
# cp Frederic Peugny, LTS
import numpy as num
#from matplotlib import pyplot as plt

diff_nodes=num.zeros((7,7))
#high precision values
diff_nodes[0,:]=[1E9] # Inf
diff_nodes[1,:]=[0.,0.,0.,0.,0.,0.,0.]
diff_nodes[2,:]=[-0.57735,0.57735,0.,0.,0.,0.,0.]
diff_nodes[3,:]=[-0.7746, 0., 0.7746, 0.,0.,0.,0.]
diff_nodes[4,:]=[-0.86114,-0.340,0.340,0.86114, 0.,0.,0.]
diff_nodes[5,:]=[-0.90618,-0.53857,0.,0.53857, 0.90618, 0.,0.]
#low precision values
diff_nodes[6,:]=[-0.92,-0.65,-0.24,0.24,0.65,0.92, 0.]
#diff_nodes[7,:]=[-0.945,-0.74,-0.39,0.,0.39,0.74,0.945]

class MatLagrKernel(object):
    def __init__(x,nx):
        KM=num.zeros(nx,nx-1)
        self.__fillKernel()

    def fillKernel():
        for i in range([0,nx]):
            for j in range([0,i-1]):
                KM[i,j]=xext[i]-xext[j]
            for j in range([i+1,nx-1]):
                KM[i,j]=xext[i]-xext[j]

    def get_KM(self):
        return self.__KM  # test this!
#end class MatLagrKernel

class MatDiffPoints(object):
    def __init__(m,x,param):
    # m: order of derivative
    # x: numpy array of points (the mesh)
    # param(1): (tuple) no. of collocation points per subinterval w/o extrems
    # param(2): left value
    # param(3): right value
    # param(4): no. of subintervals in the initial mesh
    # param(...) think about other params..
        xext=[-1.0,x,1.0]
        nx=num.size(xext)
        DM = num.zeros((nx+1, nx-1))

    def setallDM():    
    # here implement the lagrange difference and products
    # to have DM(i,j)=deriv(Lagr{j}, dt)(ti)
     #first initialize a Kernel of ALL possible differences combination of the x
        instKernel = MatLagrKernel(x,nx, 0., 1.0)
        Kernel = instKernel.get_KM()
     #second run a double loop for DM(i,j)
        for i in range([0,nx+1]):
         # set a common factor for D(1)(i,j):
            alphac_i = 1.0
            for l in range([0,nx]):
               if (l != i): # note that alpha_c means completed by (ti-tj)
                    alphac_i *= Kernel(i,l)
            for j in range([0,nx+1]):
               if (j != i):
         # now we have to compute alpha_j to complete DM^(1)(i,j)=al_j/alc_i
                    alpha_j = 1.0
                    for l in range([0,nx]):
                        if (l != j) and (l != i): 
         # this test is necessary to have not simplification by (ti-tj)
                            alpha_j *= Kernel(j,l)
                    DM[i,j] = alpha_j / alphac_i
               else: 
                 if (j==i):
         # next we have to complete DM^(1)(i,i)=sum 1/(ti-tl)
                    term_ii = 0.0
                    for l in range([0,nx]):
         #               if (l != i):
         # this test is not necessary to not divide by zero because Kernel(i,l) is already safe
                        term_ii += 1.0/Kernel(i,l)
                    DM[i,i] = term_ii
 #end method setallDM() of class MatDiffPoints.

    def __str__():
        for i in range([0,nx+1]):
            print(i+"    ")
            for j in range([0,nx-1]):
                print(DM[i,j]+"    ");
            print("\n")

# end class MatDiffPoints


# test program, to comment if package is vaidated:
xtest=num.zeros((4))
for i in range(0,4):
    xtest[i]=diff_nodes[3,i]
testparam = (3,-1.0,1.0,1) #tuple

instMatDiff = MatDiffPoints(1,xtest) # tried with args (1,xtest,testparam). failed
instMatDiff.setallDM
instMatDiff.__str__
