#package of Bsplines primitives (nodes, polynom eval) for preview of colsys
# cp Frederic Peugny, LTS
# SIMPLIFICATION: a single interval shall be tested (no internal points guess)
import numpy as num

class ColsysBPoints(object):

    def __init__(self,param,x):	
    # x: numpy array of points (the mesh)
    # param(0): (list) no. of collocation points per subinterval w/o extrems
    # param(1): left value
    # param(2): right value
    # param(3): choice of nodes: 'norm'=in [-1,1] 'scal'=[p1,p2]    
    # param(4) no. of subintervals in the initial mesh.
        nsub=param[4]
        self.__nsub=nsub
        nxe=param[0]+1  # num of inner colloc points + one for one end
        # nxep1=nxe+1 
        self.__nxe=nxe
        nxext=nxe*nsub+1  # local only
        self.__xext=num.zeros(nxe+1) # SIMPLIFICATION: a single interval
        # self.__xext=num.zeros(nxe*nsub+1)
        self.__xext[0]=param[1] # -1.0 default left value
        self.__xext[nxext-1]=param[2] # 1.0 default right value
        if (nsub==1):  # param[3]=='norm'
            for i in range(1,nxe-1):  # a single interval
                self.__xext[i]=self.__xext[0]+i*(self.__xext[nxext-1] - self.__xext[0])/nxe  # yes nxe=nx+1, as nx is num of inner points
        elif (param[3]=='scal'):  # preview for a complete mesh (more sub interv)
            for k in range(1,nsub-1):
                self.__xext[k*nxe]=self.__xext[0]+k*(self.__xext[nxext-1] - self.__xext[0])/nxe 
                for i in range(1,nxe-1):	
                   self.__xext[(k-1)*nxe+i]=self.__xext[(k-1)*nxe]+i*(self.__xext[k*nxe] - self.__xext[(k-1)*nxe])/nxe  # yes (k-1)*nxe is for the k-th interval beginning

    def displayresult(self):
        print("This is correct")
# end class ColsysBPoints
