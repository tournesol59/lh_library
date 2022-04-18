import numpy as num
import packcollocpts as cpts

#test third-order
xtest=num.zeros((3))
for i in range(0,3):
    xtest[i]=cpts.diff_nodes[3,i]
#list
params = [3,-1.0,1.0,1] 
# first derivative
instMatDiff = cpts.MatDiffPoints(params,xtest)
instMatDiff.setallDiffMat()
instMatDiff.__str__()

instMatPseudo = cpts.MatPseudoPoints(params,xtest)
instMatPseudo.setButcherMat()
#instMatPseudo.checkBtable2()
instMatPseudo.__str__()

