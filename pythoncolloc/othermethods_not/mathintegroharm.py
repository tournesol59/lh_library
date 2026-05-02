import numpy as np
import math as math

#copy from packcollocpts.py to compute the harmonic 2-order finite element matrix
diff_nodes4=np.array([-0.86114,-0.340,0.340,0.86114, 0.,0.,0.])
diff_quad4=np.array([0.34785, 0.65215, 0.65215, 0.34785])

herm_coeffs=np.array([[1., 0., -3.0, 2.0],
                      [0., 1.0, -2.0, 1.0],
                      [0., 0., 3.0, -2.0],
                      [0., 0., -1.0, 1.0]])

herm_dercoeffs=np.array([0., -6., 6., 0.],
                     [1.0, -4., 3., 0.], 
                     [0., 6., -6., 0.],
                     [0., -2., 3., 0.])

def calc_FEmatrix(deltah, stiff):
    # deltah is the mesh step
   matrixK=num.zeros((4,4))
   # phi,psi,eta,rho = 4 hermite polynoms
   # nodes i,i+1
   matrixD=num.zeros((4,4))
   matrixP=num.zeros((4,4))
# phi(i)*phi(i)
   prodpoly=num.zeros((7))
   for l in range(0,4):
       #express coeff of square(phi)
       for r in range(0,4):
          prodpoly[l+r]=prodpoly[l+r]+herm_coeff[0,l]*herm_coeff[0,r]
   # integrate with gaussian quadrature
   for i in range(0,4);
       matrixP[0,0]=matrixP[0,0]+num.polyval(prodpol, diff_nodes4[i]) *diff_quad4[i]
# phi(i)*psi(i)
   prodpoly=num.zeros((7))
   for l in range(0,4):
       #express coeff of (phi)(psi)
       for r in range(0,4):
          prodpoly[l+r]=prodpoly[l+r]+herm_coeff[0,l]*herm_coeff[1,r]
   # integrate with gaussian quadrature
   for i in range(0,4);
       matrixP[0,2]=matrixP[0,0]+num.polyval(prodpol, diff_nodes4[i]) *diff_quad4[i]
# phi(i)*eta(i)
   prodpoly=num.zeros((7))
   for l in range(0,4):
       #express coeff of (phi)(eta)
       for r in range(0,4):
          prodpoly[l+r]=prodpoly[l+r]+herm_coeff[0,l]*herm_coeff[2,r]
   # integrate with gaussian quadrature
   for i in range(0,4);
       matrixP[2,2]=matrixP[0,0]+num.polyval(prodpol, diff_nodes4[i]) *diff_quad4[i]
# phi(i)*rho(i)
   prodpoly=num.zeros((7))
   for l in range(0,4):
       #express coeff of (phi)(rho)
       for r in range(0,4):
          prodpoly[l+r]=prodpoly[l+r]+herm_coeff[0,l]*herm_coeff[3,r]
   # integrate with gaussian quadrature
   for i in range(0,4);
       matrixP[0,3]=matrixP[0,0]+num.polyval(prodpol, diff_nodes4[i]) *diff_quad4[i]

# psi(i)*psi(i)
   prodpoly=num.zeros((7))
   for l in range(0,4):
       #express coeff of square(psi)
       for r in range(0,4):
          prodpoly[l+r]=prodpoly[l+r]+herm_coeff[1,l]*herm_coeff[1,r]
   # integrate with gaussian quadrature
   for i in range(0,4);
       matrixP[1,1]=matrixP[1,1]+num.polyval(prodpol, diff_nodes4[i]) *diff_quad4[i]
# psi(i)*eta(i)
   prodpoly=num.zeros((7))
   for l in range(0,4):
       #express coeff of (psi)(eta)
       for r in range(0,4):
          prodpoly[l+r]=prodpoly[l+r]+herm_coeff[1,l]*herm_coeff[2,r]
   # integrate with gaussian quadrature
   for i in range(0,4);
       matrixP[1,2]=matrixP[1,2]+num.polyval(prodpol, diff_nodes4[i]) *diff_quad4[i]
# psi(i)*rho(i)
   prodpoly=num.zeros((7))
   for l in range(0,4):
       #express coeff of (psi)(rho)
       for r in range(0,4):
          prodpoly[l+r]=prodpoly[l+r]+herm_coeff[1,l]*herm_coeff[3,r]
   # integrate with gaussian quadrature
   for i in range(0,4);
       matrixP[1,3]=matrixP[1,3]+num.polyval(prodpol, diff_nodes4[i]) *diff_quad4[i]
# eta(i)*eta(i)
   prodpoly=num.zeros((7))
   for l in range(0,4):
       #express coeff of square(eta)
       for r in range(0,4):
          prodpoly[l+r]=prodpoly[l+r]+herm_coeff[2,l]*herm_coeff[2,r]
   # integrate with gaussian quadrature
   for i in range(0,4);
       matrixP[2,2]=matrixP[2,2]+num.polyval(prodpol, diff_nodes4[i]) *diff_quad4[i]
       # eta(i)*eta(i)
   prodpoly=num.zeros((7))
   for l in range(0,4):
       #express coeff of square(eta)
       for r in range(0,4):
          prodpoly[l+r]=prodpoly[l+r]+herm_coeff[2,l]*herm_coeff[2,r]
   # integrate with gaussian quadrature
   for i in range(0,4);
       matrixP[2,2]=matrixP[2,2]+num.polyval(prodpol, diff_nodes4[i]) *diff_quad4[i]


