import numpy as np
g2= lambda x,y: np.sin(x)/(x+.1)+x*np.exp(y)+.25*y
P1.name="linear"
V2=Fd(S1,P1, name="2D")
h2=fd.inFd(V2)
tau=[]
tau.append( S1.cmptau() )
tau.append([1., 4.])
Ni=[len(t) for t in tau]
N=np.prod(Ni)
#
# Use FORTRAN order for continuous elements in memory
# when creating arrays. Double precision arithmetics
# is used as default in the fortran libraries. If simple
# precision is requested, specify dtype=np.float32. However,
# this will cause copies of data due to element type
# mismatch.
gtau=np.zeros((N,), order='F')
 for x,i in zip(itt.product(*tau), \
... itt.product(*[range(nx) for nx in Ni])):
... gtau[getidx(i,Ni)]=g2(*x)
h2.cmpcoef(tau,gtau)
# verify interpolation at data sites
for x in itt.product(*tau):
... fx,gx=h2(x),g2(*x)
... print("%+13.6e %+13.6e %+7.2f"%(fx,gx,(1.-fx/gx)*100.))
#+3.733255e+00 +3.733255e+00 +0.00
#+5.636312e+01 +5.636312e+01 +0.00
#+6.119562e+00 +6.119562e+00 +0.00
#+1.106293e+02 +1.106293e+02 +0.00
#+8.450368e+00 +8.450368e+00 +0.00
#+1.648400e+02 +1.648400e+02 +0.00
#+1.093854e+01 +1.093854e+01 +0.00
#+2.192080e+02 +2.192080e+02 +0.00
print(h2)
#fd \in Fd(S(k=2,t=[1.0, 1.0, 2.0, 3.0, 4.0, 4.0]), PP(k=2,xi=[1.0, 4.0],nu=[]))
#coef =
#[ 3.73325545 6.11956243 8.45036807 10.93854134 17.5432894
#34.8365788 52.1298682 69.42315761]
#tau = [array([ 1., 2., 3., 4.]), [1.0, 4.0]]
#gtau = [ 3.73325545 6.11956243 8.45036807 10.93854134 56.36312366
#110.62929884 164.83997268 219.20801416]
#---
print("%13.6e"%h2([1.89,3.05]))
#7.337242e+01
