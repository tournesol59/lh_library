# package legendre polynomials
# cp Frederic Peugny, LTS

#package legendre polynomial
import numpy as num
#from matplotlib import pyplot as plt


class LegendPoly(object):
    def __init__(self, n, a, b):
        self.deg=n
        self.xa=a
        self.xb=b
        self.coeffs=[]
        cn=1
        for j in range(1,n+1):
            cn=(j)/(2*j)
        for k in range(0,deg):
            if (n-n/2==n/2): 
                ck=1
            else:
                ck=-1
            for j in range(1,k+1):
                ck=ck * (-1)*(n-j)*(n+k-j)/(j)/(j) # formula of binomial for legendre
            self.coeffs.append(ck*cn)


    def print(self):
        for i in range(0,deg+1):
            print self.coeffs[i]
            print(" x")
            print(i)
            print(" + ")

