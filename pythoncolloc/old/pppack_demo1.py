from pppack import *
k,t=2,[1., 1., 2., 3., 4., 4.]
S1 = S(k,t)
xi = [1., 4.]
P1 = PP(k,xi)
F1 = S1 * P1
print(F1)
#F^(2) = F0 x F1
#with,
#F0 = S_{ k=2, t }, with
#t = [1.000,1.000,2.000,3.000,4.000,4.000]
#F1 = PP_{ k=2, xi, nu }, with
#xi = [1.000,4.000], nu = [-,,-]
s0 = s.inS(S1) # = S(k,t) = s(k,t)
print(s0)
#s \in S(k=2,t=[1.0, 1.0, 2.0, 3.0, 4.0, 4.0])
#coef (to be determined)
