15c15
< def newmesh(xi,xiold,fixpnt):
---
> def newmesh(xi,xiold,fixpnt,dims):
18,21c18,21
<     Npold=8 #args['Np']
<     Nfixpnt=3 #args['Nfixpnt']
<     Np=2*Npold-1
< 
---
>     Npold=dims[0] #8 #args['Np']
>     Nfixpnt=dims[1] #3 #args['Nfixpnt']
>     Np=2*Npold+2
>     dims[2]=Np
42,43c42,44
<     # tip : fix points are part of the mesh halving
<         while (i < Npold-1) and (xiold[i] < xright-PRECIS): 
---
>     # fix points should not be part of the mesh halving
>         while (i < Npold-1) and (xiold[i+1] < xright-PRECIS):
> #            print("loop while here: "+str(xright))
45a47
>             dx=xiold[i+1]-xiold[i]
48d49
<                 dx=xiold[i+1]-xiold[i]
52c53
<             else: 
---
>             elif (i >= Npold-1):
53a55,78
>         #now treat the case a mid point must be inserted just before or 
>         # just after a fix pnt:
>         dx=xiold[i+2]-xiold[i]
>         #just before:
>  #       print("just before just before: "+str(xiold[i+2])+" "+ str(xright))
>         if (i < Npold-1) and (xiold[i+2] > xright+PRECIS) and ((xiold[i]+dx/2) < xright-PRECIS):
> #            print("if just before here "+str(i)+" "+str(xright)) #debug
>             xi[k]=xiold[i]
>             k=k+1
>             xi[k]=xiold[i]+dx/2
>             k=k+1
>             xi[k]=xright
>             k=k+1
>             i=i+1
>          #just after 
>         elif (i < Npold-1) and (xiold[i+2] > xright+PRECIS) and ((xiold[i]+dx/2) > xright+PRECIS):
>  #           print("if just after here: "+str(i)+" "+str(xright)) #debug
>             xi[k]=xiold[i]
>             k=k+1
>             xi[k]=xright
>             k=k+1
>             xi[k]=xiold[i]+dx/2
>             k=k+1
>             i=i+1
67a93
>   #          print("k: "+str(k)+" xi[k-1]: "+str(xi[k-1]))
69,70c95
<             k=k+1 
<         
---
>             k=k+1  
72c97
<             print("we are here: last point "+str(k)+" "+str(xi[k]))
---
>            # print("we are here: last point "+str(k)+" "+str(xi[k])) #debug
92,93c117,118
< 
< newmesh(XI,XIOLD,FIXPNT)
---
> dimsize=[8,3,1]
> newmesh(XI,XIOLD,FIXPNT,dimsize)
96,97c121,123
< print("Np= "+str(NP))
< for i in range(0,NP):
---
> if (NP != dimsize[2]):
>     print("error Np: "+str(NP)+" != "+str(dimsize[2]))
> for i in range(0,dimsize[2]):
