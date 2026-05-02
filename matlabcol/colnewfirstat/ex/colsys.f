c
c----------------------------------------------------------------
c
c  problem 1 - see companion paper
c
      implicit real*8 (a-h,o-z)
      double precision fspace(2000), zeta(4), tol(2), z(4), u(4), err(4)
      integer ispace(200), m(1), ipar(11), ltol(2)
      external fsub, dfsub, gsub, dgsub
c
      write (6,99)
   99 format(1h1, 35h example of a simple problem setup.
     .         /  46h  uniformly loaded beam of variable stiffness,
     .         /  32h  simply supported at both ends. /)
c
c     one differential equation of order 4.
      m(1) = 4
c     give location of boundary conditions
      zeta(1) = 1.d0
      zeta(2) = 1.d0
      zeta(3) = 2.d0
      zeta(4) = 2.d0
c     set up parameter array.
c     use default values for all parameters except for initial
c     mesh size, no. of tolerances and sizes of work arrays
      do 10 i=1,11
   10   ipar(i) = 0
      ipar(3) = 1
      ipar(4) = 2
      ipar(5) = 2000
      ipar(6) = 200
c     two error tolerances (on u and its second derivative)
      ltol(1) = 1
      ltol(2) = 3
      tol(1) = 1.d-7
      tol(2) = 1.d-7
c
      call colsys (1, m, 1.d0, 2.d0, zeta, ipar, ltol, tol,
     .             dummy, ispace, fspace, iflag, fsub,
     .             dfsub, gsub, dgsub, dummy)
c
      if (iflag .ne. 1)  stop
c     calculate the error at 101 points using the known
c     exact solution
      x = 1.d0
      do 20 i=1,4
   20    err(i) = 0.d0
      do 40 j=1,101
         call appsln (x, z, fspace, ispace)
         call exact (x, u)
         do 30 i=1,4
   30       err(i) = dmax1(err(i), dabs(u(i)-z(i)))
   40    x = x + .01
      write(6,100) (err(i),i=1,4)
  100 format(/27h error tolerances satisfied//22h the exact errors are,
     .       / 7x,4d12.4)
      stop
      end
c................................................................
      subroutine fsub (x, z, f)
      double precision z(4), f(1), x
      f(1) = (1.d0 - 6.d0*x**2*z(4) - 6.d0*x*z(3)) / x**3
      return
      end
c................................................................
      subroutine dfsub (x, z, df)
      double precision z(4), df(1,4), x
      df(1,1) = 0.d0
      df(1,2) = 0.d0
      df(1,3) = -6.d0/x**2
      df(1,4) = -6.d0/x
      return
      end
c................................................................
      subroutine gsub (i, z, g)
      double precision z(4), g
      go to (1, 2, 1, 2), i
    1 g = z(1) - 0.d0
      return
    2 g = z(3) - 0.d0
      return
      end
c................................................................
      subroutine dgsub (i, z, dg)
      double precision z(4), dg(4)
      do 10 j=1,4
   10    dg(j) = 0.d0
      go to (1, 2, 1, 2), i
    1 dg(1) = 1.d0
      return
    2 dg(3) = 1.d0
      return
      end
c................................................................
      subroutine exact(x, u)
      implicit real*8 (a-h,o-z)
      double precision u(4)
c     exact solution
      u(1) = .25d0* (10.d0*dlog(2.d0)-3.d0) * (1.d0-x) +
     .      .5d0* (1.d0/x+ (3.d0+x)*dlog(x) - x)
      u(2) = -.25d0* (10.d0*dlog(2.d0) - 3.d0) + .5d0 *
     .       (-1.d0/x/x + dlog(x) + (3.d0+x)/x - 1.d0)
      u(3) = .5d0 * (2.d0/x**3 + 1.d0/x -3.d0/x/x)
      u(4) = .5d0 * (-6.d0/x**4 - 1.d0/x/x + 6.d0/x**3)
      return
      end
c
c----------------------------------------------------------------
c
c  problem 2 - see companion paper
c
      implicit real*8 (a-h,o-z)
      double precision zeta(4), fspace(40000), tol(4), z(4)
      integer m(2), ipar(11), ispace(2500), ltol(4)
      common eps, dmu, eps4mu, gamma, xt
      external solutn,fsub,dfsub,gsub,dgsub
c     define constants, print a heading.
      gamma = 1.1d0
      eps = .01d0
      dmu = eps
      eps4mu = eps**4/dmu
      xt = dsqrt(2.d0*(gamma-1.d0)/gamma)
      write (6,100) gamma, xt, eps, dmu, eps4mu
  100 format(1h1,27hdimpling of spherical caps.
     ./8h gamma =,f7.2/ 6h  xt =,d12.5/ 6h eps =,d12.5/ 6h  mu =,d12.5/
     .       12h eps**4/mu =,d12.5)
c     define no. of differential equations.
      ncomp = 2
c     orders
      m(1) = 2
      m(2) = 2
c     interval ends
      aleft = 0.d0
      aright = 1.d0
c     locations of side conditions
      zeta(1) = 0.d0
      zeta(2) = 0.d0
      zeta(3) = 1.d0
      zeta(4) = 1.d0
c     ipar  values
c     a nonlinear problem
      ipar(1) = 1
c     4 collocation points per subinterval
      ipar(2) = 4
c     initial uniform mesh of 10 subintervals
      ipar(3) = 10
      ipar(8) = 0
c     dimension of real work array  fspace  is 40000
      ipar(5) = 40000
c     dimension of integer work array  ispace  is 2500
      ipar(6) = 2500
c     (these dimensions of  fspace  and  ispace
c      enable  colsys  to use meshes of up to 192 intervals.)
c     print full output.
      ipar(7) = -1
c     initial approximation for nonlinear iteration is provided
c     in  solutn
      ipar(9) = 1
c     a regular problem
      ipar(10) = 0
c     no fixed points in the mesh
      ipar(11) = 0
c     tolerances on  all components
      ipar(4) = 4
      do 10 i=1,4
        ltol(i) = i
   10   tol(i) = 1.d-5
c     call  colsys
      call colsys (ncomp, m, aleft, aright, zeta, ipar, ltol,
     .             tol, fixpnt, ispace, fspace, iflag,
     .             fsub,dfsub,gsub,dgsub,solutn)
c     print values of the obtained approximate solution at points
c     x = 0,.05, ..., 1.
      x = 0.d0
      write (6,201)
  201 format (1h1,44h       x             phi           dphi      ,
     .      23h      psi          dpsi /)
  202 format (6x, f5.2, 4x, 6d15.5)
      np1 = 21
      do 555 iii=1,np1
        call appsln (x,z,fspace,ispace)
        write (6,202) x, z
        x = x + .05d0
  555 continue
      stop
      end
c................................................................
      subroutine solutn (x, z, dmval)
      implicit real*8 (a-h,o-z)
      common eps, dmu, eps4mu, gamma, xt
      dimension z(4) , dmval(2)
      cons = gamma * x * (1.d0-.5d0*x*x)
      dcons = gamma * (1.d0 - 1.5d0*x*x)
      d2cons = -3.d0 * gamma * x
      if (x .gt. xt) go to 10
      z(1) = 2.d0 * x
      z(2) = 2.d0
      z(3) = -2.d0*x + cons
      z(4) = -2.d0 + dcons
      dmval(2) = d2cons
      go to 20
   10 z(1) = 0.d0
      z(2) = 0.d0
      z(3) = -cons
      z(4) = -dcons
      dmval(2) = -d2cons
   20 dmval(1) = 0.d0
      return
      end
c................................................................
      subroutine fsub (x, z, f)
      implicit real*8 (a-h,o-z)
      dimension z(4), f(2)
      common eps, dmu, eps4mu, gamma, xt
      f(1) = z(1)/x/x - z(2)/x + (z(1) - z(3)*(1.d0-z(1)/x) -
     .       gamma*x*(1.d0-x*x/2.)) / eps4mu
      f(2) = z(3)/x/x - z(4)/x + z(1)*(1.d0-z(1)/2.d0/x) / dmu
      return
      end
c................................................................
      subroutine dfsub (x, z, df)
      implicit real*8 (a-h,o-z)
      dimension z(4), df(2,4)
      common eps, dmu, eps4mu, gamma, xt
      df(1,1) = 1.d0/x/x +(1.d0 + z(3)/x) / eps4mu
      df(1,2) = -1.d0/x
      df(1,3) = -(1.d0-z(1)/x) / eps4mu
      df(1,4) = 0.d0
      df(2,1) = (1.d0 - z(1)/x) / dmu
      df(2,2) = 0.d0
      df(2,3) = 1.d0/x/x
      df(2,4) = -1.d0/x
      return
      end
c................................................................
      subroutine gsub (i, z, g)
      implicit real*8 (a-h,o-z)
      dimension z(4)
      go to (1, 2, 1, 3), i
    1 g = z(1)
      return
    2 g = z(3)
      return
    3 g = z(4) - .3d0*z(3) + .7d0
      return
      end
c................................................................
      subroutine dgsub (i, z, dg)
      implicit real*8 (a-h,o-z)
      dimension z(4), dg(4)
      do 10 j=1,4
   10    dg(j) = 0. d0
      go to (1, 2, 1, 3), i
    1 dg(1) = 1.d0
      return
    2 dg(3) = 1.d0
      return
    3 dg(4) = 1.d0
      dg(3) = -.3d0
      return
      end
c
c----------------------------------------------------------------
c
c  problem 3 - see companion paper
c
      implicit real*8 (a-h,o-z)
      double precision zeta(5), fspace(40000), tol(2), sval(3), elval(3)
      integer m(2), ispace(2500), ltol(2), ipar(11)
      double precision z(5)
      common en, s, el, cons
      external fsub, dfsub, gsub, dgsub, solutn
      data sval/.2d0, .1d0, .05d0/, elval/60.d0, 120.d0, 200.d0/
c
      en = .2d0
      cons = .5d0 * (3.d0-en)
      ncomp = 2
      m(1) = 2
      m(2) = 3
      aleft = 0.d0
      aright = 1.d0
c
      zeta(1) = 0.d0
      zeta(2) = 0.d0
      zeta(3) = 0.d0
      zeta(4) = 1.d0
      zeta(5) = 1.d0
c
      ipar(1) = 1
      ipar(2) = 4
      ipar(3) = 10
      ipar(4) = 2
      ipar(5) = 40000
      ipar(6) = 2500
      ipar(7) = -1
      ipar(8) = 0
      ipar(9) = 1
      ipar(10) = 0
      ipar(11) = 0
c
      ltol(1) = 1
      ltol(2) = 3
      tol(1) = 1.d-5
      tol(2) = 1.d-5
c
c     solve a chain of 3 problems
      do 777 ijk = 1,3
         s = sval(ijk)
         el = elval(ijk)
         if (ijk .eq. 1)  go to 701
c        set continuation parameters
         ipar(9) = 3
         ipar(3) = ispace(1)
  701    continue
         write (6,100) en, s, el
  100    format(1h1,38h rotating flow over a stationary disk.
     .          /19h  parameters -  n =,f5.2, 6h   s =,f5.2,
     .          6h   l =,f6.1/)
c
         call colsys (ncomp, m, aleft, aright, zeta, ipar, ltol,
     .                tol, fixpnt, ispace, fspace, iflag,
     .                fsub, dfsub, gsub, dgsub, solutn)
c
         if (iflag .ne. 1)  stop
c        print values of the obtained approximate solution at points
c        x = 0,1,2, ..., l.
         is6 = ispace(6)
         is5 = ispace(1) + 2
         x = 0.d0
         write (6,201)
  201    format (1h1,44h       x              g              dg      ,
     .         38h       h             dh            d2h/)
  202    format (6d15.5)
         np1 = el + 1.5d0
         do 555 iii=1,np1
           call approx (ii,x,z,fspace(is6),fspace(1),ispace(1),
     .                  fspace(is5),ispace(2),ncomp,m,ispace(4),1,dm,0)
           xl = x * el
           z(2) = z(2) / el
           z(4) = z(4) / el
           z(5) = z(5) / el / el
           write (6,202) xl, z
           x = x + 1.d0 / el
  555    continue
  777 continue
      stop
      end
c................................................................
      subroutine solutn (x, z, dmval)
      implicit real*8 (a-h,o-z)
      common en, s, el, cons
      double precision z(5), dmval(2)
      ex = dexp(-el*x)
      z(1) = 1.d0 - ex
      z(2) = el * ex
      z(3) = -el**2 * x**2 * ex
      z(4) = (el**3 *x**2 - 2.d0 * el**2 * x) * ex
      z(5) = (-el**4 * x**2  + 4.d0 * el**3 * x - 2.d0 * el**2)*ex
      dmval(1) = -el * z(2)
      dmval(2) = (el**5*x*x - 6.d0*el**4*x + 6.d0*el**3) * ex
      return
      end
c................................................................
      subroutine fsub (x, z, f)
      implicit real*8 (a-h,o-z)
      double precision z(1), f(1)
      common en, s, el, cons
      f(1) = -el * (cons * z(3) * z(2) + (en-1.d0) * z(4) * z(1))
     .       + el**2 * s * (z(1)-1.d0)
      f(2) = -el * (cons * z(3) * z(5) + en * z(4)**2) +
     .       el**2 * s * z(4) + el**3 * (1.d0-z(1)**2)
      return
      end
c................................................................
      subroutine dfsub (x, z, df)
      implicit real*8 (a-h,o-z)
      double precision z(1), df(2,1)
      common en, s, el, cons
      df(1,1) = -el * (en-1.d0) * z(4) + el**2 * s
      df(1,2) = -el * cons * z(3)
      df(1,3) = -el * cons * z(2)
      df(1,4) = -el * (en-1.d0) * z(1)
      df(1,5) = 0.d0
      df(2,1) = -el**3 * 2.d0 * z(1)
      df(2,2) = 0.d0
      df(2,3) = -el * cons * z(5)
      df(2,4) = -el * en * 2.d0 * z(4) + el**2 * s
      df(2,5) = -el * cons * z(3)
      return
      end
c................................................................
      subroutine gsub (i, z, g)
      double precision z(1), g
      go to (1, 2, 3, 4, 3), i
    1 g = z(1)
      return
    2 g = z(3)
      return
    3 g = z(4)
      return
    4 g = z(1) - 1.d0
      return
      end
c................................................................
      subroutine dgsub (i, z, dg)
      double precision z(1), dg(1)
      do 10 j=1,5
   10    dg(j) = 0.d0
      go to (1, 2, 3, 1, 3), i
    1 dg(1) = 1.d0
      return
    2 dg(3) = 1.d0
      return
    3 dg(4) = 1.d0
      return
      end
