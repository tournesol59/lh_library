c................................................................
      subroutine fsub (x, z, f)
      double precision z(2), f(2), x
c     Harmonic oscillator
      f(1) = z(2)
      f(2) = -0.1*z(2) - 16*z(1)
      return
      end
c................................................................
      subroutine dfsub (x, z, df)
      double precision z(2), df(2,2), x
      df(1,1) = 0.d0
      df(1,2) = 1.d0
      df(2,1) = 16.d0
      df(2,2) = -0.1
      return
      end
c     TO COMPLETE BELOW!
c................................................................
      subroutine gsub (i, z, g)
      double precision z(2), g(2)
      go to (1, 2), i
    1 g = z(1) - 0.d0
      return
    2 g = z(1) - 0.d0
      return
      end
c................................................................
      subroutine dgsub (i, z, dg)
      double precision z(2), dg(2,2)
      do 10 j=1,2
   10      do 20 k=1,2
   20    dg(j,k) = 0.d0
      go to (1, 2), i
    1 dg(1,1) = 1.d0
      dg(1,2) = 0.0
      return
    2 dg(2,1) = 1.d0
      dg(2,2) = 0.0
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
