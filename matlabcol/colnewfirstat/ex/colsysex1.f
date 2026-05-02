c  fred created this file and modified the write to output results
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
      write(6,200)
  200 format(1h1, 35h solution calculated values 1,4 ix.
     .         /  46h ---- x   u1  z1  u2  z2  u3  z3  u4  z4 ----/)
c     calculate the error at 101 points using the known
c     exact solution
      x = 1.d0
      do 20 i=1,4
   20    err(i) = 0.d0
      do 40 j=1,101
         call appsln (x, z, fspace, ispace)
         call exact (x, u)
         write(6,210) x, u(1), z(1), u(2), z(2), u(3), z(3), u(4), z(4)
  210    format(/ 7x,9d12.4)
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
      subroutine dgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(1),info
      double precision a(lda,1)
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
      double precision t
      integer idamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end
      subroutine dgesl(a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(1),job
      double precision a(lda,1),b(1)
c
c     dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c
c     internal variables
c
      double precision ddot,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end
      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),da
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
   20 continue
      do 30 i = 1,n
        dy(i) = dy(i) + da*dx(i)
   30 continue
      return
      end
      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),dtemp
      integer i,incx,incy,ix,iy,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
   20 continue
      do 30 i = 1,n
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      ddot = dtemp
      return
      end
      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     jack dongarra, linpack, 3/11/78.
c
      double precision da,dx(1)
      integer i,incx,n,nincx
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
   20 continue
      do 30 i = 1,n
        dx(i) = da*dx(i)
   30 continue
      return
      end
      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. dabsolute value.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n .lt. 1 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
      double precision function epslon (x)
      double precision x
c
c     estimate unit roundoff in quantities of size x.
c
      double precision a,b,c,eps
c
c     this program should function properly on all systems
c     satisfying the following two assumptions,
c        1.  the base used in representing dfloating point
c            numbers is not a power of three.
c        2.  the quantity  a  in statement 10 is represented to 
c            the accuracy used in dfloating point variables
c            that are stored in memory.
c     the statement number 10 and the go to 10 are intended to
c     force optimizing compilers to generate code satisfying 
c     assumption 2.
c     under these assumptions, it should be true that,
c            a  is not exactly equal to four-thirds,
c            b  has a zero for its last bit or digit,
c            c  is not exactly equal to one,
c            eps  measures the separation of 1.0 from
c                 the next larger dfloating point number.
c     the developers of eispack would appreciate being informed
c     about any systems where these assumptions do not hold.
c
c     *****************************************************************
c     this routine is one of the auxiliary routines used by eispack iii
c     to avoid machine dependencies.
c     *****************************************************************
c
c     this version dated 4/6/83.
c
      a = 4.0d0/3.0d0
   10 b = a - 1.0d0
      c = b + b + b
      eps = dabs(c-1.0d0)
      if (eps .eq. 0.0d0) go to 10
      epslon = eps*dabs(x)
      return
      end
      subroutine dmxpy (n1, y, n2, ldm, x, m)
      double precision y(*), x(*), m(ldm,*)
c
c   purpose:
c     multiply matrix m times vector x and add the result to vector y.
c
c   parameters:
c
c     n1 integer, number of elements in vector y, and number of rows in
c         matrix m
c
c     y double precision(n1), vector of length n1 to which is added 
c         the product m*x
c
c     n2 integer, number of elements in vector x, and number of columns
c         in matrix m
c
c     ldm integer, leading dimension of array m
c
c     x double precision(n2), vector of length n2
c
c     m double precision(ldm,n2), matrix of n1 rows and n2 columns
c
c ----------------------------------------------------------------------
c
c   cleanup odd vector
c
      j = mod(n2,2)
      if (j .ge. 1) then
         do 10 i = 1, n1
            y(i) = (y(i)) + x(j)*m(i,j)
   10    continue
      endif
c
c   cleanup odd group of two vectors
c
      j = mod(n2,4)
      if (j .ge. 2) then
         do 20 i = 1, n1
            y(i) = ( (y(i))
     $             + x(j-1)*m(i,j-1)) + x(j)*m(i,j)
   20    continue
      endif
c
c   cleanup odd group of four vectors
c
      j = mod(n2,8)
      if (j .ge. 4) then
         do 30 i = 1, n1
            y(i) = ((( (y(i))
     $             + x(j-3)*m(i,j-3)) + x(j-2)*m(i,j-2))
     $             + x(j-1)*m(i,j-1)) + x(j)  *m(i,j)
   30    continue
      endif
c
c   cleanup odd group of eight vectors
c
      j = mod(n2,16)
      if (j .ge. 8) then
         do 40 i = 1, n1
            y(i) = ((((((( (y(i))
     $             + x(j-7)*m(i,j-7)) + x(j-6)*m(i,j-6))
     $             + x(j-5)*m(i,j-5)) + x(j-4)*m(i,j-4))
     $             + x(j-3)*m(i,j-3)) + x(j-2)*m(i,j-2))
     $             + x(j-1)*m(i,j-1)) + x(j)  *m(i,j)
   40    continue
      endif
c
c   main loop - groups of sixteen vectors
c
      jmin = j+16
      do 60 j = jmin, n2, 16
         do 50 i = 1, n1
            y(i) = ((((((((((((((( (y(i))
     $             + x(j-15)*m(i,j-15)) + x(j-14)*m(i,j-14))
     $             + x(j-13)*m(i,j-13)) + x(j-12)*m(i,j-12))
     $             + x(j-11)*m(i,j-11)) + x(j-10)*m(i,j-10))
     $             + x(j- 9)*m(i,j- 9)) + x(j- 8)*m(i,j- 8))
     $             + x(j- 7)*m(i,j- 7)) + x(j- 6)*m(i,j- 6))
     $             + x(j- 5)*m(i,j- 5)) + x(j- 4)*m(i,j- 4))
     $             + x(j- 3)*m(i,j- 3)) + x(j- 2)*m(i,j- 2))
     $             + x(j- 1)*m(i,j- 1)) + x(j)   *m(i,j)
   50    continue
   60 continue
      return
      end
      DOUBLE PRECISION FUNCTION RAN( ISEED )
*
*     modified from the LAPACK auxiliary routine 10/12/92 JD
*  -- LAPACK auxiliary routine (version 1.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
*     ..
*
*  Purpose
*  =======
*
*  DLARAN returns a random real number from a uniform (0,1)
*  distribution.
*
*  Arguments
*  =========
*
*  ISEED   (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd.
*          On exit, the seed is updated.
*
*  Further Details
*  ===============
*
*  This routine uses a multiplicative congruential method with modulus
*  2**48 and multiplier 33952834046453 (see G.S.Fishman,
*  'Multiplicative congruential random number generators with modulus
*  2**b: an exhaustive analysis for b = 32 and a partial analysis for
*  b = 48', Math. Comp. 189, pp 331-344, 1990).
*
*  48-bit integers are stored in 4 integer array elements with 12 bits
*  per element. Hence the routine is portable across machines with
*  integers of 32 bits or more.
*
*     .. Parameters ..
      INTEGER            M1, M2, M3, M4
      PARAMETER          ( M1 = 494, M2 = 322, M3 = 2508, M4 = 2549 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      INTEGER            IPW2
      DOUBLE PRECISION   R
      PARAMETER          ( IPW2 = 4096, R = ONE / IPW2 )
*     ..
*     .. Local Scalars ..
      INTEGER            IT1, IT2, IT3, IT4
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MOD
*     ..
*     .. Executable Statements ..
*
*     multiply the seed by the multiplier modulo 2**48
*
      IT4 = ISEED( 4 )*M4
      IT3 = IT4 / IPW2
      IT4 = IT4 - IPW2*IT3
      IT3 = IT3 + ISEED( 3 )*M4 + ISEED( 4 )*M3
      IT2 = IT3 / IPW2
      IT3 = IT3 - IPW2*IT2
      IT2 = IT2 + ISEED( 2 )*M4 + ISEED( 3 )*M3 + ISEED( 4 )*M2
      IT1 = IT2 / IPW2
      IT2 = IT2 - IPW2*IT1
      IT1 = IT1 + ISEED( 1 )*M4 + ISEED( 2 )*M3 + ISEED( 3 )*M2 +
     $      ISEED( 4 )*M1
      IT1 = MOD( IT1, IPW2 )
*
*     return updated seed
*
      ISEED( 1 ) = IT1
      ISEED( 2 ) = IT2
      ISEED( 3 ) = IT3
      ISEED( 4 ) = IT4
*
*     convert 48-bit integer to a real number in the interval (0,1)
*
      RAN = R*( DBLE( IT1 )+R*( DBLE( IT2 )+R*( DBLE( IT3 )+R*
     $         ( DBLE( IT4 ) ) ) ) )
      RETURN
*
*     End of RAN
*
      END

