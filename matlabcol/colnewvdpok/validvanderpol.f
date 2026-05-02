C     EXAMPLE PROGRAM 
C     FROM THESE PACKAGES FROM ASHER, BADER (CHRISTIANSEN, RUSSELL)

C     THESE SUBROUTINES ARE COMPRISED IN COLNEW.F:

C      SUBROUTINE COLNEW (NCOMP, M, ALEFT, ARIGHT, ZETA, IPAR, LTOL,
C     1                   TOL, FIXPNT, ISPACE, FSPACE, IFLAG,
C     2                   FSUB, DFSUB, GSUB, DGSUB, GUESS)
C    EQUIVALENT ENTRY:
C      SUBROUTINE COLNEW (NCOMP, M, ALEFT, ARIGHT, ZETA, IPAR, LTOL,
C     1                   TOL, FIXPNT, ISPACE, FSPACE, IFLAG,
C     2                   FSUB, DFSUB, GSUB, DGSUB, GUESS)

C     THESE SUBROUTINES ARE PROGRAMMED HERE AND DECLARED "EXTERN" IN 
C     COLNEW.F

      SUBROUTINE FSUB (XCOL, ZVAL, F)
      IMPLICIT NONE
      REAL,POINTER :: XCOL(:)
      REAL,POINTER :: ZVAL(:)
      REAL,POINTER :: F(:)
      F(1)=ZVAL(2)
      F(2)=(1-ZVAL(1)*ZVAL(1))*ZVAL(2)-4*ZVAL(1)
      END SUBROUTINE FSUB

      SUBROUTINE DFSUB (XCOL, ZVAL, DF)
      IMPLICIT NONE
C      REAL,DIMENSION(:) :: XCOL
C      REAL,DIMENSION(:) :: ZVAL
C      REAL,DIMENSION(:) :: DF
      REAL,POINTER :: XCOL(:)
      REAL,POINTER :: ZVAL(:)
      REAL,POINTER :: DF(:)
      DF(1)=0.0
      DF(2)=1.0
      DF(3)=-4.0-2*ZVAL(1)*ZVAL(2)
      DF(4)=(1.0-ZVAL(1)*ZVAL(1))
      END SUBROUTINE DFSUB

      SUBROUTINE GSUB (IZETA, ZVAL, GVAL)
      IMPLICIT NONE
      INTEGER, DIMENSION(:) :: IZETA
      REAL, DIMENSION(:) :: ZVAL
      REAL, DIMENSION(:) :: GVAL
      ! normally GVAL(1)=ZVAL(1)-0.0
      ! AND      GVAL(2)=ZVAL(2)-1.0
      ! BUT how can it know the time values?
      END SUBROUTINE GSUB


      SUBROUTINE DGSUB (IZETA, ZVAL, DG)
      IMPLICIT NONE
      INTEGER, DIMENSION(:) :: IZETA

      REAL, DIMENSION(:) :: ZVAL
      REAL, DIMENSION(:) :: DG
      END SUBROUTINE DGSUB

      SUBROUTINE GUESS (XII, ZVAL, DMVAL)
      IMPLICIT NONE
      REAL, DIMENSION(:) :: XII
      REAL, DIMENSION(:) :: ZVAL
      REAL, DIMENSION(:) :: DMVAL
      END SUBROUTINE GUESS

C  MAIN PROGRAM:
      PROGRAM VALIDVANDERPOL

        EXTERNAL FSUB, DFSUB, GSUB, DGSUB, GUESS
        EXTERNAL COLNEW
C      IMPLICIT NONE
C      VARIABLES FOR THE PROBLEM VAN DER POL
C      	X1'=X2
C       X2'=(1-X1^2)*X2-4*X1

C       GENERIC VARS
C        IMPLICIT REAL*8 (A-H,O-Z)  ! NOT ALLOWED WITH ALLOCATE !
        INTEGER,DIMENSION(:),ALLOCATABLE    :: DUMMY

        INTEGER :: K, MMAX,MSTAR
C      VARIABLES TO BE PASSED TO COLNEW:
        INTEGER :: NCOMP=2   ! no of diff. equation
        INTEGER,DIMENSION(1) :: M= (/ 2 /)
        REAL :: ALEFT=0.0
        REAL :: ARIGHT=1.0
        REAL,DIMENSION(2) :: NINTV=(/ 0.0, 1.0 /)
        INTEGER, DIMENSION(11) :: IPAR
C       IPAR(1)=1  ! non linear
C       IPAR(2)=4  ! no collocation pts per subinterval 
C       IPAR(3)=4  ! no subintervals
C       IPAR(4)=2  ! no solution and derivative tolerances
C       IPAR(5)=396 ! fspace dim see after
C       IPAR(6)=21  ! ispace dim see after
C       IPAR(7)=-1  ! output control
C       IPAR(8)=1   ! iread
C       IPAR(9)=0   ! iguess
C       IPAR(10)=0  ! pb regular
C       IPAR(11)=3  ! no points in interval alongside bounds
C       INTEGER,DIMENSION(11) :: IPAR=(/ 1, 4, 4, 2, 396, 21, -1, 1, 0, 0, 3 /)
C      4 intervals, each has 4 points..
C       21=3+((4*4)+2), 396=4+3*2+(5+4*4)*(4*4+2)+(2*2-2)*2*2 
        INTEGER,DIMENSION(1) :: LTOL=(/ 2 /)
        REAL,DIMENSION(2) :: TOL=(/ 1.0E-4, 1.0E-4 /)
        REAL,DIMENSION(3) :: FIXPNT=(/ 0.25, 0.50, 0.75 /)
        INTEGER :: ALLOCSTAT
        INTEGER :: I
C       REAL,DIMENSION(21) :: ISPACE
        INTEGER,POINTER :: ISPACE(:)
        REAL,POINTER :: FSPACE(:)
C       REAL,DIMENSION(:),ALLOCATABLE :: FSPACE
        INTEGER :: IFLAG=1

C       VARIABLES TO BE PASSED TO OTHER SUBROUTINES THAT ARE NORMALLY CALLED BY COLNEW:
        REAL,POINTER :: ZETA(:)
        INTEGER :: IZCOUNT
C       REAL,DIMENSION(:),ALLOCATABLE :: ZETA
        INTEGER        :: KSTORE, NOLD, LXIOLD, LXI, LG, NREC, LVAL
        REAL,DIMENSION(:),ALLOCATABLE :: VALSTR      
        REAL,DIMENSION(28,4)          :: ASAVE  
        REAL,DIMENSION(:),ALLOCATABLE :: XIOLD
        REAL,DIMENSION(:),ALLOCATABLE :: Z, DMZ
        REAL           :: X
C       VARIABLES FOR EVALUATION OF THE SOLUTION
        INTEGER,DIMENSION(:),ALLOCATABLE :: ISPACECPY
        REAL,DIMENSION(:),ALLOCATABLE :: FSPACECPY
        REAL     :: XVAR
        REAL,DIMENSION(:),ALLOCATABLE :: XSOL
        REAL,DIMENSION(2) :: ZVAR
        REAL,DIMENSION(:),ALLOCATABLE :: ZSOL

        ALLOCATE(ISPACE(21), ISPACECPY(21), STAT=ALLOCSTAT)
        DO 20 I=1,21
          ISPACE(I)=0
   20   CONTINUE

        ALLOCATE(FSPACE(396), FSPACECPY(396), STAT=ALLOCSTAT)
        DO 30 I=1,396
           FSPACE(I)=0.0
   30   CONTINUE

        ALLOCATE(ISPACE(21))

        IZCOUNT=4+3*M(1)+(5+IPAR(2)*IPAR(3))*(IPAR(2)*IPAR(3)+2)
     1                +(2*M(1)-2)

        ALLOCATE(ZETA(IZCOUNT), STAT=ALLOCSTAT) ! -2 at the end means 2 boundary

C        CALL SUBROUTINE COLNEW (NCOMP, M, ALEFT, ARIGHT, ZETA, IPAR, LTOL,
C      1                   TOL, FIXPNT, ISPACE, FSPACE, IFLAG,
C      2                   FSUB, DFSUB, GSUB, DGSUB, GUESS)

C      TRY THE ROUTINE APPROX: JUST TO SEE IF THE CALLS TO FSUB,DFSUB..GUESS WORKS
C
C...  save in  valstr  the values of the old solution
C...  at the relative positions 1/6 and 5/6 in each subinterval.
C
      LVAL = M(1)*4
      NMAX=2
      MSTAR=2
      NREC=0  ! seems to be zero
      LXI = 1
      LG = LXI + NMAX + 1
      LXIOLD = LG + 2*MSTAR * (NMAX * (2*MSTAR-NREC) + NREC)
      IZ = (NOLD-1) * MSTAR + 1
      IDMZ = (NOLD-1) * K * NCOMP

      ALLOCATE(VALSTR(LVAL), XIOLD(LXIOLD), Z(IZ), DMZ(IDMZ),
     1         STAT=ALLOCSTAT) 
      ! index token from loops of routines of colnew
C     THE PORTION TO TEST AS A FIRST TRY
      KSTORE = 1
      DO 130 I = 1, NOLD
          HD6 = (XIOLD(I+1) - XIOLD(I)) / 6.D0
          X = XIOLD(I) + HD6
          CALL APPROX (I, X, VALSTR(KSTORE), ASAVE(1,1), DUMMY, XIOLD,
     1         NOLD, Z, DMZ, K, NCOMP, MMAX, M, MSTAR, 4, DUMMY, 0)
          X = X + 4.D0 * HD6
          KSTORE = KSTORE  +  3 * MSTAR
          CALL APPROX (I, X, VALSTR(KSTORE), ASAVE(1,4), DUMMY, XIOLD,
     1         NOLD, Z, DMZ, K, NCOMP, MMAX, M, MSTAR, 4, DUMMY, 0)
          KSTORE = KSTORE  +  MSTAR
  130 CONTINUE
C      END OF BLOCK TRY THE ROUTINE APPROX

C       AT THE END EVALUATE THE SOLUTION USING COEFFICIENTS

        DO 40,I=1,21
          ISPACECPY(I)=ISPACE(I)
   40   CONTINUE
        DO 50,I=1,396
          FSPACECPY(I)=FSPACE(I)
   50   CONTINUE

        ALLOCATE(XSOL(51), ZSOL(102), STAT=ALLOCSTAT) ! -2 at the end means 2 boundary

        DO 60,I=0,50
          XVAR = (I)*0.02
          XSOL(I+1)=(I)*0.02
C          CALL APPSLN (XVAR, ZVAR, FSPACECPY, ISPACECPY)
          ZSOL(2*I+1)=ZVAR(1)  ! X1
          ZSOL(2*I+2)=ZVAR(2)  ! X2=DERIVATE(X1)          
  60    CONTINUE

C      DEALLOCATE MEMORY
       DEALLOCATE(VALSTR, STAT=ALLOCSTAT)
       DEALLOCATE(ISPACE,ISPACECPY,FSPACE,FSPACECPY,STAT=ALLOCSTAT)
       DEALLOCATE(ZETA,XSOL,ZSOL,STAT=ALLOCSTAT)
       DEALLOCATE(VALSTR, XIOLD, Z, DMZ, STAT=ALLOCSTAT) 

      END PROGRAM VALIDVANDERPOL


