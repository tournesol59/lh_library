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

#define PI 3.14157
#define PI4 0.7854
      
      SUBROUTINE FSUB (XCOL, ZVAL, F)
      IMPLICIT NONE
      REAL         :: XCOL(:)
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
      REAL         :: XCOL
      REAL,POINTER :: ZVAL(:)
      REAL,POINTER :: DF(:)
      DF(1)=0.0
      DF(2)=1.0
      DF(3)=-4.0-2*ZVAL(1)*ZVAL(2)
      DF(4)=(1.0-ZVAL(1)*ZVAL(1))
      END SUBROUTINE DFSUB

      SUBROUTINE GSUB (IZETA, ZVAL, GVAL)
      IMPLICIT NONE
      INTEGER            :: IZETA
      REAL, DIMENSION(:) :: ZVAL
      REAL               :: GVAL
      ! normally GVAL(1)=ZVAL(1)-0.0
      ! AND      GVAL(2)=ZVAL(2)-1.0
      ! BUT how can it know the time values?
      ! Re: IT SEEMS FROM LINE 2081.. FROM COLNEW.f THAT IZETA
      ! IS AN INDEX OF THE BC AND ZVAL IS THE CORRESP VECTOT AT TIME
      ! OF IZETA, GVAL IS A SCALAR
      IF (IZETA.eq.1) THEN
         GVAL=ZVAL(1)-0.0 ! left
      ELSE IF (IZETA.eq.2) THEN ! right
         GVAL=ZVAL(1)-0.897
      ENDIF
      END SUBROUTINE GSUB


      SUBROUTINE DGSUB (IZETA, ZVAL, DG)
      IMPLICIT NONE
      INTEGER            :: IZETA
      REAL,POINTER       :: ZVAL(:)
      REAL,POINTER      :: DG(:)
      DG(1)=1.0  ! VALID FOR BOTH LEFT AND RIGHT BC CONDITIONS
      DG(2)=1.0
      END SUBROUTINE DGSUB

      SUBROUTINE GUESS (XII, ZVAL, DMVAL)
      IMPLICIT NONE
      REAL               :: XII
      REAL,POINTER       :: ZVAL(:)
      REAL,POINTER       :: DMVAL(:)
C   XII IS AN INPUT: IE A SCALAR, BUT IT IS USEFUL TO REMIND THAT
C     THE INITIAL MESH FURNISHED HAS SIZE NMAX*(K)+1
C   ZVAL HAS MSTAR COMPONENTS AND DMVAL HAS NCOMP (ALL Mjth DERIVs)
C   COMPONENTS
C      INTEGER            :: I
         ZVAL(1)=SIN(2*XII)  ! A REASONABLE APPROX TO VDPOL(1,2*2)
                             ! BY A SINE SOLUTION
         ZVAL(2)=2*COS(2*XII)

         DMVAL(1)=-4*SIN(2*XII)
      END SUBROUTINE GUESS

C  MAIN PROGRAM:
      PROGRAM VALIDVANDERPOL

        IMPLICIT NONE
        EXTERNAL FSUB, DFSUB, GSUB, DGSUB, GUESS
        EXTERNAL COLNEW
C      VARIABLES FOR THE PROBLEM VAN DER POL
C      	X1'=X2
C       X2'=(1-X1^2)*X2-4*X1

C       GENERIC VARS
C        IMPLICIT REAL*8 (A-H,O-Z)  ! NOT ALLOWED WITH ALLOCATE !
        INTEGER,DIMENSION(:),ALLOCATABLE    :: DUMMY

        INTEGER :: K,KD,NMAX,MMAX,MSTAR,NREC
        INTEGER :: NCOMP
        INTEGER,DIMENSION(1) :: M= (/ 2 /)
        REAL :: ALEFT
        REAL :: ARIGHT
        REAL,DIMENSION(2) :: NINTV
        INTEGER, DIMENSION(11) :: IPAR
        INTEGER IPAR5,IPAR6
        INTEGER,DIMENSION(1) :: LTOL
        REAL,DIMENSION(2) :: TOL
        REAL,DIMENSION(3) :: FIXPNT
        INTEGER :: ALLOCSTAT
        INTEGER :: I
C       REAL,DIMENSION(21) :: ISPACE
        INTEGER,POINTER :: ISPACE(:)
        REAL,POINTER :: FSPACE(:)
C       REAL,DIMENSION(:),ALLOCATABLE :: FSPACE
        INTEGER :: IFLAG

C       VARIABLES TO BE PASSED TO OTHER SUBROUTINES THAT ARE NORMALLY CALLED BY COLNEW:
        REAL,POINTER :: ZETA(:)
        INTEGER :: IZCOUNT
C       REAL,DIMENSION(:),ALLOCATABLE :: ZETA
        INTEGER        :: KSTORE, NOLD, LVAL
        REAL,DIMENSION(:),ALLOCATABLE :: VALSTR      
        REAL,DIMENSION(28,4)          :: ASAVE  
        REAL,DIMENSION(:),ALLOCATABLE :: XIOLD
        REAL,DIMENSION(:),ALLOCATABLE :: Z, DMZ
        REAL           :: X
C       VARIABLES FOR EVALUATION OF THE SOLUTION
        INTEGER,DIMENSION(:),ALLOCATABLE :: ISPACECPY
        REAL,DIMENSION(:),ALLOCATABLE :: FSPACECPY
        REAL     :: XVAR, HD6
        REAL,DIMENSION(:),ALLOCATABLE :: XSOL
        REAL,DIMENSION(2) :: ZVAR
        REAL,DIMENSION(:),ALLOCATABLE :: ZSOL
        INTEGER    ::  LXI, LG, LXIOLD, LW, LV, LZ, LDMZ, LDELZ 
        INTEGER    ::  LDELDZ, LDQZ, LDQDMZ, LRHS, LVALST, LSLOPE
        INTEGER    ::  LACCUM, LSCL, LDSCL, LPVTG, LPVTW, LINTEG 
C 
        IFLAG=-1
        K=4
        KD=K*1
        NMAX=10
        MSTAR=2
        NREC=2
C      VARIABLES TO BE PASSED TO COLNEW:
        NCOMP=2   ! no of diff. equation
        ALEFT=0.0
        ARIGHT=2*PI/2/4
        NINTV(1)= 0.0
        NINTV(2)= PI4
      IPAR5=K+3*MSTAR+(5+K*NCOMP)*(K*NCOMP+MSTAR)+(2*MSTAR-NREC)*2*MSTAR
      IPAR5=NMAX*IPAR5
      IPAR6=3+(K*NCOMP+MSTAR)
      IPAR6=IPAR6*NMAX
       IPAR(1)=1  ! non linear
       IPAR(2)=K  ! 4 no collocation pts per subinterval 
       IPAR(3)=NMAX  ! 4 no subintervals
       IPAR(4)=2  ! no solution and derivative tolerances
       IPAR(5)=IPAR5  ! 720 ispace dim see after
       IPAR(6)=IPAR6  ! 90 ispace dim see after
       IPAR(7)=-1  ! output control
       IPAR(8)=1   ! iread
       IPAR(9)=0   ! iguess
       IPAR(10)=0  ! pb regular
       IPAR(11)=3  ! no points in interval alongside bounds
C      IPAR=(/ 1,4,4,2,720,90,-1,1,1,0,3 /)
C      4 intervals, each has 4 points..
C       90=NMAX*sizei wh sizei:=3+((4*1)+2), 720=NMAX*sizef wh. sizef:=72=4+3*2+(5+4*1)*(4*1+2)+(2*2-1)*2*2 
        LTOL(1)= 2 
        TOL(1)= 1.0E-4
        TOL(2)= 1.0E-4 
        FIXPNT(1)= 0.25
        FIXPNT(2)= 0.50
        FIXPNT(3)= 0.75
       ALLOCATE(ISPACE(IPAR6), ISPACECPY(IPAR6), STAT=ALLOCSTAT)
        DO 20 I=1,IPAR6
          ISPACE(I)=0
   20   CONTINUE

        ALLOCATE(FSPACE(IPAR5), FSPACECPY(IPAR5), STAT=ALLOCSTAT)
        DO 30 I=1,IPAR5
           FSPACE(I)=0.0
   30   CONTINUE
C   From documentation, the FSPACE(ith) first numbers are shall be equal
C   from the x(n+1) points of the initial mesh (with 1=ALEFT and
C   n+1=ARIGHT   
        FSPACE(1)=ALEFT        
        DO 31 I=1,NMAX-1
           FSPACE(I+1)=I*(ARIGHT-ALEFT)/NMAX
   31   CONTINUE
        FSPACE(NMAX+1)=ARIGHT

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
      LW = LXIOLD + NMAX + 1
      LV = LW + KD**2 * NMAX
      LZ = LV + MSTAR * KD * NMAX     !(NOLD-1) * MSTAR + 1
      LDMZ = LZ + MSTAR *(NMAX + 1)  !(NOLD-1) * K * NCOMP
      LDELZ = LDMZ + KD * NMAX
      LDELDZ = LDELZ + MSTAR * (NMAX + 1)
      LDQZ = LDELDZ + KD * NMAX
      LDQDMZ = LDQZ + MSTAR * (NMAX+ 1)
      LRHS = LDQDMZ + KD * NMAX
      LVALST = LRHS + KD * NMAX + MSTAR
      LSLOPE = LVALST + 4 * MSTAR * NMAX
      LACCUM = LSLOPE + NMAX
      LSCL = LACCUM + NMAX + 1
      LDSCL = LSCL + MSTAR * (NMAX + 1)
      LPVTG = 1
      LPVTW = LPVTG + MSTAR * (NMAX + 1)
      LINTEG = LPVTW + KD * NMAX

      ALLOCATE(VALSTR(LVAL), XIOLD(LXIOLD), Z(LZ), DMZ(LDMZ),
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


