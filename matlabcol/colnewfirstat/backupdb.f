
        DO 20 I=1,21
          ISPACE(I)=0
   20   CONTINUE

        REAL,POINTER :: FSPACE(:)
C       REAL,DIMENSION(:),ALLOCATABLE :: FSPACE
        ALLOCATE(FSPACE(396), STAT=ALLOCSTAT)
        DO 30 I=1,396
           FSPACE(I)=0.0
   30   CONTINUE

        INTEGER :: IFLAG=1

        REAL,POINTER :: ZETA(:)
C       REAL,DIMENSION(:),ALLOCATABLE :: ZETA
        ALLOCATE(ZETA( 4+3*K(1)+(5+PAR(2)*PAR(3))*(PAR(2)*PAR(3)+2)+ &
                       (2*K(1)-2) ), STAT=ALLOCSTAT) ! -2 at the end means 2 boundary
