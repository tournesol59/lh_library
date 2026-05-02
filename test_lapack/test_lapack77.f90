!----------------------------------------------------------------------
!
! test_lapack77.f90
!
! simple demo for using LAPACK77 routines in Fortran90
!
! Author: Bernhard Seiwald
! Date: 12.01.2002
!
!----------------------------------------------------------------------
PROGRAM test_lapack77
!----------------------------------------------------------------------
!----------------------------------------------------------------------
 IMPLICIT NONE

 INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
 INTEGER, PARAMETER :: DP = KIND(1.0D0)
 INTEGER, PARAMETER :: w_us = 6
 REAL(DP), PARAMETER :: EPS = 1.0D-9
 CHARACTER(len=255) :: filename, msg
 INTEGER(I4B) :: iunit, i_alloc, istat
 INTEGER(I4B) :: dim, info
 INTEGER(I4B) :: n, nrhs, lda, ldb
 INTEGER(I4B) :: i, j
 REAL(DP) :: erg
 INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: ipiv
 REAL(DP), DIMENSION(:), ALLOCATABLE :: inhom, loes
 REAL(DP), DIMENSION(:,:), ALLOCATABLE :: mat, mat1
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! open file containing info about equation system
!---------------------------------------------------------------------
 WRITE(filename, '(A)') 'test_lapack_in.dat'
 iunit = 7
 OPEN(unit=iunit, file=TRIM(filename), status='unknown', form='formatted', &
    iostat=istat)
 IF (istat /= 0) THEN
    WRITE(w_us,*) 'test_lapack77: Error on opening file ', filename
    STOP
 END IF

!----------------------------------------------------------------------
! read dimension of matrix and allocate arrays
!----------------------------------------------------------------------
 READ(iunit,*) dim
 ALLOCATE( mat(dim,dim), mat1(dim,dim), inhom(dim), loes(dim), &
   stat = i_alloc)
 IF (i_alloc /= 0) STOP 'test_lapack77: Allocation for arrays1 failed!'

!----------------------------------------------------------------------
! read matrix
!----------------------------------------------------------------------
 DO i=1, dim
    READ(iunit, *) (mat(i,j), j=1,dim)
 END DO
!----------------------------------------------------------------------
! read right hand side
!---------------------------------------------------------------------
 DO i =1, dim
    READ(iunit,*) inhom(i)
 END DO

 CLOSE(unit=iunit)

 loes = inhom
 mat1 = mat
 n = SIZE(loes,1)
 nrhs = 1
 lda = MAX(1, SIZE(mat,1))
 ldb = MAX(1, SIZE(loes,1))

 ALLOCATE( ipiv(n), stat = i_alloc)
 IF (i_alloc /= 0) STOP 'test_lapack77: Allocation for arrays2 failed'

!----------------------------------------------------------------------
! solve system
!----------------------------------------------------------------------
 CALL dgesv(n, nrhs, mat1, lda, ipiv, loes, ldb, info)
! error handling: messages are copies from man pages
 IF (info < 0) THEN
    WRITE(msg,*) 'test_lapack77: Lapack routine dgesv did not exit ', &
       'successful! The', info, '-th argument had an illegal value'
    PRINT *, msg
 END IF
 IF ( info > 0) THEN
    WRITE(msg,*) 'test_lapack77: Lapack routine dgesv did not exit ', &
      ' successful! U(', info, ',' , info, ') is exactly zero. ',  &
      'The factorization has been completed, but the factor U ',  &
      'is exactly singular, so the solution could not be computed.'
    PRINT *,msg
 END IF


!----------------------------------------------------------------------
! check result
!----------------------------------------------------------------------
  DO i = 1, dim
     erg = SUM(mat(i,:)*loes)
     IF ( ABS(inhom(i) - erg) > EPS ) THEN
 ! report error in file 'FORTRAN_INT_ERROR'
        WRITE(filename, '(A)') 'FORTRAN_INT_ERROR'
        OPEN(unit=iunit, file=TRIM(filename), status='unknown', &
           form='formatted', iostat=istat)

        IF (istat /= 0) THEN
          WRITE(w_us,*) 'test_lapack77: Error on opening file ', filename
          STOP
        END IF
        WRITE(iunit,*) 'problem in ', i, erg
        CLOSE(unit=iunit)
     ! stop, when first error occurs
        STOP
     END IF
  END DO
!----------------------------------------------------------------------
! deallocate arrays
!----------------------------------------------------------------------
DEALLOCATE( mat, mat1, inhom, loes, stat = i_alloc )
IF (i_alloc /= 0) STOP 'test_lapack77: Deallocation for arrays1 failed!'
DEALLOCATE( ipiv, stat = i_alloc )
IF (i_alloc /= 0) STOP 'test_lapack77: Deallocation for arrays2 failed!'

END PROGRAM test_lapack77

