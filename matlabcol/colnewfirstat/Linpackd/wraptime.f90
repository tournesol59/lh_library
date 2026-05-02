 module wraptime


  use ISO_C_BINDING, only : C_CHAR, C_INT, C_FLOAT, C_PTR
  type, bind(C) :: cel
    real(kind=C_FLOAT)  :: r
    integer(kind=C_INT) :: n
  end type cel
  interface

    !
    subroutine Cprocsec(res) bind(C,name='C_wraptime')
      import C_INT, C_FLOAT, C_PTR
      real(kind=C_FLOAT) :: res
    end subroutine Cprocsec

  end interface

 public :: second
 contains
 integer function second()
   implicit none
   real(kind=4) :: sec
   call Cprocsec(sec)
   second = sec
 end function second

 end module wraptime

PROGRAM appel_c
  use ISO_C_BINDING, only : C_CHAR, C_INT, C_FLOAT, C_PTR

  IMPLICIT NONE
  INTEGER           ::  I, K
  REAL(kind=4)      :: before, after, elapsed
  write(*,100) "starting test second:", before
  before = second()
  write(100,*) "starting test second:", before
  do i=1,100000000
     K=I
  enddo
  after = second()
  elapsed = after-before
  write(*,100) "museconds elapsed: ", elapsed
100 format(a,f10.5)
END PROGRAM appel_c
