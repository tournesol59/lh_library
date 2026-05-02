 module wraptime
  use ISO_C_BINDING, only : C_INT
 interface 
   type(C_INT) function second( ) bind(C,name='C_wraptime')
      import C_CHAR, C_INT, C_FLOAT
   end function second
 end interface
 end module wraptime
