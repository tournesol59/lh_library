--Subspace Deterministic (as opposed to stochastic) package, which
--has, among other subroutine, the "sysorder_determ" procedure, which
--computes the order of a system prior the identification of matrices A,B,C,D.
--It shall be based on singular value decomposition but simplest call to linear
--inverse of a matrix will be tested through as it a problem of mixing fortran extern procedures and Ada.

with Ada.Numerics.Generic_Elementary_Functions;
with Text_IO;
--with ada_fortran; use ada_fortran;

package body SubspaceDeterm is

  package Maths is new Ada.Numerics.Generic_Elementary_Functions(Real);-- for Sqrt
  use Maths;


  --------------
  -- Product --
  --------------
  
  --  Matrix Vector multiplication

  function Product
    (A              : in Matrix;
     X              : in Col_Vector;
     Final_Index    : in Index := Index'Last;
     Starting_Index : in Index := Index'First) 
     return Col_Vector
  is
     Result : Col_Vector := (others => Zero);
     Sum : Real;
  begin   
     for i in Starting_Index .. Final_Index loop 
        Sum := Zero; 
        for j in Starting_Index .. Final_Index loop
           Sum := Sum + A(i, j) * X(j);  
        end loop;
        Result(I) := Sum;
     end loop;

     return Result;
        
  end Product;
  
   --------------------
   --A transpose
   --------------------
   procedure A_transpose
     (A            : in out Matrix;
     A_row	   : in Index;
     A_col	   : in Index) is
     Item : Real;
   begin
	   -- first case square matrices:
      if (A_row = A_col) then
         for Row in 1..A_row loop
            for Col in Row..A_col loop
		Item :=A(Row, Col);
		A(Row, Col) :=A(Col, Row);
		A(Col, Row) := Item;
	    end loop;
	 end loop;
      end if;
      -- second case: more rows than colums
   if (A_row > A_col) then
         for Row in 1..A_row loop
            for Col in Row..A_row loop
		Item :=A(Row, Col);
		A(Row, Col) :=A(Col, Row);
		A(Col, Row) := Item;
	    end loop;
	 end loop;      
   end if;	
      -- third case: more columns than rows
   if (A_row < A_col) then
         for Row in 1..A_col loop
            for Col in Row..A_col loop
		Item :=A(Row, Col);
		A(Row, Col) :=A(Col, Row);
		A(Col, Row) := Item;
	    end loop;
	 end loop;      
   end if;   
   end A_transpose;

   ----------------------------
   -- A_transpose_times_Vec  --
   ----------------------------
  
   function A_transpose_times_Vec 
     (A            : in     Matrix;
      Vector       : in     Col_Vector;
    Final_Row    : in     Index; 
    Final_Col    : in     Index;
    Starting_Row : in     Index; 
    Starting_Col : in     Index)
      return Col_Vector is
      Result : Col_Vector := (others => Zero);
      Sum : Real; --_Extended;
   begin

      for Col in Index range Starting_Col .. Final_Col loop
         Sum := +0.0;
         for Row in Index range Starting_Row .. Final_Row loop
            Sum := Sum + (A(Row, Col)) * (Vector(Row));
         end loop;
	 Result (Col) := Real (Sum);
      end loop;

      return Result;

   end A_transpose_times_Vec;

   ------------------------------
   -- Multiplication 2 matrices with first is transposated
   ------------------------------
   function MultiplicationA_trans_X
	   (A      : in Matrix;
	    X      : in Matrix;
	    A_row  : in Index;
	    A_col  : in Index;  -- must be equal to X_row
	    X_col  : in Index)
	    return Matrix is
	    Result : Matrix :=(others => ( others => Zero));
	    Vector, ResultVec : Col_Vector;
   begin
      for Col in 1..A_col loop	  
	 for Row in 1..A_row loop
	    Vector(Row) :=X(Row, Col);
	 end loop;
	 ResultVec :=A_transpose_times_Vec(A, Vector, A_row, A_col, 1, 1);
	 for Row in 1..A_row loop
	    Result(Row, Col) :=ResultVec(Row);
	 end loop;
      end loop;

      return Result;
   end MultiplicationA_trans_X;

   -----------------
   -- A_times_Vec --
   -----------------
  
   function A_times_Vec 
     (A            : in     Matrix;
      Vector       : in     Col_Vector;
      Final_Row    : in     Index; 
      Final_Col    : in     Index;
      Starting_Row : in     Index; 
      Starting_Col : in     Index)
      return Col_Vector is
      Result : Col_Vector := (others => Zero);
      Sum : Real; --_Extended;
   begin

      for Row in Index range Starting_Row .. Final_Row loop
         Sum := +0.0;
         for Col in Index range Starting_Col .. Final_Col loop
            Sum := Sum + (A(Row, Col)) * (Vector(Col));
         end loop;
	 Result (Row) := Real (Sum);
      end loop;

      return Result;

   end A_times_Vec;

   ------------------------------
   -- Multiplication 2 matrices (classical)
   ------------------------------
   function MultiplicationXA
	   (X      : in Matrix;
	    A      : in Matrix;
	    X_row  : in Index;
	    A_row  : in Index;  -- must be equal to X_col
	    A_col  : in Index)
	    return Matrix is
	    Result : Matrix :=(others => (others => Zero));
	    Vector, ResultVec : Col_Vector;
   begin
      for Col in 1..A_col loop	  
	 for Row in 1..A_row loop
	    Vector(Row) :=A(Row, Col);
	 end loop;
	 ResultVec :=A_times_Vec(X, Vector, A_row, A_col, 1 ,1);
	 for Row in 1..A_row loop
	    Result(Row, Col) :=ResultVec(Row);
	 end loop;
      end loop;

      return Result;
   end MultiplicationXA;
 
   -- This procedure sets an Hankel matrix filled with a subset of index
   -- of data(Index) with spans over i- and j-times as shown below
   -- X(k)      X(k+1)  ..   X(k+j-1)
   -- X(k+1)    X(k+2)  ..   X(k+j)
   -- ..
   -- X(k+i-1)  X(k+i+1) ..  X(k+j+i-2)
   -- --------------------------------
   -- X(k+i)    X(k+i+2) ..  X(k+j+i-1)
   -- ..
   -- X(k+2i-1) X(k+2i)  ..  X(k+2i+j-2)
   -- --------------------------------
   --
   procedure WriteData2Hankel(
	  X : in DataVec;
	  Starting_Index : in Index;
	  FinalI_Index : in Index; -- Final - Starting must be positive even
	  FinalJ_Index : in Index;
	  HX : in out Matrix ) is
      n,m : Index;

   begin
      n := (FinalI_Index-Starting_Index)/2; --a verification here would be fine
      m := FinalJ_Index-Starting_Index;
      for i in 1..(2*n) loop
	 for j in 1..(m) loop
            HX(i,j):= X(Starting_Index+i+j-2); -- the Hankel matrix is not square
	 -- but has even rows as in prevision to the next stage
	 end loop;
      end loop;
   End WriteData2Hankel;


  procedure ConcatInOutputs
	  (U : in DataVec;
	  Y : in DataVec;
	  Starting_Index : in Index;
	  FinalI_Index : in Index; --Final - Starting must be positive even number
	  FinalJ_Index : in Index;
	  W : in out Matrix) is 
  n,m : Index;
  begin
      n := (FinalI_Index-Starting_Index); --a verification here would be fine
      m := FinalJ_Index-Starting_Index;
      for i in 1..n loop
         for j in 1..m loop
	     W(i,j):= U(Starting_Index+i-1+j-1);
	     W(n+i,j):= Y(Starting_Index+i-1+j-1);
         end loop;
      end loop;
  End ConcatInOutputs;
   -- procedure in construction previewed for calling Singular Value Decomposition
   procedure calcSVDfromFortran(
	   A : in Row_Vector;  -- type vector defined in ada_fortran module
	   A_row : in Index;
	   A_col : in Index;
--           B : Vector_Type;
--           C : Vector_Type;
	   U : in out Matrix;
	   D : in out Matrix;
	   V : in out Matrix) is
   begin
	   -- first basic test of adding two matrices with a module from fortran
--       C:=ada_fortran.add(A,B);
--	 C:=add(A,B);
   D(1,1):=1.0;
   End calcSVDfromFortran;

End SubspaceDeterm;
