-- Big program to test the kalman filter matrix recursive programming
-- NOT COMPLETED NOW!
--
-- this program performs in four steps:
-- 1st step: charging the package Cholesky_LU from Jonathan Parker to solve the inverse of a
-- Positive Definite Symmetric matrix
-- 2nd step: read the input file "input.txt" which has the system matrices F, G, H and the 
-- covariance matrices Q, R in the following format and initialization of the problem
-- (nrow, ncol)
-- (real numbers in table format space-separated)
--
-- 3rd step: recursion algorithm, after each step i_iter, the matrix K(i_iter) of the
-- kalman filter and the matrix P(i_iter) are displayed, rows per rows on the screen. 
-- After each one, the Frobenius norm of ||K(i_iter)-K(i_iter-1)|| and
-- ||P(i_iter)-P(i_iter)|| are displayed.
--
-- 4th step: after the difference norm is lower than minimal bound or maximal number of
-- iterations i_iter_max is reached, writes the final matrices K and P in "output.txt"

with Ada.Text_IO; use Ada.Text_IO;
with cholesky_LU;
with Crout_LU;
with Ada.Numerics;
with Ada.Numerics.Generic_Elementary_Functions;
--with Ada.Float_IO;
with Ada.Integer_Text_IO;

procedure kalman_mainpgm is
---------------------------------------------
	--Step 1: charging vars and packages
---------------------------------------------
   type Real is digits 15;

   subtype Index is Integer range 1..191;
   
   subtype Row_Index is Index;
   subtype Col_Index is Index;

   Starting_Row : constant Row_Index := Index'First + 0;
   Starting_Col : constant Col_Index := Index'First + 0;
   Final_Row : constant Row_Index := Index'Last- 0;
   Final_Col : constant Col_Index := Index'Last- 0;

   type Matrix is array(Row_Index, Col_Index) of Real;
   type Row_Vect is array(Index) of Real;

   subtype Col_Vect is Row_Vect;
   --pragma Convention (Fortran, Matrix); --No! This QR prefers Ada convention.

   package math is new Ada.Numerics.Generic_Elementary_Functions (Real);
   use math;

   package ch_LU is new cholesky_LU
     (Real     => Real, 
      Index  => Index, 
      Matrix => Matrix);
   use ch_LU;

   package cr_LU is new Crout_LU
     (Real     => Real, 
      Index  => Index, 
      Matrix => Matrix);
   use cr_LU;

   package rio is new Float_IO(Real);
   use rio;

   package iio is new Integer_IO(Integer);
   use iio;

   subtype Real_Extended  is Real;     -- general case, and for best speed
 --type Real_Extended  is digits 18;   -- 18 ok on intel
   Zero : constant Real := +0.0;
   One  : constant Real := +1.0;
   Two  : constant Real := +2.0;

   --------------------
   --A transpose
   --------------------
   procedure A_transpose
     (A            : in out Matrix;
     A_row	   : in Integer;
     A_col	   : in Integer) is
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
      Vector       : in     Col_Vect)
    --Final_Row    : in     Index; 
    --Final_Col    : in     Index;
    --Starting_Row : in     Index; 
    --Starting_Col : in     Index)
      return Col_Vect is
      Result : Col_Vect := (others => Zero);
      Sum : Real_Extended;
   begin

      for Col in Col_Index range Starting_Col .. Final_Col loop
         Sum := +0.0;
         for Row in Row_Index range Starting_Row .. Final_Row loop
            Sum := Sum + Real_Extended (A(Row, Col)) * Real_Extended (Vector(Row));
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
	    A_row  : in Integer;
	    A_col  : in Integer;  -- must be equal to X_row
	    X_col  : in Integer)
	    return Matrix is
	    Result : Matrix :=(others => ( others => Zero));
	    Vector, ResultVec : Col_Vect;
   begin
      for Col in 1..A_col loop	  
	 for Row in 1..A_row loop
	    Vector(Row) :=X(Row, Col);
	 end loop;
	 ResultVec :=A_transpose_times_Vec(A, Vector);
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
      Vector       : in     Col_Vect)
    --Final_Row    : in     Index; 
    --Final_Col    : in     Index;
    --Starting_Row : in     Index; 
    --Starting_Col : in     Index)
      return Col_Vect is
      Result : Col_Vect := (others => Zero);
      Sum : Real_Extended;
   begin

      for Row in Row_Index range Starting_Row .. Final_Row loop
         Sum := +0.0;
         for Col in Col_Index range Starting_Col .. Final_Col loop
            Sum := Sum + Real_Extended (A(Row, Col)) * Real_Extended (Vector(Col));
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
	    X_row  : in Integer;
	    A_row  : in Integer;  -- must be equal to X_col
	    A_col  : in Integer)
	    return Matrix is
	    Result : Matrix :=(others => (others => Zero));
	    Vector, ResultVec : Col_Vect;
   begin
      for Col in 1..A_col loop	  
	 for Row in 1..A_row loop
	    Vector(Row) :=A(Row, Col);
	 end loop;
	 ResultVec :=A_times_Vec(X, Vector);
	 for Row in 1..A_row loop
	    Result(Row, Col) :=ResultVec(Row);
	 end loop;
      end loop;

      return Result;
   end MultiplicationXA;
  

   --------------------
   -- Frobenius_Norm --
   --------------------
   function Frobenius_Norm 
     (A : in Matrix)
      return Real
    is
      Max_A_Val : Real := Zero;
      Sum, Scaling, tmp : Real := Zero;
    begin
 
      Max_A_Val := Zero;
      for Row in Starting_Row .. Final_Row  loop
      for Col in Starting_Col .. Final_Col  loop
         if Max_A_Val < Abs A(Row, Col) then Max_A_Val := Abs A(Row, Col); end if;
      end loop;
      end loop;

      Max_A_Val := Max_A_Val + Two ** (Real'Machine_Emin + 4);
      Scaling := One / Max_A_Val;

      Sum := Zero;
      for Row in Starting_Row .. Final_Row  loop
      for Col in Starting_Col .. Final_Col  loop
         tmp := Scaling * A(Row, Col);
         Sum := Sum + tmp * tmp;
      end loop;
      end loop;

      return Sqrt (Sum) * Max_A_Val;

   end Frobenius_Norm;


  MaxName      : CONSTANT Positive := 80;
  SUBTYPE NameRange IS Positive RANGE 1..MaxName;   
  DecomposeTypName   : String(NameRange):= (others => '#');
  InFileName   : String(NameRange):= (others => '#');
  OutFileName  : String(NameRange):= (others => '#'); 
  TypNameLength : NameRange; 
  InNameLength : NameRange;
  OutNameLength: NameRange;
  InData       : Ada.Text_IO.File_Type;
  OutData      : Ada.Text_IO.File_Type;
  
  DecCholesky : Boolean;
 -- NextCh       : Character;
  Sum, Number       : Real;
  R_Index, C_Index : Integer;
  DP_norm, DPSym_norm      : Real;
  Max_DP_norm : constant Real :=0.0001;
  Col_Result_1, Col_Source_1      : Ch_LU.Col_Vector;
  Col_Result_2, Col_Source_2      : Cr_LU.Col_Vector;
  Result_1, Result_2              : Matrix;
  F_IN, G_IN, H_IN, Q_IN, R_IN    : Matrix;  --inputs
  Kal_Iter, P_Iter, P_Next, DP_Iter, DPSym  : Matrix;  -- outputs
  Eps_cov, Eps_cov_tr, Sig_res, Sig_res_tr  : Matrix;  -- intermediate
  Eps_cov_tr_diag_inv             : Ch_LU.Col_Vector;  -- internal data for cholesky
  Eps_cov_tr_scale_factors	  : Cr_LU.Scale_Vectors;  -- internal data for crout
  Eps_cov_tr_rearrange		  : Cr_LU.Rearrangement;
  P_interm, F_KH, F_KH_tr, GQG, KRK         : Matrix;  -- intermediate
  n_test, n_test2, n_n, n_r, n_p  : Integer; -- the dimensions
  i_iter                          : Integer; -- useful for the algorithm

  
    
begin
  ---------------------------------------------
  -- Step 2: read input file and initialization
  ---------------------------------------------		

   -- get input file name and open it
  Ada.Text_IO.Put(Item => "Please enter name of file to copy >");
  Ada.Text_IO.Get_Line(Item => InFileName, Last => InNameLength);
  Ada.Text_IO.Open(File => InData, Mode => Ada.Text_IO.In_File, 
		   Name => InFileName(1..InNameLength));
   -- get output file name and create it
  Ada.Text_IO.Put(Item => "Please enter name of the new file >");
  Ada.Text_IO.Get_Line(Item => OutFileName, Last => OutNameLength);
  Ada.Text_IO.Create(File=>OutData, Mode=>Ada.Text_IO.Out_File, 
                     Name=>OutFileName(1..OutNameLength));
  -- get cholesky or crout:
  Ada.Text_IO.Put(Item => "Please enter 'ch' for cholesky or 'cr' for crout >");
  Ada.Text_IO.Get_Line(Item => DecomposeTypName, Last => TypNameLength);


  R_Index :=1;
  C_Index :=1;
  -- Ada.Integer_IO.Get(File => InData, Item => n_p, Width => 0);  -- Linux
  Get(File => InData, Item => n_n, Width => 0); -- Windows
  -- Ada.Integer_IO.Get(File => InData, Item => n_p, Width => 0); --Linux, A_IN is square mat
  Get(File => InData, Item => n_test, Width => 0); --Windows, F_IN is square mat
  If (n_test /= n_n) then
	Ada.Text_IO.Put("Error, F_IN must be a square matrix.");
	Ada.Text_IO.New_Line;
--        exit;  -- exit not possible, try goto in future
  end if;
 
  loop
     exit when (R_Index > n_n);  
     C_Index :=1;
     loop
--	exit when Ada.Text_IO.End_Of_Line(InData);
	exit when (C_Index > n_n);
--	Ada.Float_IO.Get(File => InData, Item => Number, Width => 0);	-- Linux
	Get(File => InData, Item => Number, Width => 0);  -- Windows

	F_IN(R_Index, C_Index) := Number;	
	C_Index :=C_Index+1;
     end loop;
     R_Index :=R_Index+1;
  end loop;

  Get(File => InData, Item => n_test, Width => 0);  -- Windows
  -- Ada.Integer_IO.Get(File => InData, Item => n_test, Width => 0);  -- Linux  
  If (n_test /= n_n) then
	Ada.Text_IO.Put("Error, F_IN and G_IN must have same number of rows.");
	Ada.Text_IO.New_Line;
--        exit;  -- exit not possible, try goto in future
  end if;  
  Get(File => InData, Item => n_r, Width => 0);  -- Windows
  R_Index :=1;
  C_Index :=1;
  loop
	  exit when (R_Index > n_n);
	  C_Index :=1;
	  loop
		  exit when (C_Index > n_r);
		  Get(File => InData, Item => Number, Width => 0); --Windows
		  G_IN(R_Index, C_Index) :=Number;
		  C_Index :=C_Index+1;
	  end loop;
	  R_Index :=R_Index+1;
  end loop;

  Get(File => InData, Item => n_p, Width => 0);  -- Windows  
  Get(File => InData, Item => n_test, Width => 0);  -- Windows
  -- Ada.Integer_IO.Get(File => InData, Item => n_test, Width => 0);  -- Linux  
  If (n_test /= n_n) then
	Ada.Text_IO.Put("Error, F_IN and H_IN must have same number of columns.");
	Ada.Text_IO.New_Line;
--        exit;  -- exit not possible, try goto in future
  end if;  
  R_Index :=1;
  C_Index :=1;
  loop
	  exit when (R_Index > n_p);
	  C_Index :=1;
	  loop
		  exit when (C_Index > n_n);
		  Get(File => InData, Item => Number, Width => 0); --Windows
		  H_IN(R_Index, C_Index) :=Number;
		  C_Index :=C_Index+1;
	  end loop;
	  R_Index :=R_Index+1;
  end loop;

  Get(File => InData, Item => n_test, Width => 0);  -- Windows  
  Get(File => InData, Item => n_test2, Width => 0);  -- Windows
  -- Ada.Integer_IO.Get(File => InData, Item => n_test, Width => 0);  -- Linux  
  If (n_test /= n_p) then
	Ada.Text_IO.Put("Error, H_IN and Q_IN must have same number of rows.");
	Ada.Text_IO.New_Line;
--        exit;  -- exit not possible, try goto in future
  end if;  
  If (n_test /= n_test2) then
	Ada.Text_IO.Put("Error, Q_IN must be square matrix.");
	Ada.Text_IO.New_Line;
--        exit;  -- exit not possible, try goto in future
  end if;
  R_Index :=1;
  C_Index :=1;
  loop
	  exit when (R_Index > n_p);
	  C_Index :=1;
	  loop
		  exit when (C_Index > n_p);
		  Get(File => InData, Item => Number, Width => 0); --Windows
		  Q_IN(R_Index, C_Index) :=Number;
		  C_Index :=C_Index+1;
	  end loop;
	  R_Index :=R_Index+1;
  end loop;

  Get(File => InData, Item => n_test, Width => 0);  -- Windows  
  Get(File => InData, Item => n_test2, Width => 0);  -- Windows
  -- Ada.Integer_IO.Get(File => InData, Item => n_test, Width => 0);  -- Linux  
  If (n_test /= n_r) then
	Ada.Text_IO.Put("Error, G_IN and R_IN must have same number of columns.");
	Ada.Text_IO.New_Line;
--        exit;  -- exit not possible, try goto in future
  end if;  
  If (n_test /= n_test2) then
	Ada.Text_IO.Put("Error, R_IN must be square matrix.");
	Ada.Text_IO.New_Line;
--        exit;  -- exit not possible, try goto in future
  end if;
  R_Index :=1;
  C_Index :=1;
  loop
	  exit when (R_Index > n_r);
	  C_Index :=1;
	  loop
		  exit when (C_Index > n_r);
		  Get(File => InData, Item => Number, Width => 0); --Windows
		  R_IN(R_Index, C_Index) :=Number;
		  C_Index :=C_Index+1;
	  end loop;
	  R_Index :=R_Index+1;
  end loop;

  -- ! close Input file
  Ada.Text_IO.Close(File => InData);

  for R_Index in 1..n_n loop
     for C_Index in 1..n_n loop
	     P_iter(R_Index, C_Index) :=Zero;
     end loop;
     P_iter(R_index, R_Index) :=One;  -- identity matrix
  end loop;

  ---------------------------------
  -- Step 3: iterative algorithm
  ---------------------------------
   DP_norm :=2.0 * Max_DP_norm;
   DPSym_norm := 0.0;
   i_iter :=1;
   Kal_iter :=(others => (others => Zero));
   P_Next :=(others => (others => Zero));
   for Row in 1..n_n loop
      P_Next(Row, Row) := One;  -- Initialization has to be positive-definite
   end loop;
   loop
	   exit when (DP_norm < Max_DP_norm) or (i_iter> 1000);

  Put(OutData, "Number of iterations: ");
  Put(File=> OutData, 
                   Item => i_iter);
  Ada.Text_IO.New_Line(OutData);
  Ada.Text_IO.Put(OutData, "Norm of the Difference of the Matrix P(i+1) and P(i): ");
  Put(File=> OutData, 
  		   Item => DP_norm,
                   Fore => Default_Fore,
		   Aft => Default_Aft,
		   Exp => Default_Exp);
  Ada.Text_IO.New_Line(OutData);
  Ada.Text_IO.Put(OutData, "degree of symmetricness of matrix P(i+1)/ ||P-P'||/||P||: ");
  Put(File=> OutData, 
  		   Item => DPSym_norm,
                   Fore => Default_Fore,
		   Aft => Default_Aft,
		   Exp => Default_Exp);
  Ada.Text_IO.New_Line(OutData);
-- Writes Kal_Iter in Output File
  Put(File => OutData, Item => n_n);  -- Windows
  Put(File => OutData, Item => n_p);  -- Windows
  Ada.Text_IO.New_Line(OutData);
  R_Index :=1;
  C_Index :=1;
  loop 
	exit when (R_Index > n_n);
	C_Index :=1;
	loop 
		exit when (C_Index > n_p);
		Number :=Kal_Iter(R_Index, C_Index);
--		Ada.Float_IO.Put(File => OutData, Item => Number,
--			 Fore => Default_Fore,			 
--			 Aft => Default_Aft,			 
--			 Exp => Default_Exp);  -- Linux
		Put(File => OutData, Item => Number,
			 Fore => Default_Fore,			 
			 Aft => Default_Aft,			 
			 Exp => Default_Exp);  -- Windows
		C_Index :=C_Index+1;
	end loop;
	Ada.Text_IO.New_Line(OutData);
	R_Index :=R_Index+1;
  end loop;
-- Writes P_Next in Output File
  Put(File => OutData, Item => n_n);  -- Windows
  Put(File => OutData, Item => n_n);  -- Windows
  Ada.Text_IO.New_Line(OutData);
  R_Index :=1;
  C_Index :=1;
  loop 
	exit when (R_Index > n_n);
	C_Index :=1;
	loop 
		exit when (C_Index > n_n);
		Number :=P_Next(R_Index, C_Index);
--		Ada.Float_IO.Put(File => OutData, Item => Number,
--			 Fore => Default_Fore,			 
--			 Aft => Default_Aft,			 
--			 Exp => Default_Exp);  -- Linux
		Put(File => OutData, Item => Number,
			 Fore => Default_Fore,			 
			 Aft => Default_Aft,			 
			 Exp => Default_Exp);  -- Windows
		C_Index :=C_Index+1;
	end loop;
	Ada.Text_IO.New_Line(OutData);
	R_Index :=R_Index+1;
  end loop;	   

	   -- Calculation of Matrix Kalman Gain Ki = Fi*Pi*Hi'*(Ri+Hi*Pi*Hi')^-1
	   -- by the following way:
	   -- first, calculate Epsilon,i = Ri+Hi*Pi*Hi' and Sigma,i=Fi*Pi*Hi'
	   -- second, since there is a post-multiplication by the inverse of Epsilon,i
	   -- tranpose(Ki) is the solution of the linear system:
	   --  trans(Epsilon,i) * trans(Ki) = trans(Sigma,i)
	   -- which can be solved by using the package cholesky_LU
	   -- Note(1): the indices i refers in the program to i_iter (iteration i)
	   --          and apply only for matrices Ki and Pi (Fi,Hi,Ri,Gi,Qi) are constant
	   -- Note(2): Ki is a Matrix with n_n rows and n_p colums
	   -- Note(3): Pi is a square symmetric, positive matrix with dimension n_n
	   for Row in 1..n_n loop
	      for Col in 1..n_n loop
		 P_Iter(Row, Col) :=P_Next(Row, Col);
	      end loop;
	   end loop;
	   Eps_cov :=(others => (others => Zero));
	   Result_1 := MultiplicationXA(X => H_IN,
	  				A => P_iter,
					X_row => n_p,
					A_row => n_n,
					A_col => n_n);   -- Result_1 := H_IN*P_Iter
	   for Row in 1..n_p loop
	      for Col in 1..n_n loop
		 Result_2(Row, Col) :=H_IN(Row, Col);
	      end loop;
	   end loop;
	   A_transpose(Result_2, n_p, n_n);              -- Result_2 := H_IN'
	   Eps_cov := MultiplicationXA(X => Result_1,
	  				A => Result_2,
					X_row => n_p,
					A_row => n_n,
					A_col => n_p);   -- Eps_cov := H_IN*P_Iter*H_IN'
	   for Row in 1..n_p loop
	      for Col in 1..n_p loop
	   	 Eps_cov_tr(Row, Col) := Eps_cov(Row, Col);  -- Eps_cov_tr := H_IN'*P_Iter'*H_IN
	      end loop;
	   end loop;
	   A_transpose(Eps_cov_tr, n_p, n_p);

	   for Row in 1..n_p loop
	      for Col in 1..n_p loop
		 Eps_cov_tr(Row, Col) := Eps_cov_tr(Row, Col) + R_IN(Col, Row); -- transpose R_IN here
	      end loop;
	   end loop;   -- Eps_cov_tr := R_IN' + H_IN'*P_Iter'*H_IN
	   
	   Result_1 := MultiplicationXA(X => F_IN,
	   				A => P_Iter,
					X_row => n_n,
					A_row => n_n,
					A_col => n_n);  -- Result_1 := F_IN*P_Iter
	   for Row in 1..n_p loop
	      for Col in 1..n_n loop
		 Result_2(Row, Col) :=H_IN(Row, Col);
	      end loop;
	   end loop;
	   A_transpose(Result_2, n_p, n_n);             -- Result_2 := H_IN'
	   Sig_res := MultiplicationXA(X => Result_1,
	   				A => Result_2,
					X_row => n_n,
					A_row => n_n,
					A_col => n_p);  -- Sig_res := F_IN*P_Iter*H_IN'
	   for Row in 1..n_n loop
	       for Col in 1..n_p loop
		  Sig_res_tr(Row, Col) := Sig_res(Row, Col);
	       end loop;
	   end loop;
	   A_transpose(Sig_res_tr, n_n, n_p);           -- Sig_res_tr := H_IN*P_Iter*F_IN'


-- Writes Sig_res in Output File
  Put(File => OutData, Item => n_n);  -- Windows
  Put(File => OutData, Item => n_p);  -- Windows
  Ada.Text_IO.New_Line(OutData);
  R_Index :=1;
  C_Index :=1;
  loop 
	exit when (R_Index > n_n);
	C_Index :=1;
	loop 
		exit when (C_Index > n_p);
		Number :=Sig_res(R_Index, C_Index);
		Put(File => OutData, Item => Number,
			 Fore => Default_Fore,			 
			 Aft => Default_Aft,			 
			 Exp => Default_Exp);  -- Windows
		C_Index :=C_Index+1;
	end loop;
	Ada.Text_IO.New_Line(OutData);
	R_Index :=R_Index+1;
  end loop;
           -- selected Decomposition Typ:
           If String'Image(DecomposeTypName(2)) = "h" then
	      Ch_LU.LU_Decompose(A => Eps_cov_tr,
                               Diag_Inverse => Eps_cov_tr_diag_inv,
                              Final_Index => n_p, 
			      Starting_Index => 1);
   --                         Final_Index => Index'Last,
   --                         Starting_Index => Index'First);			      
	                      -- Eps_cov_tr is over-written with LU. 

	   	for Col in 1..n_n loop
	      	   for Row in 1..n_p loop
	         	Col_Source_1(Row) := Sig_res_tr(Row, Col);
		 	Col_Result_1(Row) := Zero;
	      	   end loop;                       -- Col_Source_1(:) := Sig_res_tr(:, Col)   (n_p x 1)

	           Ch_LU.Solve (X => Col_Result_1,
                            B => Col_Source_1,
                            A_LU => Eps_cov_tr,
			    Diag_Inverse => Eps_cov_tr_diag_inv, 
			    Final_Index => n_p,
			    Starting_Index => 1);
   --                         Final_Index => Index'Last,
   --                         Starting_Index => Index'First);
	            for Row in 1..n_p loop
	          	 Kal_Iter(Col, Row) := Col_Result_1(Row); 
		  	-- Important: trans(Ki) has been solved => transpose the result
	      	    end loop;
	        end loop;
	   else if String'Image(DecomposeTypName(2)) = "r" then
		cr_LU.LU_Decompose(A => Eps_cov_tr,    -- A is overwritten with L and U
      				Scalings => Eps_cov_tr_scale_factors,
      				Row_Permutation => Eps_cov_tr_rearrange,
      				Final_Index => n_p,
      				Starting_Index => 1,
      				Scaling_Desired => False);
	   	for Col in 1..n_n loop
	      	   for Row in 1..n_p loop
	         	Col_Source_2(Row) := Sig_res_tr(Row, Col);
		 	Col_Result_2(Row) := Zero;
	      	   end loop;                       -- Col_Source_1(:) := Sig_res_tr(:, Col)   (n_p x 1)
		      cr_LU.LU_Solve(X => Col_Result_2,
      				B => Col_Source_2,
      				A_LU => Eps_cov_tr,
      				Scalings => Eps_cov_tr_scale_factors,
      				Row_Permutation => Eps_cov_tr_rearrange,
      				Final_Index => n_p,
      				Starting_Index => 1);
       
	      	    for Row in 1..n_p loop
	          	 Kal_Iter(Col, Row) := Col_Result_2(Row); 
		  	-- Important: trans(Ki) has been solved => transpose the result
	      	    end loop;
	        end loop;
	      else 
		Ada.Text_IO.Put("unrecognizable character");
		Ada.Text_IO.New_Line;
	      end if;
	   end if;	-- Now the matrix Kal_Iter has been obtained by one or either method


	   -- Calculation of Matrix Norm Covariance Pi+1 =(Fi-Ki*Hi)*Pi*(Fi-Ki*Hi)' +Gi*Qi*Gi' + Ki*Ri*Ki'
	   -- no complicated operations (MultiplicationXA and transpose)
	   --
	   for Row in 1..n_n loop
	       for Col in 1..n_n loop
	           F_KH(Row, Col) := Zero;
	       end loop;
	   end loop;
	   F_KH := MultiplicationXA( X => Kal_Iter,
	                             A => H_IN,
				     X_row => n_n,
				     A_row => n_p,
				     A_col => n_n);      -- F_KH := Kal_Iter*H_IN
	  
	   for Row in 1..n_n loop
	       for Col in 1..n_n loop
	           F_KH(Row, Col) := F_IN(Row, Col) - F_KH(Row, Col);
		   F_KH_tr(Col, Row) := F_KH(Row, Col);  -- transpose
	       end loop;
	   end loop;                                     -- F_KH := F_IN - Kal_Iter*H_IN
	   P_interm := MultiplicationXA( X => F_KH,
	                                 A => P_Iter,
					 X_row => n_n,
					 A_row => n_n,
					 A_col => n_n);  -- P_iterm := F_KH*P_Iter
	   P_Next := MultiplicationXA( X => P_interm,
	                               A => F_KH_tr,
				       X_row => n_n,
				       A_row => n_n,
				       A_col => n_n);	 -- P_Next := F_KH*P_Iter*F_KH'
	   for Row in 1..n_n loop
	       for Col in 1..n_p loop
	           Result_1(Row, Col) := Zero;
		   Result_2(Col, Row) := G_IN(Row, Col);  -- transpose
	       end loop;
	   end loop;
	   Result_1 := MultiplicationXA( X => G_IN,
	                                 A => Q_IN,
					 X_row => n_n,
					 A_row => n_p,
					 A_col => n_p);   -- Result_1 := G_IN*Q_IN
	   GQG := MultiplicationXA( X => Result_1,
	                            A => Result_2,
				    X_row => n_n,
				    A_row => n_p,
				    A_col => n_n);        -- QGQ := G_IN*Q_IN*G_IN'
	   for Row in 1..n_n loop
	       for Col in 1..n_p loop
	           Result_1(Row, Col) := Zero;
		   Result_2(Col, Row) := Kal_Iter(Row, Col);  -- transpose
	       end loop;
	   end loop;
	   Result_1 := MultiplicationXA( X => Kal_Iter,
	                                 A => R_IN,
					 X_row => n_n,
					 A_row => n_p,
					 A_col => n_p);     -- Result_1 := Kal_Iter*R_IN
	   KRK := MultiplicationXA( X => Result_1,
	                            A => Result_2,
				    X_row => n_n,
				    A_row => n_p,
				    A_col => n_n);         -- Result_2 := Kal_Iter*R_IN*Kal_Iter'
	   for Row in 1..n_n loop
	       for Col in 1..n_n loop
		   P_Next(Row, Col) := P_Next(Row, Col) + GQG(Row, Col) + KRK(Row, Col);
		   DP_Iter(Row, Col) := P_Next(Row, Col) - P_Iter(Row, Col);
	       end loop;
	   end loop;
	   for Row in 1..n_n loop
	       for Col in 1..n_n loop
		   DPSym(Row, Col) := P_Next(Row, Col) - P_Next(Col, Row);
	       end loop;
	   end loop;	       
	   DP_norm := Frobenius_Norm(DP_Iter);  -- calculate the norm of the difference of the matrix Pi+1 and Pi between each iteration to determine the exit of the loop
	   DPSym_norm := Frobenius_Norm(DPSym) / Frobenius_Norm(P_Next);  -- evaluate the degree of symmetricness of the Matrix P_Next
	   i_iter := i_iter+1;                  -- or max number of iterations ...
	     Ada.Text_IO.Put("Number of iterations: ");
 		Put(Item => i_iter);
	  	Ada.Text_IO.New_Line;
   end loop;  -- i_iter


  ---------------------------------
  --Step 4: write results
  ---------------------------------
-- Writes Kal_Iter in output file (test only)
  --
  Put(File => OutData, Item => n_n);  -- Windows
  Put(File => OutData, Item => n_p);  -- Windows
  Ada.Text_IO.New_Line(OutData);
  R_Index :=1;
  C_Index :=1;
  loop 
	exit when (R_Index > n_n);
	C_Index :=1;
	loop 
		exit when (C_Index > n_p);
		Number :=Kal_Iter(R_Index, C_Index);
--		Ada.Float_IO.Put(File => OutData, Item => Number,
--			 Fore => Default_Fore,			 
--			 Aft => Default_Aft,			 
--			 Exp => Default_Exp);  -- Linux
		Put(File => OutData, Item => Number,
			 Fore => Default_Fore,			 
			 Aft => Default_Aft,			 
			 Exp => Default_Exp);  -- Windows
		C_Index :=C_Index+1;
	end loop;
	Ada.Text_IO.New_Line(OutData);
	R_Index :=R_Index+1;
  end loop;
-- Writes P_Next in Output File
  Put(File => OutData, Item => n_n);  -- Windows
  Put(File => OutData, Item => n_n);  -- Windows
  Ada.Text_IO.New_Line(OutData);
  R_Index :=1;
  C_Index :=1;
  loop 
	exit when (R_Index > n_n);
	C_Index :=1;
	loop 
		exit when (C_Index > n_n);
		Number :=P_Next(R_Index, C_Index);
--		Ada.Float_IO.Put(File => OutData, Item => Number,
--			 Fore => Default_Fore,			 
--			 Aft => Default_Aft,			 
--			 Exp => Default_Exp);  -- Linux
		Put(File => OutData, Item => Number,
			 Fore => Default_Fore,			 
			 Aft => Default_Aft,			 
			 Exp => Default_Exp);  -- Windows
		C_Index :=C_Index+1;
	end loop;
	Ada.Text_IO.New_Line(OutData);
	R_Index :=R_Index+1;
  end loop;
  Ada.Text_IO.Close(File => OutData);
  
  Ada.Text_IO.Put("Number of iterations: ");
  Put(Item => i_iter);
  Ada.Text_IO.New_Line;
  Ada.Text_IO.Put("Norm of the Difference of the Matrix P(i+1) and P(i): ");
  Put(Item => DP_norm,
                   Fore => Default_Fore,
		   Aft => Default_Aft,
		   Exp => Default_Exp);
  Ada.Text_IO.New_Line;


end kalman_mainpgm;

  
