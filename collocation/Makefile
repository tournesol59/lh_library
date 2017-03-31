CC=/usr/bin/g++
CANSI=/usr/bin/gcc
CXXFLAGS= -std=c++11 -Wall -O2
CANSIFLAGS= -Wall -O2
#INCCLAPK=/home/frederic/CLAPACK/INCLUDE
INCLAPACKE=/home/frederic/gitrep/ThirdParty/lapack/LAPACKE/include
LOADER=/usr/bin/gcc
LOADOPTS='-O3'

TMGLIB=/home/frederic/gitrep/ThirdParty/lapack/libtmglib.a
LAPACKLIB=/home/frederic/gitrep/ThirdParty/lapack/liblapack.a
BLASLIB=/home/frederic/gitrep/ThirdParty/lapack/librefblas.a
#F2CLIB=/home/frederic/CLAPACK/F2CLIBS/libf2c.a
LAPACKELIB=/home/frederic/gitrep/ThirdParty/lapack/liblapacke.a

# std C++11 thus
CINCLUDE=../include/
#../../../Ipopt-3.12.6/include/coin/
LIBS=-lm \
     -lc

OBJ=../collocation/ident05_coll.o \
  main_example.o \
  ../MathFunctions/mySolveLinearLapack.o

main_example: main_example.cpp 
	${CC} ${CXXFLAGS} -o main_example.o -c $< -I${CINCLUDE}

ident05_coll: ../collocation/ident05_coll.cpp
	${CC} ${CXXFLAGS} -o ../collocation/ident05_coll.o -c $< -I${CINCLUDE}

mySolveLinearLapack: ../MathFunctions/mySolveLinearLapack.c
	${CANSI} ${CANSIFLAGS} -o ../MathFunctions/mySolveLinearLapack.o -c $< -I${INCLUDE} -I${INCLAPACKE}

#main no ref to std=c++11 !
main: ${OBJ}
	${LOADER} ${LOADOPTS} -o main_example -I${CINCLUDE} ${OBJ} \ 
	${LAPACKELIB}  ${LAPACKLIB} ${BLASLIB} ${TMGLIB} -lgfortran \ 
	${LIBS}

clean:
	rm -rf ../collocation/*.o ../MathFunctions/mySolveLinearLapack.o main_example 

