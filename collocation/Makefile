CC=/c/mingw/bin/mingw32-g++
CANSI=/c/mingw/bin/mingw32-gcc
CXXFLAGS= -std=c++11 -Wall -O2
CANSIFLAGS= -Wall -O2
LOADER=/c/mingw/bin/mingw32-g++
LOADOPTS='-O3'

LAPACKELIB=/home/frede/lapack-3.7.1/liblapacke.a

INCLAPACKE=/home/frede/lapack-3.7.1/LAPACKE/include
TMGLIB=/lib/local/libtmglib.a
LAPACKLIB=/lib/local/liblapack.a
BLASLIB=/lib/local/librefblas.a

# std C++11 thus
CINCLUDE=../include/
#../../../Ipopt-3.12.6/include/coin/
LIBS=-lm \

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
	${LOADER} ${LOADOPTS} -o main_example -I${CINCLUDE} ${OBJ} ${LAPACKELIB}  ${LAPACKLIB} ${BLASLIB} ${TMGLIB} -lgfortran ${LIBS}

clean:
	rm -rf ../collocation/*.o ../MathFunctions/mySolveLinearLapack.o main_example 

