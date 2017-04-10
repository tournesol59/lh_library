CC=/usr/bin/g++
CANSI=/usr/bin/gcc
CXXFLAGS= -std=c++11 -Wall -O2
CANSIFLAGS= -Wall -O2
LOADER=/usr/bin/g++
LOADOPTS='-O3'

LAPACKELIB=/home/frederic/lapack-3.7.0/liblapacke.a

INCLAPACKE=/home/frederic/lapack-3.7.0/LAPACKE/include
TMGLIB=/usr/local/lib/libtmglib.a
LAPACKLIB=/usr/local/lib/liblapack.a
BLASLIB=/usr/local/lib/libblas.a

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
	${LOADER} ${LOADOPTS} -o main_example -I${CINCLUDE} ${OBJ} ${LAPACKELIB}  ${LAPACKLIB} ${BLASLIB} ${TMGLIB} -lgfortran ${LIBS}

clean:
	rm -rf ../collocation/*.o ../MathFunctions/mySolveLinearLapack.o main_example 

