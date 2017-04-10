CC=/usr/bin/g++
CANSI=/usr/bin/gcc
CXXFLAGS= -std=c++11 -Wall -O2
CANSIFLAGS= -Wall

#INCLAPACK=/home/frederic/Ipopt-3.12.6/ThirdParty/Lapack

LOADER=/usr/bin/g++
LOADOPTS='-O3'
#library downloaded from netlib.org:
#INCLAPACKE=/home/frederic/lapack-3.7.0/LAPACKE/include
#TMGLIB=/home/frederic/lapack-3.7.0/libtmglib.a
#LAPACKLIB=/home/frederic/lapack-3.7.0/liblapack.a
#BLASLIB=/home/frederic/lapack-3.7.0/libcblas.a
LAPACKELIB=/home/frederic/lapack-3.7.0/liblapacke.a

#library installed on linux system:
INCLAPACKE=/home/frederic/lapack-3.7.0/LAPACKE/include
TMGLIB=/usr/local/lib/libtmglib.a
LAPACKLIB=/usr/local/lib/liblapack.a
BLASLIB=/usr/local/lib/libblas.a

CINCLUDE=../include/
LIBS=-lm \
     -lc

OBJ=math_example.o \
  mySolveLinearLapack.o

math_example: math_example.cpp
	${CC} ${CXXFLAGS} -o math_example.o -c $<

mySolveLinearLapack: mySolveLinearLapack.c
	${CANSI} ${CANSIFLAGS} -o mySolveLinearLapack.o -c $< -I${INCLAPACKE}

main: ${OBJ}
#	${CANSI} ${CANSIFLAGS} -o math_example ${OBJ} ${LAPACKELIB} ${LIBS}
	${CC} ${CXXFLAGS} -o math_example ${OBJ}  ${LAPACKELIB} ${LAPACKLIB} ${BLASLIB} ${TMGLIB} -lgfortran ${LIBS}
clean:
	rm -rf mySolveLinearLapack.o math_example 

