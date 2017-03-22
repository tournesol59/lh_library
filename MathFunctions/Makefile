#CC=/usr/bin/g++
CANSI=/usr/bin/gcc
#CXXFLAGS= -std=c++11 -Wall -O2
CANSIFLAGS= -Wall
INCCLAPK=/home/frederic/CLAPACK/INCLUDE
#LOADER=/usr/bin/gcc
#LOADOPTS='-O3'
#ALINEX=clapk_example.o
TMGLIB=/home/frederic/CLAPACK/tmglib_LINUX.a
LAPACKLIB=/home/frederic/CLAPACK/lapack_LINUX.a
BLASLIB=/home/frederic/CLAPACK/blas_LINUX.a
F2CLIB=/home/frederic/CLAPACK/F2CLIBS/libf2c.a

CINCLUDE=../include/
LIBS=-lm \
     -lc

OBJ=math_example.o \
  mySolveLinearLapack.o

math_example: math_example.cpp 
	${CANSI} ${CANSIFLAGS} -o math_example.o -c $<

mySolveLinearLapack: mySolveLinearLapack.c
	${CANSI} ${CANSIFLAGS} -o mySolveLinearLapack.o -c $< -I${INCCLAPK}

main: ${OBJ}
	${CANSI} ${CANSIFLAGS} -o math_example ${OBJ} ${TMGLIB} ${LAPACKLIB} ${BLASLIB} ${F2CLIB} ${LIBS}

clean:
	rm -rf mySolveLinearLapack.o math_example 

