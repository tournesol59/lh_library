#CC=/usr/bin/g++
CANSI=/usr/bin/gcc
#CXXFLAGS= -std=c++11 -Wall -O2
CANSIFLAGS= -Wall

INCLAPACKE=/home/frederic/gitrep/ThirdParty/lapack/LAPACKE/include
LOADER=/usr/bin/g++
LOADOPTS='-O3'

TMGLIB=/home/frederic/gitrep/ThirdParty/lapack/libtmglib.a
LAPACKLIB=/home/frederic/gitrep/ThirdParty/lapack/liblapack.a
BLASLIB=/home/frederic/gitrep/ThirdParty/lapack/librefblas.a
#F2CLIB=/home/frederic/CLAPACK/F2CLIBS/libf2c.a
LAPACKELIB=/home/frederic/gitrep/ThirdParty/lapack/liblapacke.a
CINCLUDE=../include/
LIBS=-lm \
     -lc

OBJ=math_example.o \
  mySolveLinearLapack.o

math_example: math_example.c
	${CANSI} ${CANSIFLAGS} -o math_example.o -c $<

mySolveLinearLapack: mySolveLinearLapack.c
	${CANSI} ${CANSIFLAGS} -o mySolveLinearLapack.o -c $< -I${INCLAPACKE}

main: ${OBJ}
#	${CANSI} ${CANSIFLAGS} -o math_example ${OBJ} ${LAPACKELIB} ${LIBS}
	${CANSI} ${CANSIFLAGS} -o math_example ${OBJ}  ${LAPACKELIB} ${LAPACKLIB} ${BLASLIB} ${TMGLIB} -lgfortran ${LIBS}
clean:
	rm -rf mySolveLinearLapack.o math_example 

