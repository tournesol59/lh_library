CC=/c/MinGW/bin/mingw32-g++
CANSI=/c/MinGW/bin/mingw32-gcc
CXXFLAGS= -std=c++11 -Wall -O2
CANSIFLAGS= -Wall


LOADER=/c/MinGW/bin/mingw32-g++
LOADOPTS='-O3'
HOME=/home/frede/lh_library
#library downloaded from netlib.org: one must link ALL of them
INCLAPACKE=${HOME}/lapack_inst/include
#library downloaded from netlib.org: one must link ALL of them
INCLAPACKE=${HOME}/lapack_inst/LAPACKE/include
#TMGLIB=/home/utilisateur/lapack_inst/libtmglib.a
LAPACKLIB=${HOME}/lapack_inst/liblapack.a
BLASLIB=${HOME}/lapack_inst/libblas.a
LAPACKELIB=${HOME}/lapack_inst/liblapacke.a

#library installed on linux system:
#INCLAPACKE=/home/frederic/lapack-3.8.0/LAPACKE/include
#TMGLIB=/usr/local/lib/libtmglib.a
#LAPACKLIB=/usr/local/lib/liblapack.a
#BLASLIB=/usr/local/lib/libblas.a

CINCLUDE=../include/
LIBS=-lm 

OBJ=math_example.o \
  mySolveLinearLapack.o

math_example: math_example.c
	${CANSI} ${CANSIFLAGS} -o math_example.o -c $< -I${INCLAPACKE}

mySolveLinearLapack: mySolveLinearLapack.c
	${CANSI} ${CANSIFLAGS} -o mySolveLinearLapack.o -c $< -I${INCLAPACKE}

main: ${OBJ}
	${CANSI} ${CANSIFLAGS} -o math_example  -I${INCLAPACKE} ${OBJ} ${BLASLIB} ${LAPACKELIB} ${LAPACKLIB} -lgfortran ${LIBS}
#	${CC} ${CXXFLAGS} -o math_example ${OBJ}  ${LAPACKELIB} ${LAPACKLIB} ${BLASLIB} ${TMGLIB} -lgfortran ${LIBS}
clean:
	rm -rf mySolveLinearLapack.o math_example 

