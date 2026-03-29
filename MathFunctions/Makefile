CC=/c/MinGW/bin/mingw32-g++
CANSI=/c/MinGW/bin/mingw32-gcc
CXXFLAGS= -std=c++11 -Wall -O2
CANSIFLAGS= -Wall


LOADER=/c/MinGW/bin/mingw32-g++
LOADOPTS='-O3'
#library downloaded from netlib.org: one must link ALL of them
INCLAPACKE=/home/utilisateur/lapack-3.8.0/LAPACKE/include
#TMGLIB=/home/utilisateur/lapack_inst/libtmglib.a
LAPACKLIB=/home/utilisateur/lapack_inst/liblapack.a
BLASLIB=/home/utilisateur/lapack_inst/libblas.a
LAPACKELIB=/home/utilisateur/lapack_inst/liblapacke.a

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
	${CANSI} ${CANSIFLAGS} -o math_example.o -c $<

mySolveLinearLapack: mySolveLinearLapack.c
	${CANSI} ${CANSIFLAGS} -o mySolveLinearLapack.o -c $< -I${INCLAPACKE}

main: ${OBJ}
	${CANSI} ${CANSIFLAGS} -o math_example ${OBJ} ${BLASLIB} ${LAPACKELIB} ${LAPACKLIB} -lgfortran ${LIBS}
#	${CC} ${CXXFLAGS} -o math_example ${OBJ}  ${LAPACKELIB} ${LAPACKLIB} ${BLASLIB} ${TMGLIB} -lgfortran ${LIBS}
clean:
	rm -rf mySolveLinearLapack.o math_example 

