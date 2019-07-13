#tip

CC=/c/mingw/bin/mingw32-g++
CANSI=/c/mingw/bin/mingw32-gcc
CXXFLAGS= -g -std=c++11 -Wall -O2
CANSIFLAGS= -Wall -O2
LOADER=/c/mingw/bin/mingw32-g++
LOADOPTS='-g -O3'

#library downloaded from netlib.org: one must link ALL of them
INCLAPACKE=/home/Utilisateur/lapack-3.8.0/LAPACKE/include
TMGLIB=/home/Utilisateur/lapack-3.8.0/libtmglib.a
LAPACKLIB=/home/Utilisateur/lapack-3.8.0/liblapack.a
BLASLIB=/home/Utilisateur/lapack-3.8.0/librefblas.a
LAPACKELIB=/home/Utilisateur/lapack-3.8.0/liblapacke.a

# std C++11 thus
CINCLUDE=../include/
#../../../Ipopt-3.12.6/include/coin/
LIBS=-lm

OBJ=../collocation/ident05_coll.o \
  main_example.o \
  ../MathFunctions/mySolveLinearLapack.o

main_example: main_example.cpp 
	${CC} ${CXXFLAGS} -o main_example.o -c $< -I${CINCLUDE}

ident05_coll: ../collocation/ident05_coll.cpp
	${CC} ${CXXFLAGS} -o ../collocation/ident05_coll.o -c $< -I${CINCLUDE}

mySolveLinearLapack: ../MathFunctions/mySolveLinearLapack.c
	${CANSI} ${CANSIFLAGS} -o ../MathFunctions/mySolveLinearLapack.o -c $< -I${INCLAPACKE}

#main no ref to std=c++11 !
main: ${OBJ}
	${LOADER} ${LOADOPTS} -o main_example -I${CINCLUDE} ${OBJ} ${LAPACKELIB}  ${LAPACKLIB} ${BLASLIB} ${TMGLIB} -lgfortran ${LIBS}

clean:
	rm -rf ../collocation/*.o ../MathFunctions/mySolveLinearLapack.o main_example 

