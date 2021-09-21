#tip

CC=/c/mingw/bin/mingw32-g++
CANSI=/c/mingw/bin/mingw32-gcc
CXXFLAGS= -g -std=c++11 -Wall -O2 -Wstringop-truncation
CANSIFLAGS= -g -Wall -O2
LOADER=/c/mingw/bin/mingw32-g++
LOADOPTS='-g'

#library downloaded from netlib.org: one must link ALL of them
INCLAPACKE=/home/utilisateur/lapack-3.8.0/LAPACKE/include
#TMGLIB=/home/utilisateur/lapack_inst/libtmglib.a
LAPACKLIB=/home/utilisateur/lapack_inst/liblapack.a
BLASLIB=/home/utilisateur/lapack_inst/libblas.a
LAPACKELIB=/home/utilisateur/lapack_inst/liblapacke.a

# std C++11 thus
CINCLUDE=../include/
#../../../Ipopt-3.12.6/include/coin/
LIBS=-lm

OBJ=../collocation/ident05_coll.o \
    ../main_cpp/ident05_fftdata.o \
    ../MathFunctions/mySolveLinearLapack.o \
    main_example.o 

main_example: main_example.cpp 
	${CC} ${CXXFLAGS} -o main_example.o -c $< -I${CINCLUDE}

ident05_coll: ../collocation/ident05_coll.cpp
	${CC} ${CXXFLAGS} -o ../collocation/ident05_coll.o -c $< -I${CINCLUDE}

ident05_fftdata: ../main_cpp/ident05_fftdata.cpp
	${CC} ${CXXFLAGS} -o ../main_cpp/ident05_fftdata.o -c $< -I${CINCLUDE}

mySolveLinearLapack: ../MathFunctions/mySolveLinearLapack.c
	${CANSI} ${CANSIFLAGS} -o ../MathFunctions/mySolveLinearLapack.o -c $< -I${INCLAPACKE}

#main no ref to std=c++11 !
main: ${OBJ}
	${LOADER} ${LOADOPTS} -o main_example -I${CINCLUDE} ${OBJ} ${LAPACKELIB}  ${LAPACKLIB} ${BLASLIB} -lgfortran ${LIBS}

clean:
	rm -rf ../collocation/*.o ../main_cpp/*.o ../MathFunctions/mySolveLinearLapack.o main_example 

