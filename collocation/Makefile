#tip

CC=g++
CANSI=gcc
CXXFLAGS= -g -std=c++11 -Wall -O2 # -Wstringop-truncation
CANSIFLAGS= -g -Wall -O2
LOADER=g++
LOADOPTS='-g'

HOME=/home/frede/lh_library/
#library downloaded from netlib.org: one must link ALL of them
INCLAPACKE=${HOME}/lapack_inst/include
#TMGLIB=/home/utilisateur/lapack_inst/libtmglib.a
LAPACKLIB=$(HOME)/lapack_inst/liblapack.a
# change here libblas to librefblas
BLASLIB=${HOME}/lapack_inst/librefblas.a
LAPACKELIB=${HOME}/lapack_inst/liblapacke.a

# std C++11 thus
CINCLUDE=../include/ \

#../../../Ipopt-3.12.6/include/coin/
LIBS=-lm

OBJ=../collocation/ident05_collinput.o \
    ../collocation/ident05_coll.o \
    ../collocation/ident05_colloutput.o \
    ../main_cpp/ident05_fftdata.o \
    ../MathFunctions/mySolveLinearLapack.o \
    main_example.o 

main_example: main_example.cpp 
	${CC} ${CXXFLAGS} -o main_example.o -c $< -I${CINCLUDE} -I${INCLAPACKE}

ident05_collinput: ../collocation/ident05_collinput.cpp
	${CC} ${CXXFLAGS} -o ident05_collinput.o -c $< -I${CINCLUDE} -I${INCLAPACKE}
	
ident05_coll: ../collocation/ident05_coll.cpp
	${CC} ${CXXFLAGS} -o ident05_coll.o -c $< -I${CINCLUDE} -I${INCLAPACKE}

ident05_colloutput: ../collocation/ident05_colloutput.cpp
	${CC} ${CXXFLAGS} -o ident05_colloutput.o -c $< -I${CINCLUDE} -I${INCLAPACKE}

ident05_fftdata: ../main_cpp/ident05_fftdata.cpp
	${CC} ${CXXFLAGS} -o ident05_fftdata.o -c $< -I${CINCLUDE} -I${INCLAPACKE}

mySolveLinearLapack: ../MathFunctions/mySolveLinearLapack.c
	${CANSI} ${CANSIFLAGS} -o ../MathFunctions/mySolveLinearLapack.o -c $< -I${INCLAPACKE}

#main no ref to std=c++11 !
main: ${OBJ}
	${LOADER} ${LOADOPTS} -o main_example -I${CINCLUDE} -I${INCLAPACKE} ${OBJ} ${LAPACKELIB}  ${LAPACKLIB} ${BLASLIB} -lgfortran ${LIBS}

clean:
	rm -rf ../collocation/*.o ../main_cpp/*.o ../MathFunctions/mySolveLinearLapack.o main_example 

#g++ -std=c++11 -o main_example -I..\include main_example.o  -I..\..\lapack_inst\include ..\main_cpp\ident05_fftdata.o ..\collocation\ident05_coll.o ..\collocation\ident05_collinput.o ..\MathFunctions\mySolveLinearLapack.o ..\..\lapack_inst\liblapacke.a ..\..\lapack_inst\liblapack.a ..\..\lapack_inst\librefblas.a -lgfortran -lm
#-lgfortran raises an error, what should i do
