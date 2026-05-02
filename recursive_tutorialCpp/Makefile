#tip
CC=g++
CANSI=gcc
CXXFLAGS= -g -std=c++11 -Wall -O2
CANSIFLAGS= -Wall -O2
LOADER=g++
LOADOPTS= -g 

#CC=/c/mingw/bin/mingw32-g++
#CANSI=/c/mingw/bin/mingw32-gcc
#CXXFLAGS= -g -std=c++11 -Wall -O2
#CANSIFLAGS= -Wall -O2
#LOADER=/c/mingw/bin/mingw32-g++
#LOADOPTS= -g 

#library downloaded from netlib.org: one must link ALL of them
#INCLAPACKE=../../lapack-3.8.0/LAPACKE/include
#TMGLIB=/home/frede/lapack_inst/libtmglib.a
#LAPACKLIB=../../lapack_inst/liblapack.a
#BLASLIB=../../lapack_inst/librefblas.a
#LAPACKELIB=../../lapack_inst/liblapacke.a

# std C++11 thus
CINCLUDE=./ 

#../../../Ipopt-3.12.6/include/coin/
LIBS=-lm

OBJ=../recursive_tutorialCpp/test_Computation.o \
    ../recursive_tutorialCpp/impl_Computation.o \

test_Computation: ./test_Computation.cpp 
	${CC} ${CXXFLAGS} -o test_Computation.o -c $< -I${CINCLUDE}

impl_Computation: ./impl_Computation.cpp
	${CC} ${CXXFLAGS} -o impl_Computation.o -c $< -I${CINCLUDE}

#main no ref to std=c++11 !
main: ${OBJ}
	${LOADER} ${LOADOPTS} -o test_Computation -I${CINCLUDE} ${OBJ}    
#	${LAPACKELIB}  ${LAPACKLIB} ${BLASLIB} -lgfortran ${LIBS}

clean:
	rm -rf *.o test_Computation.exe
	
