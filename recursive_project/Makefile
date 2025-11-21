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
CINCLUDE=../include/ 

#../../../Ipopt-3.12.6/include/coin/
LIBS=-lm

OBJ=../recursive_LSQ/generate.o \
    ../recursive_project/impl_iogenerate.o \
    ../recursive_project/impl_basicrelsq.o \
    ../main_cpp/ident05_fftdata.o \
  test_basiclsq.o

test_basiclsq: ./test_basiclsq.cpp 
	${CC} ${CXXFLAGS} -o test_basiclsq.o -c $< -I${CINCLUDE}

ident05_fftdata: ../main_cpp/ident05_fftdata.cpp
	${CC} ${CXXFLAGS} -o ../main_cpp/ident05_fftdata.o -c $< -I${CINCLUDE}

generate: ../recursive_LSQ/generate.cpp
	${CC} ${CXXFLAGS} -o ../recursive_LSQ/generate.o -c $< -I${CINCLUDE}

impl_basicrelsq: ./impl_basicrelsq.cpp
	${CC} ${CXXFLAGS} -o impl_basicrelsq.o -c $< -I${CINCLUDE}

imp_iogenerate: ./impl_iogenerate.cpp
	${CC} ${CXXFLAGS} -o impl_iogenerate.o -c $< -I${CINCLUDE}

#main no ref to std=c++11 !
main: ${OBJ}
	${LOADER} ${LOADOPTS} -o test_basiclsq -I${CINCLUDE} ${OBJ}    
#	${LAPACKELIB}  ${LAPACKLIB} ${BLASLIB} -lgfortran ${LIBS}

clean:
	rm -rf ../recursive_LSQ/*.o test_basiclsq.exe
	
