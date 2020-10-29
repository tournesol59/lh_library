#tip
CC=/c/mingw/bin/mingw32-g++
CANSI=/c/mingw/bin/mingw32-gcc
CXXFLAGS= -g -std=c++11 -Wall -O2
CANSIFLAGS= -Wall -O2
LOADER=/c/mingw/bin/mingw32-g++
LOADOPTS= -g 

#library downloaded from netlib.org: one must link ALL of them
INCLAPACKE=/home/Utilisateur/lapack-3.8.0/LAPACKE/include
#TMGLIB=/home/Utilisateur/lapack_inst/libtmglib.a
LAPACKLIB=/home/Utilisateur/lapack_inst/liblapack.a
BLASLIB=/home/Utilisateur/lapack_inst/librefblas.a
LAPACKELIB=/home/Utilisateur/lapack_inst/liblapacke.a

# std C++11 thus
CINCLUDE=../include/ 

#../../../Ipopt-3.12.6/include/coin/
LIBS=-lm

OBJ=../recursive_LSQ/generate.o \
    ../recursive_LSQ/impl_relsq.o \
    ../main_cpp/ident05_fftdata.o \
  test_relsq.o

test_relsq: test_relsq.cpp 
	${CC} ${CXXFLAGS} -o test_relsq.o -c $< -I${CINCLUDE}
	
test_reest: test_reest.cpp 
	${CC} ${CXXFLAGS} -o test_reest.o -c $< -I${CINCLUDE}

ident05_fftdata: ../main_cpp/ident05_fftdata.cpp
	${CC} ${CXXFLAGS} -o ../main_cpp/ident05_fftdata.o -c $< -I${CINCLUDE}

generate: ../recursive_LSQ/generate.cpp
	${CC} ${CXXFLAGS} -o ../recursive_LSQ/generate.o -c $< -I${CINCLUDE}

impl_relsq: ../recursive_LSQ/impl_relsq.cpp
	${CC} ${CXXFLAGS} -o ../recursive_LSQ/impl_relsq.o -c $< -I${CINCLUDE}

#main no ref to std=c++11 !
main: ${OBJ}
#	${LOADER} ${LOADOPTS} -o test_reest -I${CINCLUDE} ${OBJ} ${LIBS}
	${LOADER} ${LOADOPTS} -o test_relsq -I${CINCLUDE} ${OBJ} ${LAPACKELIB}  ${LAPACKLIB} ${BLASLIB} -lgfortran ${LIBS}

clean:
	rm -rf ../recursive_LSQ/*.o test_relsq.exe
