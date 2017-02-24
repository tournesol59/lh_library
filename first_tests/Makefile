CC=/usr/bin/gcc
CXXFLAGS=-Wall -O2
CINCLUDE=../include/
#../../../Ipopt-3.12.6/include/coin/
LIBS=-lm

OBJ=../MathFunctions/myComputeScalarProduct.o \
  ../galerkin/ident05_gakm.o \
  ident05_test.o

ident05_test: ident05_test.cpp
	${CC} ${CXXFLAGS} -o ident05_test.o -c $< -I${CINCLUDE}

myComputeScalarProduct: ../MathFunctions/myComputeScalarProduct.cpp
	${CC} ${CXXFLAGS} -o ../MathFunctions/myComputeScalarProduct.o -c $<

ident05_gakm: ../galerkin/ident05_gakm.cpp
	${CC} ${CXXFLAGS} -o ../galerkin/ident05_gakm.o -c $< -I${CINCLUDE}

main: ${OBJ}
	${CC} ${CXXFLAGS} -o ident05_test -I${CINCLUDE} ${OBJ} ${LIBS}

clean:
	rm -rf ../MathFunctions/*.o ../galerkin/*.o ident05_test

