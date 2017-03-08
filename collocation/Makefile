CC=/usr/bin/g++
CXXFLAGS= -std=c++11 -Wall -O2
# std C++11 thus
CINCLUDE=../include/
#../../../Ipopt-3.12.6/include/coin/
LIBS=-lm

OBJ=../collocation/ident05_coll.o \
  main_example.o

main_example: main_example.cpp 
	${CC} ${CXXFLAGS} -o main_example.o -c $< -I${CINCLUDE}

ident05_coll: ../collocation/ident05_coll.cpp
	${CC} ${CXXFLAGS} -o ../collocation/ident05_coll.o -c $< -I${CINCLUDE}

main: ${OBJ}
	${CC} ${CXXFLAGS} -o main_example -I${CINCLUDE} ${OBJ} ${LIBS}

clean:
	rm -rf ../collocation/*.o main_example 

