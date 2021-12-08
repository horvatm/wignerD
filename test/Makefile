CXX = g++
#CXX= icpc

DEPS = $(wildcard *.h)
INCLUDE = -I../include -I.    
LDFLAGS = 
CXXFLAGS = -O3 -Wall -fopenmp -std=c++17 

.PHONY: all clean

all: test_wignerD.e

%.e : %.cpp ${DEPS}
	${CXX} ${CXXFLAGS} ${INCLUDE} $< -o $@ ${LDFLAGS}

clean:
	rm -rf *.e *.dat 
