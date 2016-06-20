
ARCHIVES=../lib/c++hdl.a ../lib/sparse.a
CC=g++
CXX=g++
CXXFLAGS=-std=c++0x -O -Wall -Wno-deprecated -UDEBUG -I ../src
LDLIBS= ../lib/c++hdl.a ../lib/sparse.a
SPICE=spice
SPFLAGS= -g

TESTS=d1 q1
SPICERUNS=d1.dat

all: $(TESTS) $(SPICERUNS)

d1: d1.o ../lib/c++hdl.a ../lib/sparse.a
q1: q1.o ../lib/c++hdl.a ../lib/sparse.a

%.o : %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

%.dat : %.spi
	$(SPICE) $(SPFLAGS) $<

debug:
	make "CXXFLAGS=-g3 -I ../src -O0 -Wall -DDEBUG" "LDLIBS= ../lib/c++hdl.a ../lib/sparse.a"

dep:
	$(CXX) -I ../src -MM *.cpp > depend

include depend

clean:
	rm -f $(TESTS) *.o *.out *.tran *.dc
