CC=g++
CXX=g++
CXXFLAGS=-std=c++11 -g3 -O0 -Wall -Wno-deprecated -UDEBUG -I ../src
#LDLIBS= ../lib/c++hdl.a ../lib/sparse.a
LD=ld
LDFLAGS=-std=c++11 -g3 -O0 -Wall
#LDFLAGS=-pg

EXAMPLES=rc

all: $(EXAMPLES)

rc: rc.o ../lib/c++hdl.a ../lib/sparse.a
rcckt: rcckt.o ../lib/c++hdl.a ../lib/sparse.a
diockt: diockt.o ../lib/c++hdl.a ../lib/sparse.a
bipckt: bipckt.o ../lib/c++hdl.a ../lib/sparse.a
example2: example2.o ../lib/c++hdl.a ../lib/sparse.a
example4: example4.o ../lib/c++hdl.a ../lib/sparse.a
example5: example5.o ../lib/c++hdl.a ../lib/sparse.a
pll: pll.o ../src/pll.h ../lib/c++hdl.a ../lib/sparse.a
example7: example7.o ../lib/c++hdl.a ../lib/sparse.a
example8: example8.o ../lib/c++hdl.a ../lib/sparse.a
diffpair: diffpair.o ../lib/c++hdl.a ../lib/sparse.a
ua709: ua709.o ../lib/c++hdl.a ../lib/sparse.a
ua741: ua741.o ../lib/c++hdl.a ../lib/sparse.a
ua741x: ua741x.o ../lib/c++hdl.a ../lib/sparse.a
ua727: ua727.o ../lib/c++hdl.a ../lib/sparse.a
rtlinv: rtlinv.o ../lib/c++hdl.a ../lib/sparse.a
ttlinv: ttlinv.o ../lib/c++hdl.a ../lib/sparse.a
rca3040: rca3040.o ../lib/c++hdl.a ../lib/sparse.a
vco: ../src/pll.h vco.o ../lib/c++hdl.a ../lib/sparse.a

%.o : %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

build:
	make "CXXFLAGS=-O0 -Wall -UDEBUG -I ../src"

debug:
	make "CXXFLAGS=-std=c++0x -g3 -O0 -Wall -Wno-deprecated -DEBUG -I ../src"

prof:
	make "CXXFLAGS= -pg -I ../src -O0 -Wall -UDEBUG"

dep:
	$(CXX) -I ../src -MM *.cpp > depend

include depend

clean:
	rm -f $(EXAMPLES) *.o *.exe *.out *.tran *.dc 

