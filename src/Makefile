AR=ar
CXX=g++
CXXFLAGS=-std=c++11 -g3 -O0 -Wall -Wno-deprecated -Winline -UDEBUG
ARCHIVES=../lib/c++hdl.a

SRCS=util.cpp global.cpp node.cpp quantity.cpp equation.cpp model.cpp \
parameter.cpp port.cpp matrix.cpp sim.cpp sta.cpp util.cpp indent_facet.cpp
OBJS=util.o global.o node.o quantity.o equation.o model.o parameter.o port.o \
sim.o sta.o matrix.o util.o indent_facet.o
LIBRARY=../lib/c++hdl.a

all: $(LIBRARY)

$(LIBRARY): $(OBJS)
	$(AR) r ../lib/c++hdl.a $?
	ranlib ../lib/c++hdl.a

%.o : %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

dep:
	$(CXX) -MM $(SRCS) > depend

include depend

build:
	make "CXXFLAGS=-O2 -Wall -UDEBUG"

debug:
	make "CXXFLAGS=-std=c++0x -g3 -O0 -Wall -Wno-deprecated -DDEBUG"

prof:
	make "CXXFLAGS= -pg -O0 -Wall -UDEBUG"

opt:
	make "CXXFLAGS=-O"

ci: $(HEADERS) $(SRCS)
	ci -l $(HEADERS) $(SRCS) Makefile

clean:
	rm -f $(OBJS) core 
