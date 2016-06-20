CFLAGS=-O -w -Wno-traditional
LINTFLAGS=-lc -lm
SHELL=/bin/sh
CC=gcc
AR=ar

HFILES=spConfig.h spDefs.h spMatrix.h
CFILES=spAllocate.c spBuild.c spFactor.c spOutput.c spSolve.c spUtils.c
OFILES=spAllocate.o spBuild.o spFactor.o spOutput.o spSolve.o spUtils.o
LIBRARY=../lib/sparse.a
DESTINATION=sparse

SOURCE = $(HFILES) $(CFILES)

all: $(LIBRARY)

$(LIBRARY): $(OFILES)
	$(AR) r   $(LIBRARY) $?
	ranlib $(LIBRARY)

debug:
	make "CFLAGS=-g3 -O0 -Wno-traditional"

%.o : %.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $< -o $@

dep:
	$(CC) -MM $(CFILES) > depend

lint		:
	@lint $(LINTFLAGS) $(CFILES) $(TESTC) | grep -v "but never used"

clean		:
	rm -f $(OFILES) $(LIBRARY) $(TESTO) $(DESTINATION) spSMP.o core

touch		:
	touch -c  $(OFILES) $(LIBRARY)
	ranlib $(LIBRARY)

tags		: $(SOURCE) $(TESTC)
	ctags -t -w $(SOURCE) $(TESTC)
