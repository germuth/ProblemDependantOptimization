#
# Makefile for program "cover" by Kari J. Nurmela.
#

CXX = g++
LIBS= -lm

#FLAGS= -ggdb -W # unoptimized, debuggable
FLAGS= -ffast-math  -W -Ofast -D_GLIBCXX_PARALLEL# optimized, not debuggable

#if windows
#REMOVE= del cover.exe
#DELETE = del
#if UNIX
REMOVE= rm -f cover
DELETE= rm -f

HEADERS= cover.h bincoef.h tables.h setoper.h solcheck.h exp.h arg.h pdo.h
OBJECTS= cover.o bincoef.o tables.o setoper.o solcheck.o exp.o arg.o pdo.o
CS=      cover.c bincoef.c tables.c setoper.c solcheck.c exp.c arg.c pdo.c
AUTOBACKUPS= *~

.c.o :	$(CS)
	$(CXX) $(FLAGS) -c $<

cover:	$(OBJECTS)
	$(REMOVE)
	$(CXX) $(FLAGS) -lgomp -o cover $(OBJECTS) \
	$(LIBS)

distrib:
	tar cvf cover.tar arg.c bincoef.c cover.c exp.c setoper.c solcheck.c tables.c arg.h bincoef.h cover.h exp.h setoper.h solcheck.h tables.h Makefile README
	compress cover.tar
	zip cover.zip arg.c bincoef.c cover.c exp.c setoper.c solcheck.c tables.c arg.h bincoef.h cover.h exp.h setoper.h solcheck.h tables.h Makefile README
	echo "NOTE: distribution contains files cover.tar.Z, cover.zip, and README."

clean:
	$(DELETE) $(OBJECTS) $(BCS) cover
	$(DELETE) $(AUTOBACKUPS)
	$(DELETE) core
	$(DELETE) cover.tar
	$(DELETE) cover.tar.Z
#	rm -f $(OBJECTS)
#	rm -f $(AUTOBACKUPS)
#	rm -f core
#	rm -f cover.tar
#	rm -f cover.tar.Z

cover.o: \
	cover.h \
	bincoef.h \
	tables.h \
	exp.h \
	solcheck.h

bincoef.o: \
	cover.h \
	bincoef.h

tables.o: \
	cover.h \
	tables.h \
	bincoef.h \
	setoper.h

setoper.o: \
	cover.h \
	setoper.h \
	bincoef.h

solcheck.o: \
	cover.h \
	setoper.h \
	bincoef.h \
	tables.h \
	solcheck.h

arg.o: \
	cover.h \
	arg.h \
	tables.h
