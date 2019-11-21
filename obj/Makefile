#
# General C++ Makefile by Anton Crombach, A.B.M.Crombach@bio.uu.nl
#

#CXX = g++4
CXX = g++


#CPPFLAGS = -w -march=x86-64 -mmmx -msse -msse2 -DDSFMT_MEXP=521 -DHAVE_SSE2 
CPPFLAGS = -w -O2 -march=x86-64 -mmmx -msse -msse2 -DDSFMT_MEXP=19937 -DHAVE_SSE2
#CPPFLAGS = -g -Wall -march=x86-64  -DDSFMT_MEXP=521 #-DHAVE_SSE2
#CPPFLAGS = -pg -Wall -march=x86-64
#CPPFLAGS = -ggdb3 -Wall -fprefetch-loop-arrays -DDSFMT_MEXP=19937
#CPPFLAGS = -Wall -pg -DDSFMT_MEXP=19937	#If you are doing this (GPROF), also put "-pg" after "$(CXX)" in the targets $(PROJECT) and $(LIBRARY).
# -ggdb3 is debug
# -g is line numbers
# -fstack-protector-all is stack protector for bt in gdb

INCDIR = -I../include -I$(MYPATH)/include -I/usr/include 
LIBDIR = -L$(MYPATH)/lib 

LIBS = -lpng -lgrace_np #-lefence 
# Source/object paths
vpath %.cc ../src ../test
vpath %.d .
vpath %.hh ../include
vpath %.o .

BINPATH = ../
LIBPATH = ../lib

# Project vars
PROJECT = Prok_Basic_Cycle
LIBRARY = mor
ALL =  Bead.o Gene.o TFBS.o Genome.o Prokaryote.o Population.o World.o dSFMT.o
OBJECTS = $(ALL)


# Targets
all: $(PROJECT)

$(PROJECT): $(OBJECTS) 
	$(CXX) $(LNFLAGS) $(LIBDIR) $^ -o $(BINPATH)/$@ $(LIBS)

$(LIBRARY): $(OBJECTS)
	$(CXX) -shared -Wl,-soname,$@.so \
	$(LIBDIR) $(LIBS) $^ -o $(LIBPATH)/$@.so


#$(OBJECTS): %.o: %.cc
%.o: %.cc
	$(CXX) -c $(CPPFLAGS) $(INCDIR) $< -o $@

%.d: %.cc
	@set -e; rm -f $@; \
	$(CXX) $(INCDIR) -MM $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),realclean)
-include $(OBJECTS:.o=.d)
endif
endif

.PHONY: clean realclean 
clean:
	@rm -f *.d.* *.o *.d 

realclean: clean 
	@cd ../src; rm -f *~ .*.swp *.bak; \
	cd ../include; rm -f *~ .*.swp *.bak
