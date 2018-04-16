# automatic settings for root

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

# Linux with egcs, gcc 2.9x, gcc 3.x (>= RedHat 5.2)
CXX           := $(shell root-config --cxx)
CXXFLAGS      = -O -Wall -fPIC -I$(ROOTSYS)/include
LD            := $(shell root-config --ld)  -Wl,-rpath,$(ROOTSYS)/lib
LDFLAGS       = -O
SOFLAGS       = -shared
#

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

# sources
SRCS		= fmake_tree.C
OBJS		= fmake_tree.o

# program
PROGRAM		= fmake_tree

all:	$(PROGRAM)
	@echo "all done"

$(PROGRAM):	$(OBJS)
		@echo "Linking $(PROGRAM) ..."
		$(LD) $(LDFLAGS) $(OBJS) $(LIBS) $(GLIBS) -o $(PROGRAM)
		@echo "done"

clean:;		@rm -f $(OBJS) core
###
# make_tree.o:	m2root.h

#MyClass.o:	MyClass.h
#mydict.C:	MyClass.h
#	@echo "Generating Dictionary ..."
#	@rootcint mydict.C -c MyClass.h
