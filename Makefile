# Makefile for building the BuildTree shared library.

INCDIR = include
SRCDIR = src
OBJDIR = objects
LIBDIR = lib
CINTDIR = cint

#VPATH = src:include:objects:cint:lib
# Set the search path for prerequisites.
VPATH = $(INCDIR):$(SRCDIR):$(OBJDIR):$(LIBDIR):$(CINTDIR)

# C++ compiler
CXX = g++

#INCLUDE = -I./include/eic
INCLUDE = -I./include -I.
#INCLUDE = -I$(shell pwd)/include -I.

# Optimisation flag
OPT = -O2

# root-config --cflags includes ROOT include directory
# Add -Wlong-long is activiated by -pedantic.
# This will cause compilation failure due to ROOT routines, e.g.
#	  warning: use of C99 long long integer constant
# so suppress these warnings with -Wno-long-long.
CFLAGS = $(shell root-config --cflags) $(INCLUDE) $(OPT) -fPIC -pedantic \
	-Wall -Wextra -Wno-long-long -DUSE_NAMESPACE_ERHIC

# Additional library libEG required for TDatabasePDG/TParticlePDG.
LIBS=$(shell root-config --libs) -lEG

# Extract the name of all .h files in the include/ directory.
# Use grep -vi to exclude the linkdef file, which we need to
# ensure is the last argument passed to CINT.
HEADERS = $(notdir $(shell ls $(INCDIR)/*.h | grep -vi '*linkdef.h'))
#HEADERS = $(shell ls $(INCDIR)/*.h)
#HEADERS += $(shell ls include/eic/*.h)

# Source file extension
SRCEXT = cxx

# Object file extension
OBJEXT = o

# Lists of all source files to compile and the object files generate

# Extract the name of all source files in the src/ directory
#SOURCE = $(notdir $(shell ls src/*$(SRCEXT)))
SOURCE = $(shell ls src/*$(SRCEXT))

# Location of LinkDef file required for CINT dictionary generation.
LINKDEF = $(CINTDIR)/LinkDef.h

# Construct object names from source file names by substituting
# extension and directory path.
# The inner subst changes extension, the outer the directory.
#OBJECTS = $(subst $(SRCEXT),$(OBJEXT),$(SOURCE))
OBJECTS = $(subst $(SRCDIR)/,$(OBJDIR)/,$(subst $(SRCEXT),$(OBJEXT),$(SOURCE)))

# Implicit rule generating an object file in the objects directory
# from the corresponding source file in the source directory.
$(OBJDIR)/%.$(OBJEXT) : $(SRCDIR)/%.$(SRCEXT)
	mkdir -p objects
	$(CXX) -o $@ -c $(CFLAGS) $^

# Implicit rule for generating object file from source file.
%.$(OBJEXT) : %.$(SRCEXT)
	$(CXX) -o $@ -c $(CFLAGS) $^

# Builds the library
$(LIBDIR)/libBuildTree.so: $(CINTDIR)/BuildTreeDict.$(OBJEXT) $(OBJECTS)
	@mkdir -p lib
	$(CXX) -shared -o $(LIBDIR)/libBuildTree.so $(OBJECTS) \
		$(CINTDIR)/BuildTreeDict.$(OBJEXT) $(LIBS)
#	@ln -s $(PWD)/$(LIBDIR)/libBuildTree.so $(PWD)/$(LIBDIR)/libeicsmear.so
	ln -fs $(PWD)/$(LIBDIR)/libBuildTree.so $(PWD)/$(LIBDIR)/BuildTree.so

# Generates the CINT dictionary
$(CINTDIR)/BuildTreeDict.$(SRCEXT): $(HEADERS) $(LINKDEF)
	rootcint -f $(CINTDIR)/BuildTreeDict.$(SRCEXT) -c $(INCLUDE) $^

#erhic_tree: erhic_tree.o BuildTreeDict.o $(OBJECTS)
#	$(CXX) $(CFLAGS) -o $@ $^ $(LIBS)

#test.so: test.cxx testLinkDef.h
#	rootcint -f testDict.$(SRCEXT) -c $(INCLUDE) $^
#	$(CXX) $(CFLAGS) -c -o test.o test.cxx
#	$(CXX) $(CFLAGS) -c -o testDict.o testDict.cxx
#	$(CXX) -shared -o test.so test.o testDict.o $(LIBS)
#	$(CXX) -shared -o test.so test.o $(LIBS)
	
#test2: test2.o objects/PdgCode.o cint/BuildTreeD
#	$(CXX) -o test2 $^ $(LIBS)

headers:
	@echo $(HEADERS)

source:
	@echo $(SOURCE)

objs:
	@echo $(OBJECTS)

incdir:
	@echo $(INCLUDE)

.PHONY: clean

clean:
	rm -f $(CINTDIR)/*Dict* $(OBJDIR)/*.$(OBJEXT) $(LIBDIR)/*.so $(LIBDIR)/*.d

