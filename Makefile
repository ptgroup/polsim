CXX ?= clang++
CXXFLAGS ?= -O2 -g
OUTNAME ?= polsim

INCLUDEDIR = include
SRCDIR = src
BUILDDIR = build

CXXFLAGS := $(CXXFLAGS) -std=c++14 -Wall -Wextra -pedantic -I$(INCLUDEDIR)

OBJS = controller.o main.o pdp.o simulation.o
OBJ_PATHS = $(OBJS:%=$(BUILDDIR)/%)

.PHONY: all clean docs run

all: $(BUILDDIR)/$(OUTNAME)

clean:
	rm -rf $(BUILDDIR)

docs:
	doxygen

run: $(BUILDDIR)/$(OUTNAME)
	$(BUILDDIR)/$(OUTNAME)

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(BUILDDIR)/$(OUTNAME): $(OBJ_PATHS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ_PATHS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<
