CXX ?= clang++
CXXFLAGS ?= -O2 -g
OUTNAME ?= polsim

CXXFLAGS := $(CXXFLAGS) -std=c++14 -Wall -Wextra -Werror

OBJS = main.o pdp.o simulation.o

.PHONY: all build clean docs

all: build

build: $(OUTNAME)

clean:
	rm -f *.o $(OUTNAME)

docs:
	doxygen

$(OUTNAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<
