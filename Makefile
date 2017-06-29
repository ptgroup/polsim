CXX ?= clang++
CXXFLAGS ?= -O2 -g
OUTNAME ?= polsim

CXXFLAGS := $(CXXFLAGS) -std=c++14 -Wall -Werror

OBJS = main.o simulation.o

.PHONY: all build clean docs

all: build

build: $(OUTNAME)

clean:
	rm *.o $(OUTNAME)

docs:
	doxygen

$(OUTNAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<
