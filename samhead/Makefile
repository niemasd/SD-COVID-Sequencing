CXX=g++
CXXFLAGS?=-Wall -pedantic -O3 -std=c++11
OUTFILES=samhead

all: $(OUTFILES)

samhead: samhead.cpp
	$(CXX) $(CXXFLAGS) -o samhead samhead.cpp

clean:
	$(RM) $(OUTFILES) *.o
