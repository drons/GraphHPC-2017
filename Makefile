# defines
CXX=icc

#-xCORE-AVX2 -mtune=core-avx2
SPEEDFLAGS=-O3 -Ofast  -funroll-loops -qopt-prefetch=5
#PROFGENFLAGS=-prof-gen=threadsafe -prof-dir=/home/sudorgin/intel
#PROFUSEFLAGS=-prof-use -prof-dir=/home/sudorgin/intel
#-debug inline-debug-info
CXXFLAGS= -Wall -qopenmp $(SPEEDFLAGS) $(PROFGENFLAGS) $(PROFUSEFLAGS)
LDFLAGS=-lrt -liomp5 $(CULIBS) $(SPEEDFLAGS)

TARGET = validation gen_RMAT reference solution

all: $(TARGET)

# your own implementation
solution: main.o solution.o graph_tools.o
	$(CXX) $^ -o $@ $(LDFLAGS)

# reference implementation
reference: main.o reference.o graph_tools.o
	$(CXX) $^ -o $@ $(LDFLAGS)

# RMAT generator
gen_RMAT: gen_RMAT.o graph_tools.o
	$(CXX) $^ -o $@ $(LDFLAGS)

# validation
validation: validation.o reference.o graph_tools.o
	$(CXX) $^ -o $@ $(LDFLAGS)

.cpp.o:
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -rf *.o $(TARGET)

