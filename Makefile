# defines
CXX=icc

#-xCORE-AVX2 -mtune=core-avx2
SPEEDFLAGS=-O3 -Ofast
#PROFGENFLAGS=-prof-gen=threadsafe -prof-dir=./
#PROFUSEFLAGS=-prof-use -prof-dir=./
#-debug inline-debug-info
CXXFLAGS= -Wall -qopenmp $(SPEEDFLAGS) $(PROFGENFLAGS) $(PROFUSEFLAGS)
LDFLAGS=-lrt -liomp5 $(CULIBS) $(SPEEDFLAGS)

TARGET = solution reference_bfs gen_RMAT gen_random gen_valid_info validation

all: $(TARGET)

# your own implementation
solution: main.o solution.o graph_tools.o
	$(CXX) $^ -o $@ $(LDFLAGS)

# reference implementation with bfs
reference_bfs: main.o reference_bfs.o graph_tools.o
	$(CXX) $^ -o $@ $(LDFLAGS)

# RMAT generator
gen_RMAT: gen_RMAT.o graph_tools.o
	$(CXX) $^ -o $@ $(LDFLAGS)

# Erdos-Renyi (random) graph generator
gen_random: gen_random.o graph_tools.o
	$(CXX) $^ -o $@ $(LDFLAGS)

# gen_valid_info
gen_valid_info: gen_valid_info.o reference_bfs.o graph_tools.o
	$(CXX) $^ -o $@ $(LDFLAGS)

# validation
validation: validation.o reference.o graph_tools.o
	$(CXX) $^ -o $@ $(LDFLAGS)

.cpp.o:
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -rf *.o $(TARGET)

