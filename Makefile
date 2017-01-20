# defines
CXX=g++
SPEEDFLAGS=-O3
CXXFLAGS=-Wall $(SPEEDFLAGS)
LDFLAGS=-lrt $(SPEEDFLAGS)

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
