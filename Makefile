# defines
CXX=icc
SPEEDFLAGS=-O3
CXXFLAGS=-Wall -qopenmp $(SPEEDFLAGS)
LDFLAGS=-lrt -liomp5 $(SPEEDFLAGS)

TARGET = solution

all: $(TARGET)

# your own implementation
solution: main.o solution.o graph_tools.o
	$(CXX) $^ -o $@ $(LDFLAGS)

.cpp.o:
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean:
	rm -rf *.o $(TARGET)
