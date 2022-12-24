CXX = g++
LDL = -lgmpxx -lgmp
CV = c++17
OMP = -fopenmp -O3
all: dixon dixon_omp

dixon: dixon_serial.cpp
	$(CXX) -o $@ $< $(LDL)
dixon_omp: dixon_omp.cpp
	$(CXX) -std=$(CV) $(OMP) -o $@ $< $(LDL)

clean:
	rm -rf dixon dixon_omp
