CXX = g++
MPIXX= mpicxx
LDL = -lgmpxx -lgmp -O3
CV = c++17
OMP = -fopenmp
all: dixon dixon_omp dixon_mpi 

dixon: dixon_serial.cpp
	$(CXX) -o $@ $< $(LDL)
dixon_omp: dixon_omp.cpp
	$(CXX) -std=$(CV) $(OMP) -o $@ $< $(LDL)
dixon_mpi: dixon_mpi.cpp
	$(MPIXX) -std=c++11 -o $@ $< $(LDL)
clean:
	rm -rf dixon dixon_omp dixon_mpi
