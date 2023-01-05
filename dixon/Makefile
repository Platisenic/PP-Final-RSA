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
run_mpi: dixon_mpi
	 mpirun -mca btl_tcp_if_include 192.168.2.0/24 -np 36  --hostfile hosts dixon_mpi  622955649784207776129637054061200033
