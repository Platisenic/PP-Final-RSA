CXX = g++
LDL = -lgmpxx -lgmp

all: sol

sol: sol.cpp
	$(CXX) -o $@ $< $(LDL)

clean:
	rm -rf sol
