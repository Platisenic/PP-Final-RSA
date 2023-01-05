CXX = g++
CFLAGS = -std=c++17 -O3 -g -Wall
LIBS = -lgmpxx -lgmp

.PHONY: fall clean

all: trialFactor

trialFactor: trialFactor.cpp CycleTimer.h
	$(CXX) $(CFLAGS) -o $@ $< $(LIBS)

clean:
	rm -rf trialFactor
