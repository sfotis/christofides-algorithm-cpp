cc=g++
cflags=-c -std=c++11

all: driver

driver: tsp.cpp tsp.h
	g++ tsp.cpp -o tsp
	
clean:
	rm -f a.out
	rm -f *.o
	rm -f driver
