all:
	g++ -Wall src/main.cpp -o ep.exe -O3 -lcuba -lm -lgsl -fopenmp