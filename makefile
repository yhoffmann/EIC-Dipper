all:
	g++ -Wall src/*.cpp -o eic -O3 -lcuba -lm -lgsl -fopenmp