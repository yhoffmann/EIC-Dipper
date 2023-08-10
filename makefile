all:
	g++ -Wall src/*.cpp Interpolation3D/src/*.cpp -o eic -O3 -lcuba -lm -lgsl -fopenmp