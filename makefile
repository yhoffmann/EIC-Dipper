all:
	g++ -Wall src/main.cpp src/coherent.cpp src/incoherent.cpp src/integration_routines.cpp src/models_and_functions.cpp src/utilities.cpp -o eic -O3 -lcuba -lm -lgsl -fopenmp