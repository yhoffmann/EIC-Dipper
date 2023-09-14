nolibs:
	g++ -Wall src/*.cpp obj/*.o -o eic -O3 -lcuba -lm -lgsl -fopenmp

libs:
	g++ -c -Wall cubature/hcubature.c -o obj/hcubature.o -O3
	g++ -c -Wall cubature/pcubature.c -o obj/pcubature.o -O3
	g++ -c -Wall Interpolation3D/src/*.cpp -o obj/interpolation3d.o -O3 -fopenmp
	g++ -c -Wall Interpolation3D/easy-progress-monitor/src/*.cpp -o obj/easy-progress-monitor.o -O3

all:
	make libs
	make nolibs