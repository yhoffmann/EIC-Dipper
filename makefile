.PHONY: nolibs libs all debug

FLAGS =

ifeq ($(QUIET),1)
	FLAGS += -D_QUIET
endif

ifeq ($(PCTWO),1)
	FLAGS += -lgslcblas
	FLAGS += -I external/cuba
	FLAGS += -L external/cuba
endif


nolibs:
	g++ $(FLAGS) -Wall src/*.cpp obj/*.o -o eic -O3 -lcuba -lm -lgsl -fopenmp

libs:
	g++ $(FLAGS) -c -Wall external/cubature/hcubature.c -o obj/hcubature.o -O3
	g++ $(FLAGS) -c -Wall external/cubature/pcubature.c -o obj/pcubature.o -O3
	g++ $(FLAGS) -c -Wall external/Interpolation3D/src/*.cpp -o obj/interpolation3d.o -O3 -fopenmp
	g++ $(FLAGS) -c -Wall external/Interpolation3D/external/easy-progress-monitor/src/*.cpp -o obj/easy-progress-monitor.o -O3
	g++ $(FLAGS) -c -Wall external/Nucleus/src/Nucleus.cpp -o obj/Nucleus.o -O3
	g++ $(FLAGS) -c -Wall external/Nucleus/src/HotspotNucleus.cpp -o obj/HotspotNucleus.o -O3

all:
	make libs
	make nolibs

debug:
	g++ $(FLAGS) -g -c -Wall external/cubature/hcubature.c -o debug/obj/hcubature.o -O0
	g++ $(FLAGS) -g -c -Wall external/cubature/pcubature.c -o debug/obj/pcubature.o -O0
	g++ $(FLAGS) -g -c -Wall external/Interpolation3D/src/*.cpp -o debug/obj/interpolation3d.o -fopenmp -O0
	g++ $(FLAGS) -g -c -Wall external/Interpolation3D/external/easy-progress-monitor/src/*.cpp -o debug/obj/easy-progress-monitor.o -O0
	g++ $(FLAGS) -g -Wall src/*.cpp debug/obj/*.o -o debug/eic -O0 -lcuba -lm -lgsl -fopenmp

pc2:
	make all QUIET=1 PCTWO=1
