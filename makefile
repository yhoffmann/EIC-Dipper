.PHONY: nolibs libs all debug

FLAGS =

ifeq ($(QUIET),1)
	FLAGS+= -D_QUIET
endif

ifeq ($(DILUTE),1) # use dilute approx
	FLAGS+= -D_DILUTE
endif

ifeq ($(QBOTTOM), 1) # use bottom quark instead of charm
	FLAGS+= -D_Q_B
endif

ifeq ($(QTOP), 1) # use bottom quark instead of charm
	FLAGS+= -D_Q_T
endif

ifeq ($(G2MU02), 1)
	FLAGS+= -D_G2MU02
endif

ifeq ($(PCTWO),1)
	FLAGS+= -lgslcblas
	FLAGS+= -I external/cuba
	FLAGS+= -L external/cuba
endif

ifneq ($(PCTWO), 1)
	FLAGS+= -D_INTERP_LOG
endif


nolibs:
	mkdir -p data/samples/b10/de data/samples/b10/di
	mkdir -p data/samples/c05/de data/samples/c05/di
	mkdir -p data/samples/c10/de data/samples/c10/di
	mkdir -p data/samples/c20/de data/samples/c20/di
	mkdir -p data/samples/g2mu02/c
	g++ $(FLAGS) -Wall src/*.cpp obj/*.o -o eic -O3 -lcuba -lm -lgsl -fopenmp

libs:
	mkdir -p obj
	g++ $(FLAGS) -c -Wall external/cubature/hcubature.c -o obj/hcubature.o -O3
	g++ $(FLAGS) -c -Wall external/cubature/pcubature.c -o obj/pcubature.o -O3
	g++ $(FLAGS) -c -Wall external/Interpolation3D/src/Interpolator3D.cpp -o obj/interpolation3d.o -O3 -fopenmp
	g++ $(FLAGS) -c -Wall external/Interpolation3D/external/easy-progress-monitor/src/*.cpp -o obj/easy-progress-monitor.o -O3
	g++ $(FLAGS) -c -Wall external/Nucleus/src/Nucleus.cpp -o obj/Nucleus.o -O3
	g++ $(FLAGS) -c -Wall external/Nucleus/src/HotspotNucleus.cpp -o obj/HotspotNucleus.o -O3
	g++ $(FLAGS) -c -Wall external/thread-pool/src/ThreadPool.cpp -o obj/ThreadPool.o -O3

all:
	make libs
	make nolibs

debug:
	echo "DEBUG DOES NOT WORK AND NEEDS TO BE UPDATED"
	g++ $(FLAGS) -g -c -Wall external/cubature/hcubature.c -o debug/obj/hcubature.o -O0
	g++ $(FLAGS) -g -c -Wall external/cubature/pcubature.c -o debug/obj/pcubature.o -O0
	g++ $(FLAGS) -g -c -Wall external/Interpolation3D/src/*.cpp -o debug/obj/interpolation3d.o -fopenmp -O0
	g++ $(FLAGS) -g -c -Wall external/Interpolation3D/external/easy-progress-monitor/src/*.cpp -o debug/obj/easy-progress-monitor.o -O0
	g++ $(FLAGS) -g -Wall src/*.cpp debug/obj/*.o -o debug/eic -O0 -lcuba -lm -lgsl -fopenmp

pc2:
	make all QUIET=1 PCTWO=1
