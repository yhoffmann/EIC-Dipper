CXX =g++
C =gcc

OPTIMOPTS =-O3
WOPTS =-Wall -Wextra -Wpedantic -pedantic-errors
CPPSTD =-std=c++20
LIBINCOPTS = -lm -lgsl -fopenmp
OPTS =$(OPTIMOPTS) $(WOPTS) $(EXTRAS)
CXX +=$(CPPSTD)


ifeq ($(QUIET),1)
	OPTS+= -D_QUIET
endif

ifeq ($(G2MU02), 1)
	OPTS+= -D_G2MU02
else
	ifeq ($(DILUTE), 1)
		OPTS+= -D_DILUTE
	endif
endif

ifeq ($(QBOTTOM), 1) # use bottom quark instead of charm
	OPTS+= -D_Q_B
endif

ifeq ($(QTOP), 1) # use top quark instead of charm
	OPTS+= -D_Q_T
endif

ifeq ($(CLUSTER),1)
	OPTS+= -lgslcblas
	OPTS+= -D_PC2
endif

ifeq ($(TEST), 1)
    OPTS+= -D_TEST
endif


local: dirs
	$(CXX) -o eic src/*.cpp obj/*.o $(OPTS) $(LIBINCOPTS)

dirs:
	bash -c 'mkdir -p interpolator-data'
	bash -c 'mkdir -p data/samples/{c05,c10,c20,b10,O16,Si28,S32,Ca40,Ni58,Cu62,Cu63,Xe129,W186,Au197,Pb207,Pb208}/{de,di}'
	bash -c 'mkdir -p data/samples/g2mu02/{c,b}/de'
	bash -c 'mkdir -p logs/{condor,slurm}'

.PHONY: external
external:
	mkdir -p obj
	$(C) -c -o obj/hcubature.o external/cubature/hcubature.c $(OPTS) 
	$(C) -c -o obj/pcubature.o external/cubature/pcubature.c $(OPTS) 
	$(CXX) -c -o obj/Nucleus.o external/Nucleus/src/Nucleus.cpp $(OPTS) 
	$(CXX) -c -o obj/HotspotNucleus.o external/Nucleus/src/HotspotNucleus.cpp $(OPTS) 

all: external local

.PHONY: debug
debug:
# mkdir -p debug
	@echo "\033[1;31mDEBUG DOES NOT WORK AND NEEDS TO BE UPDATED\033[0m"
# $(CXX) $(OPTS) -g -c external/cubature/hcubature.c -o debug/obj/hcubature.o -O0
# $(CXX) $(OPTS) -g -c external/cubature/pcubature.c -o debug/obj/pcubature.o -O0
# $(CXX) $(OPTS) -g -c external/Interpolation3D/src/*.cpp -o debug/obj/interpolation3d.o -fopenmp -O0
# $(CXX) $(OPTS) -g -c external/Interpolation3D/external/easy-progress-monitor/src/*.cpp -o debug/obj/easy-progress-monitor.o -O0
# $(CXX) $(OPTS) -g src/*.cpp debug/obj/*.o -o debug/eic -O0 -lcuba -lm -lgsl -fopenmp

cluster-single:
	make all QUIET=1 PCTWO=1

cluster:
	tar -czf interpolator-data.tar.gz interpolator-data/
	make cluster-single DILUTE=0
	mv eic eic-sde
	make cluster-single DILUTE=1
	mv eic eic-sdi

clean:
	rm obj/*
	rm eic

cleanRawData:
	rm -r data/samples/
