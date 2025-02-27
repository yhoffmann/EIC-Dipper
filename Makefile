.PHONY: local external all debug

FLAGS =-Wall -O3

ifeq ($(QUIET),1)
	FLAGS+= -D_QUIET
endif

ifeq ($(G2MU02), 1)
	FLAGS+= -D_G2MU02
else
	ifeq ($(DILUTE), 1)
		FLAGS+= -D_DILUTE
	endif
endif

ifeq ($(QBOTTOM), 1) # use bottom quark instead of charm
	FLAGS+= -D_Q_B
endif

ifeq ($(QTOP), 1) # use bottom quark instead of charm
	FLAGS+= -D_Q_T
endif

ifeq ($(PCTWO),1)
	FLAGS+= -lgslcblas
	FLAGS+= -D_PC2
	# FLAGS+= -I external/cuba
	# FLAGS+= -L external/cuba
endif

ifeq ($(TEST), 1)
    FLAGS+= -D_TEST
endif


local:
	bash -c 'mkdir -p interpolator-data'
	bash -c 'mkdir -p data/samples/{c05,c10,c20,b10}/{de,di}'
	bash -c 'mkdir -p data/samples/g2mu02/{c,b}/de'
	g++ $(FLAGS) src/*.cpp obj/*.o -o eic -lm -lgsl -fopenmp

external:
	mkdir -p obj
	g++ $(FLAGS) -c -o obj/hcubature.o  external/cubature/hcubature.c
	g++ $(FLAGS) -c -o obj/pcubature.o external/cubature/pcubature.c
	g++ $(FLAGS) -c -o obj/Nucleus.o external/Nucleus/src/Nucleus.cpp
	g++ $(FLAGS) -c -o obj/HotspotNucleus.o external/Nucleus/src/HotspotNucleus.cpp

all:
	make external
	make local

debug:
# mkdir -p debug
	@echo "\033[1;31mDEBUG DOES NOT WORK AND NEEDS TO BE UPDATED\033[0m"
# g++ $(FLAGS) -g -c external/cubature/hcubature.c -o debug/obj/hcubature.o -O0
# g++ $(FLAGS) -g -c external/cubature/pcubature.c -o debug/obj/pcubature.o -O0
# g++ $(FLAGS) -g -c external/Interpolation3D/src/*.cpp -o debug/obj/interpolation3d.o -fopenmp -O0
# g++ $(FLAGS) -g -c external/Interpolation3D/external/easy-progress-monitor/src/*.cpp -o debug/obj/easy-progress-monitor.o -O0
# g++ $(FLAGS) -g src/*.cpp debug/obj/*.o -o debug/eic -O0 -lcuba -lm -lgsl -fopenmp

pc2:
	make all QUIET=1 PCTWO=1

clean:
	rm obj/*
	rm eic

cleanRawData:
	rm -r data/samples/