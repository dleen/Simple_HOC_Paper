spike.out : main.cpp LIF_spike.cpp LIF_spike.h LIF_constants.h \
	    DG_gen_spike_matrix.cpp \
	    LIF_gen_spike_matrix.cpp \
	    EIF_gen_spike_matrix.cpp \
	    Macke_figures.cpp
#	gcc -o spike.out main.cpp LIF_spike.cpp LIF_gen_spike_matrix.cpp Macke_figures.cpp -lgsl -lgslcblas -lm -O3 -ffast-math -funroll-loops -fomit-frame-pointer \
	-msse -march=corei7
#	g++ -o spike.out main.cpp LIF_spike.cpp LIF_gen_spike_matrix.cpp Macke_figures.cpp -lgsl -lgslcblas -lm -Wall -pedantic 
#	clang++ -o spike.out main.cpp LIF_spike.cpp LIF_gen_spike_matrix.cpp Macke_figures.cpp -lgsl -lgslcblas -lm -Wall -pedantic
#	/Users/dleen/my_gcc/bin/g++-4.6.2 -o spike_med.out main.cpp LIF_spike.cpp LIF_spike.h LIF_constants.h \
	DG_gen_spike_matrix.cpp \
	LIF_gen_spike_matrix.cpp \
	EIF_gen_spike_matrix.cpp \
	Macke_figures.cpp \
	-lm -lgsl -lgslcblas -Ofast -ffast-math -msse -march=corei7 \
	-Wall -m64
	clang++ -o spike_med.out main.cpp LIF_spike.cpp \
	DG_gen_spike_matrix.cpp \
	LIF_gen_spike_matrix.cpp \
	EIF_gen_spike_matrix.cpp \
	Macke_figures.cpp \
	-lm -lgsl -lgslcblas \
	-Wall -m64 -O3

clean:
	rm figure*.dat spike.out
