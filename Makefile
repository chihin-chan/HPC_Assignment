.PHONY: clean
default: myprog

LidDrivenCavitySolver.o: LidDrivenCavitySolver.cpp
	mpicxx -o LidDrivenCavitySolver.o -c LidDrivenCavitySolver.cpp -llapack -lblas -lboost_program_options

LidDrivenCavity.o: LidDrivenCavity.cpp LidDrivenCavity.h
	mpicxx -o LidDrivenCavity.o -c LidDrivenCavity.cpp -llapack -lblas -lboost_program_options

myprog: LidDrivenCavity.o LidDrivenCavitySolver.o
	mpicxx -o myprog LidDrivenCavity.o LidDrivenCavitySolver.o -llapack -lblas -lboost_program_options

clean:
	rm -f *.o myprog
