.PHONY: clean
default: myprog


PoissonSolver.o: PoissonSolver.cpp PoissonSolver.h
	mpicxx -Wall -g -o PoissonSolver.o -c PoissonSolver.cpp -lblas -llapack

LidDrivenCavitySolver.o: LidDrivenCavitySolver.cpp
	mpicxx -Wall -g -o LidDrivenCavitySolver.o -c LidDrivenCavitySolver.cpp -lblas -llapack -lboost_program_options

LidDrivenCavity.o: LidDrivenCavity.cpp LidDrivenCavity.h
	mpicxx -Wall -g -o LidDrivenCavity.o -c LidDrivenCavity.cpp -lblas -llapack  -lboost_program_options

myprog: LidDrivenCavity.o LidDrivenCavitySolver.o PoissonSolver.o
	mpicxx -Wall -g -o myprog PoissonSolver.o LidDrivenCavity.o LidDrivenCavitySolver.o -lblas -llapack -lboost_program_options

clean:
	rm -f *.o myprog
