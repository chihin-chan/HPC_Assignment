.PHONY: clean
default: myprog

LidDrivenCavitySolver.o: LidDrivenCavitySolver.cpp
	mpicxx -o LidDrivenCavitySolver.o -c LidDrivenCavitySolver.cpp -lblas -llapack -lboost_program_options

PoissonSolver.o: PoissonSolver.cpp PoissonSolver.h
	mpicxx -o PoissonSolver.o -c PoissonSolver.cpp -lblas -llapack

LidDrivenCavity.o: LidDrivenCavity.cpp LidDrivenCavity.h PoissonSolver.h
	mpicxx -o LidDrivenCavity.o -c LidDrivenCavity.cpp -lblas -llapack  -lboost_program_options

myprog: LidDrivenCavity.o LidDrivenCavitySolver.o
	mpicxx -o myprog -PoissonSolver.o LidDrivenCavity.o LidDrivenCavitySolver.o -lblas -llapack -lboost_program_options

clean:
	rm -f *.o myprog
