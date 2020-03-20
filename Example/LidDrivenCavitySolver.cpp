#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <iterator>
using namespace std;

#include "LidDrivenCavity.h"
#include "commandarguments.h"
#include "PoissonSolver.h"
#include <mpi.h>

int main(int argc, char **argv)
{   
    // Initialise MPI
	MPI_Init(&argc, &argv);
	
	
	// Get the rank and size in the original communicator
	int world_rank, world_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	
    // Create an instance of the LidDrivenCavity class
    LidDrivenCavity* solver = new LidDrivenCavity();
	

	// Function for the different options and its arguments
	// introduced by the user in the command line 
	arguments(argc,  argv, *solver, world_rank);



    solver->Initialise(world_rank, world_size);
    
    // For each time increment
    //for (double dlt = solver->dt; dlt<solver->T; dlt += solver->dt ){
        // Run the solver
   	solver->Integrate(*solver);

    //}
	
	
	// Finalize MPI
	MPI_Finalize();
	
	return 0;
}
