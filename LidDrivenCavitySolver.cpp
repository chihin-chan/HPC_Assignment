#include <iostream>
#include <boost/program_options.hpp>
using namespace std;
namespace po = boost::program_options;
#include "LidDrivenCavity.h"
#include <mpi.h>

int main(int argc, char **argv)
{	
	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("Lx", po::value<double>(), "Length of the domain in x-direction")
		("Ly", po::value<double>(), "Length of the domain in y-direction")
		("Nx", po::value<int>(), "Number of grid points in x-direction")
		("Ny", po::value<int>(), "Number of grid points in y-direction")
		("Px", po::value<int>(), "Number of partitions in the x-direction (parallel)")
		("Py", po::value<int>(), "Numer of partitions in the y-direction (parallel)")
		("dt", po::value<double>(), "Time step size")
		("T", po::value<int>(), "Final Time")
		("Re", po::value<double>(), "Reynolds Number");
	
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);   
	
	if (vm.count("help")) {
    	cout << desc << "\n";
    	return 1;
	}
	

    	// Create a new instance of the LidDrivenCavity class
    	LidDrivenCavity* solver = new LidDrivenCavity();

	solver->SetDomainSize(vm["Lx"].as<double>(), vm["Ly"].as<double>());
	solver->SetGridSize(vm["Nx"].as<int>(), vm["Ny"].as<int>());
	solver->SetTimeStep(vm["dt"].as<double>());
	solver->SetFinalTime(vm["T"].as<int>());
	
	/*
    // Configure the solver here...
    // ...
    solver->Initialise();

    // Run the solver
    solver->Integrate();
*/
	return 0;
}
