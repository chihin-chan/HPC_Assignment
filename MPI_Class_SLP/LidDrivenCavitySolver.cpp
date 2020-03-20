#include <iostream>
#include <boost/program_options.hpp>
using namespace std;
namespace po = boost::program_options;
#include "LidDrivenCavity.h"
#include <mpi.h>

int main(int argc, char **argv)
{	
	int s, r;
	// Parsing Values from Command Line
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
		("T", po::value<double>(), "Final Time")
		("Re", po::value<double>(), "Reynolds Number");
	
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);   
	
	if (vm.count("help")) {
    	cout << desc << "\n";
    	return 1;
	}
	// Initialising MPI
    MPI_Init(&argc, &argv);
	
	// Retrieving MPI Size
	MPI_Comm_size(MPI_COMM_WORLD, &s);

	// Retreving Rank
    MPI_Comm_rank(MPI_COMM_WORLD, &r);

	// Checking if NxNy matches number of processes
	try{
		if(vm["Px"].as<int>() * vm["Py"].as<int>() != s){
			throw std::out_of_range("");
		}
	}
	catch(std::out_of_range const &e){
		if (r==0){
			cout << "Invalid Number of Processes requested" << endl;
		}
	}
    // Create a new instance of the LidDrivenCavity class
    LidDrivenCavity* solver = new LidDrivenCavity();
	solver->SetDomainSize(vm["Lx"].as<double>(), vm["Ly"].as<double>());
	solver->SetGridSize(vm["Nx"].as<int>(), vm["Ny"].as<int>());
	solver->SetPartitionSize(vm["Px"].as<int>(), vm["Py"].as<int>());
	solver->SetTimeStep(vm["dt"].as<double>());
	solver->SetFinalTime(vm["T"].as<double>());
	solver->SetReynoldsNumber(vm["Re"].as<double>());
	solver->GetRank(r);
	solver->GetSize(s);
    // Configure the solver here...
	solver->Initialise();
	
    // Run the solver
    solver->Integrate();

	// Export Results
	//solver->ExportSol();
	
	MPI_Finalize();
	return 0;
}
