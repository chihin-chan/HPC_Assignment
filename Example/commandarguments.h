#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <iterator>
using namespace std;

#include "LidDrivenCavity.h"
void arguments(int argc, char** argv, LidDrivenCavity &LDC, int my_rank){

    try {
        // Describe the options we want to have for our program
        // These are stored in a variable of type "options_description":
        po::options_description desc("Allowed options");

        // We specify the options by using the add_options() function
        // Options may take no arguments, or a single argument
        // The middle parameter specifies the type of parameter
        // The last parameter is a description displayed on the help message.
        desc.add_options()
            ("help", "produce help message")
            ("Lx", po::value<double>(), "length of the domain in the x-direction")
            ("Ly", po::value<double>(), "length of the domain in the y-direction")
            ("Nx", po::value<int>(), "number of grid points in x-direction")
            ("Ny", po::value<int>(), "number of grid points in y-direction")
            ("Px", po::value<int>(), "number of partitions in the x-direction (parallel)")
            ("Py", po::value<int>(), "number of partitions in the y-direction (parallel)")
            ("dt", po::value<double>(), "time step size")
            ("T", po::value<double>(), "final time")
            ("Re", po::value<double>(), "Reynolds number")
        ;

        // These statements tell Boost program_options to parse the command-line
        // arguments given to main().
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
        

        // If the user gives the --help argument, print the help.
        if (vm.count("help")) {
            cout << desc << "\n";
        }

        // For each user entered -- flag, print out the value and set it to the variable.
        // Otherwise notify that it was not given.
	cout << " ------------------------------------------------- Processor " << my_rank << endl;
	        
        // ----------------------------------------------------------------------
        // 															  DOMAIN SIZE
        if (vm.count("Lx")) {
            cout << "Length of the domain in the x-direction was set to " 
                 << vm["Lx"].as<double>() << ".\n";
        } else {
            cout << "Length of the domain in the x-direction was not set.\n";
        }
        if (vm.count("Ly")) {
            cout << "Length of the domain in the y-direction was set to " 
                 << vm["Ly"].as<double>() << ".\n";
        } else {
            cout << "Length of the domain in the y-direction was not set.\n";
        }
        if (vm.count("Lx") && (vm.count("Ly"))) {
        	LDC.SetDomainSize(vm["Lx"].as<double>(),vm["Ly"].as<double>());
        }
        
        // ----------------------------------------------------------------------
        // 															    GRID SIZE
        if (vm.count("Nx")) {
            cout << "Number of grid points in x-direction was set to " 
                 << vm["Nx"].as<int>() << ".\n";
        } else {
            cout << "Number of grid points in x-direction was not set.\n";
        }
        if (vm.count("Ny")) {
            cout << "Number of grid points in y-direction was set to " 
                 << vm["Ny"].as<int>() << ".\n";
        } else {
            cout << "Number of grid points in y-direction was not set.\n";
        }
        if (vm.count("Nx") && (vm.count("Ny"))) {
        	LDC.SetGridSize(vm["Nx"].as<int>(),vm["Ny"].as<int>());
        }
                
        // ----------------------------------------------------------------------
        // 													 NUMBER OF PARTITIONS
        if (vm.count("Px")) {
            cout << "Number of partitions in the x-direction was set to " 
                 << vm["Px"].as<int>() << ".\n";
        } else {
            cout << "Number of partitions in the x-direction was not set.\n";
        }
        if (vm.count("Py")) {
            cout << "Number of partitions in the y-direction was set to " 
                 << vm["Py"].as<int>() << ".\n";
        } else {
            cout << "Number of partitions in the y-direction was not set.\n";
        }
        if (vm.count("Px") && (vm.count("Py"))) {
        	LDC.SetPartitions(vm["Px"].as<int>(),vm["Py"].as<int>());
        }
        
        // ---------------------------------------------------------------------
        // 													 	   TIME SETTINGS
        if (vm.count("dt")) {
            cout << "Time step size was set to " 
                 << vm["dt"].as<double>() << ".\n";
            LDC.SetTimeStep(vm["dt"].as<double>());
        } else {
            cout << "Time step size was not set.\n";
        }
        if (vm.count("T")) {
            cout << "Final time was set to " 
                 << vm["T"].as<double>() << ".\n";
            LDC.SetFinalTime(vm["T"].as<double>());
        } else {
            cout << "Final time was not set.\n";
        }
        
        // ---------------------------------------------------------------------
        // 													 	 REYNOLDS NUMBER
        if (vm.count("Re")) {
            cout << "Reynolds number was set to " 
                 << vm["Re"].as<double>() << ".\n";
            LDC.SetReynoldsNumber(vm["Re"].as<double>());
        } else {
            cout << "Reynolds number was not set.\n";
        }       
    }
    // Catch any exceptions thrown (which derived from std::exception)
    // e.g. std::runtime_error, std::logic_error, etc.
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
    }
    // The "..." in a catch block catches absolutely any exception thrown, no
    // matter what type it has, that hasn't been caught already by the above.
    catch(...) {
        cerr << "Exception of unknown type!\n";
    }
}
