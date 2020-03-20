#pragma once

#include <string>
#include <iostream>
#include <mpi.h>

using namespace std;

class LidDrivenCavity
{
public:
    LidDrivenCavity(); // constructor
    ~LidDrivenCavity(); // destructor
    

    // Methods:
    void SetDomainSize(double xlen, double ylen);
    void SetGridSize(int nx, int ny);
    void SetPartitions(int px, int py);
    void SetTimeStep(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);

    void Initialise(int world_rank, int world_size);
    void Communication();
    void Integrate(LidDrivenCavity &LDC);


    // Friendship
    friend void arguments(int argc, char** argv, LidDrivenCavity &LDC);
    friend class PoissonSolver;

private:

	// Vorticity and streamfunction
    double* v  = nullptr;
    double* vt = nullptr;
    double* s  = nullptr;
    
    // Problem inputs
    double dt;
    double T;
    int    Nx,  Ny;
    double Lx,  Ly;

    // Problem conditions
    double Re;
    const double U = 1.0;
    
	// Algorithm variables
	double deltax;
	double deltay;
	double dlxy;
	double pdltx;
	double pdlty;
	
	// Variables for location inside arrays	
	int firstpointx;
	int firstpointy; 
	int lastpointx;
	int lastpointy; 
	// UP, DOWN, LEFT and RIGHT first and last points
	int fpUP, 	 lpUP;
	int fpDOWN,  lpDOWN;
	int fpLEFT,  lpLEFT;
	int fpRIGHT, lpRIGHT;
	
	// Cartesian communicator
  	MPI_Comm   cart_comm;
  	
	// Parallelisation variables
	int    myrank;	    
    int	   nx,  ny;  
    int	   Npx, Npy;
    int	   Np;
    int    Px,  Py;     
	int	   ndims = 2;

	// Communication between processes variables	
	//     dims[ndims], period[ndims], coord[ndims], nbrs[ndims*2];
	int    dims[2], 	period[2], 	   coord[2], 	 nbrs[2*2];
	int    UP = 0, DOWN = 1, LEFT = 2, RIGHT = 3;
	int    tag 	   = 1;
	int    reorder = 0;
	int    tasks   = 2; 
	int	   nreqs   = 0;
	
	// vectors for send and recieve
	double *svsendu = nullptr;
	double *svrecvu = nullptr;
	double *svsendd = nullptr;
	double *svrecvd = nullptr;
	double *svsendl = nullptr;
	double *svrecvl = nullptr;
	double *svsendr = nullptr;
	double *svrecvr = nullptr;
	
	// For testing the code
	int rankdisplay;
	

  	
};

