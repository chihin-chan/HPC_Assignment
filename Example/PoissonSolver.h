#pragma once

#include "LidDrivenCavity.h"

#include <string>
#include <iostream>
using namespace std;

class PoissonSolver
{
public:
    PoissonSolver(); // constructor
    ~PoissonSolver(); // destructor

    void Initializing(LidDrivenCavity& LDC);    
    void MatrixConstr(LidDrivenCavity& LDC);
    void FillingRHS(LidDrivenCavity& LDC);
    void CartesianToScalapack(LidDrivenCavity& LDC);
    void ScalapackToCartesian(LidDrivenCavity& LDC);
    void Solver(LidDrivenCavity& LDC);

private:

	double* A  	  = nullptr;	// A matrix symmetric banded
	double* bcart	= nullptr; 	// b vector in cartesian format per processor
	double* bcartnpe = nullptr;	// b vector in cartesian format for all processors
	double* bscanpe  = nullptr;	// b vector in ScaLapack format for all processors
    double* b 	  = nullptr;	// In: RHS vector, Out: Solution
    double* workf = nullptr; 	// Workspace factorization
    double* works = nullptr; 	// Workspace solve  
    int*	ipiv  = nullptr; 	// Pivoting array
    double* AF	  = nullptr;	// Auxiliary Fillin space
    
    // Indexes
	int innernx;	// Internal points in x-dir 
	int innerny;	// Internal points in y-dir 
	int innern;		// Total number of points in each processor
	
    // Packed storage variables
    int BWU;	// Upper bandwidth
    int BWL; 	// Lower bandwidth
    int NA;		// Total size problem
    int NB; 	// Blocking size
   	int nrhs;	// Number of rhs to solve
   	int JA;		// Start offset in matrix
   	int IB; 	// Start offset in vector 
   	int LA;		// Leading dimension A * NB
   	int lda;	// Leading dimension A
   	int LAF;  	// Size of Auxiliary Fillin space
   	int lworkf; // Size of workspace work (factorization)
   	int lworks; // Size of workspace work (solve)
   	int N;		// Number of points in interior processors
   	
   	// CBLACS variables for parallel linear algebra
   	int mype, npe, ctx, nrow, ncol, myrow, mycol;
   	   	   	    
    // Descriptors for matrix A and rhs vector b
    int* desca 	 = nullptr;
    int* descb	 = nullptr;
       	
   	// A matrix constants
   	double diag;
   	double offdiagi;
   	double offdiagj;
   	
   	// Status value
   	int info;
   	
   	// Gathering vectors
	int*    innernxV = nullptr;	
	int*    innernyV = nullptr;
	
	// Cartesian to Scalapack vectors and variables
	int* new_inNV     = nullptr;		
	int* new_innernxV = nullptr;
	int* new_innernyV = nullptr;
	int localny, localnx, cproc;
	int prevny, prevnx, prevn;
	int k1, k2;	
	int q, z;	
	
	

};
