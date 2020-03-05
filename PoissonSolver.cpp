#include "PoissonSolver.h"
#include <iostream>
#include <cblas.h>
#define F77NAME(x) x##_

extern "C" {
	void F77NAME(dpbtrf) (const char& UPLO, const int& n, const int& kd,
		       	      const double* ab, const int& LDAB, int& info);
	void F77NAME(dpbtrs) (const char& UPLO, const int& n, const int& KD,
			      const int& nrhs, const double* AB, const int& ldab,
			      const double* rhsrhs, const int& ldb, int& info);
}

using namespace std;

// Default Constructor
PoissonSolver::PoissonSolver(){}

// User-Define Constructor
PoissonSolver::PoissonSolver(int nx, int ny, double ddx, double ddy){
    
    // Storing Private Variables
    Nx = nx;
    Ny = ny;
    dx = ddx;
    dy = ddy;
	internal_nodes = (Nx-2)*(Ny-2);         // Ax = b solved only for internal nodes
	ku = Nx-2;			                    // No. of super diagonals
	alpha = 2.0*(1.0/dx/dx + 1.0/dy/dy);    // Coefficients of i,j
	beta_x = -1/dx/dx;		                // Coefficients of i+/-1, j
	beta_y = -1/dy/dy;		                // Coefficients of i/, j+/-1

	k = (ku+1)*internal_nodes;	            // Size of banded matrix (ku+kl+1)*N
	a_banded = new double[k];               // a_banded holds matrix A in banded format

	// Move the initialising to poisson solver contructor
	// Generating Banded Matrix in column format for Possion Solve
	// Initialising a_banded to all zeros
	for(unsigned int i=0; i<k; i++){
		a_banded[i] = 0.0;
	}
	// Filling Bottom and Bottom +1 rows with alpha and beta_x
	a_banded[ku] = alpha;
	int count = 1;
	for(unsigned int i=ku+(ku+1); i<k; i+=(Nx-1)){
		a_banded[i] = alpha;
		if (count%(Nx-2) != 0){
			a_banded[i-1] = beta_x;
		}
		else{
			a_banded[i-1] = 0;
		}
		count++;
	}

	// Filling Top Row of Banded Matrix with beta_y
	for(unsigned int i=ku*(ku+1); i<k; i+=(Nx-1)){
		a_banded[i] = beta_y;
	}
/*	
	// Printing A_banded for checking
	for(unsigned int i = 0; i < (ku+1) ; i++){
		for(unsigned int j = 0; j < internal_nodes; j++){
			cout << a_banded[i+j*(ku+1)] << " ";
		}
		cout << endl;
	}
*/
	// Caching Cholesky factorisation
	F77NAME(dpbtrf) ('u', internal_nodes, ku, a_banded, ku+1, info);		
    cout << endl << endl << "Cholesky Factorisation was called" << endl <<endl;
/*
	// Printing A_banded for checking
	for(unsigned int i = 0; i < (ku+1) ; i++){
		for(unsigned int j = 0; j < internal_nodes; j++){
			cout << a_banded[i+j*(ku+1)] << " ";
		}
		cout << endl;
	}
*/

}
PoissonSolver::~PoissonSolver(){ 
    // Clean Matrices
//    delete[] a_banded;
  //  delete[] b;
}

void PoissonSolver::CholSolve(double* rhsrhs){
    cout << "here" << endl;  
    F77NAME(dpbtrs) ('U', internal_nodes, ku, 1, a_banded, ku+1, rhsrhs, internal_nodes, info);
    cout << "Cholesky Solver was called" << endl;
    cout << "Info: " << info << endl;
}
