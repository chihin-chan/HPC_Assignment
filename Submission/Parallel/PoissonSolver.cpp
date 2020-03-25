#include "PoissonSolver.h"
#include <iostream>
#include <cblas.h>
#include <iostream>
#include <algorithm>
using namespace std;
#define UP    0
#define DOWN  1
#define LEFT  2
#define RIGHT 3
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

PoissonSolver::~PoissonSolver(){ 
  // Clean Matrices
   delete[] a_banded;
   delete[] rhs;
}

// Class Member function that stores Cholesky Factorisation
void PoissonSolver::CholFact(LidDrivenCavity &src){
	
	// Initialising subdomain sizes
	x_off = 2;
	y_off = 2;
	
	if (src.nghbrs[UP] == -2){
		y_off += 1;
	}
	if (src.nghbrs[DOWN] == -2){
		y_off += 1;
	}
	if (src.nghbrs[RIGHT] == -2){
		x_off += 1;
	}
	if (src.nghbrs[LEFT] == -2){
		x_off += 1;
	}

	// Initialising parameters for storing banded matrix a_banded
	internal_nodes = (src.loc_nx-x_off)*(src.loc_ny-y_off);
	ku = src.loc_nx-x_off;			  // No. of super diagonals
	k = (ku+1)*internal_nodes;	  // Size of banded matrix (ku+kl+1)*N
	alpha = 2.0*(1.0/src.dx/src.dx + 1.0/src.dy/src.dy); // Coefficients of i,j
	beta_x = -1.0/src.dx/src.dx;		// Coefficients of i+/-1, j
	beta_y = -1.0/src.dy/src.dy;		// Coefficients of i/, j+/-1
	a_banded = new double[k];
	rhs = new double[internal_nodes];
	src.rhs = new double[internal_nodes];
	
	// Initialising pointers
	fill_n(a_banded, k, 0.0);
	fill_n(rhs, internal_nodes, 0.0);

	// Storing Elements of a_banded in banded format
	int count = 1;
	a_banded[ku] = alpha;
	for(int i=2*(ku+1)-1; i<k; i+=(ku+1)){
		a_banded[i] = alpha;
		// Storing beta_x
		if(count % (src.loc_nx-x_off) != 0){
			a_banded[i-1] = beta_x;
		}
		else{
			a_banded[i-1] = 0.0;
		}
		// Storing beta_y
		if(i > (ku+1)*(ku)){
			a_banded[i-ku] = beta_y;
		}
        count++;
	}

	// Factorising Caching
	F77NAME(dpbtrf) ('u', internal_nodes, ku, a_banded, ku+1, info);
}

void PoissonSolver::CholSolve(double* rhs){
	// Forward and backward substitution
	F77NAME(dpbtrs) ('U', internal_nodes, ku, 1, a_banded, ku+1, rhs, internal_nodes, info);
    //cout << "Cholesky Solver was called" << endl;
    //cout << "Info: " << info << endl;
}
