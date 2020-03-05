#pragma once
#include <iostream>

using namespace std;

class PoissonSolver
{
public:
    PoissonSolver();
    ~PoissonSolver();
    PoissonSolver(int nx, int ny, double ddx, double ddy);
    void CholSolve(double* rhsrhs);

private:
    int Nx;
    int Ny;
    double dx;
    double dy;
    double* a_banded;
   
	int internal_nodes;
	int ku;
	int k;
	double alpha; 
	double beta_x;
	double beta_y;
	int info;

};
