#pragma once
#include <iostream>
#include "LidDrivenCavity.h"

using namespace std;

class PoissonSolver
{
public:
    PoissonSolver();
    ~PoissonSolver();
	
	void CholFact(LidDrivenCavity &src);
	void CholSolve(double* rhs);
private:

	int internal_nodes;
	int ku;
	int x_off;
	int y_off;
	int k;
	int info;
	double alpha;
	double beta_x;
	double beta_y;

	double* a_banded;
	double* rhs;
	
};
