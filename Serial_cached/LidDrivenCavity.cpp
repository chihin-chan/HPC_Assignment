#include "LidDrivenCavity.h"
#include <iostream>
#include <cstring>
#include <cblas.h>
#include <math.h>
#include <fstream>
#include <iomanip>
#define F77NAME(x) x##_
extern "C" {
	void F77NAME(dpbtrf) (const char& UPLO, const int& n, const int& kd,
		       	      const double* ab, const int& LDAB, int& info);
	void F77NAME(dpbtrs) (const char& UPLO, const int& n, const int& KD,
			      const int& nrhs, const double* AB, const int& ldab,
			      const double* b, const int& ldb, int& info);
}

using namespace std;

LidDrivenCavity::LidDrivenCavity()
{
}

LidDrivenCavity::~LidDrivenCavity()
{
}

void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
{ 
	Lx = xlen;
	Ly = ylen;
	cout << "Domain Size, Lx: " << Lx <<"	Ly: " << Ly << endl;
}

void LidDrivenCavity::SetGridSize(int nx, int ny)
{
	Nx = nx;
	Ny = ny;
	cout << "Grid Size, Nx: " << Nx <<"	Ny: " << Ny << endl;
}

void LidDrivenCavity::SetTimeStep(double deltat)
{
	dt = deltat;
	cout << "Time Step: " << dt << endl;
}

void LidDrivenCavity::SetFinalTime(double finalt)
{
	T = finalt;
	cout << "Final Time: " << T << endl;
}

void LidDrivenCavity::SetReynoldsNumber(double re)
{
	Re = re;
	cout << "Reynolds Number: " << Re << endl;
}

void LidDrivenCavity::MatPrint(double *x, int n)
{
	for(int i = 0 ; i<n; i++){
		for(int j = 0; j<n; j++){
			cout << x[i+j*n] <<"	" ;
		}
		cout << endl;
	}
}	

void LidDrivenCavity::Initialise()
{	
	v = new double[Nx*Ny];
	s = new double[Nx*Ny];
	v_new  = new double[Nx*Ny];
	s_new = new double[Nx*Ny];
	// Initialising Vorticity, w and Streamfunction, s to zero
	memset(v, 0.0, Nx*Ny*sizeof(double));
	memset(s, 0.0, Nx*Ny*sizeof(double));	
	memset(v_new, 0.0, Nx*Ny*sizeof(double));	
	memset(s_new, 0.0, Nx*Ny*sizeof(double));
	
	cout << endl << endl <<"Initialising solution..." << endl << endl;
	
}

void LidDrivenCavity::Integrate()
{
	double dx = double(Lx) / double((Nx-1.0));
	double dy = double(Ly) / double((Ny-1.0));
	double U = 1.0;
	double t_elapse = 0.0;

	// Initialising properties of banded matrix A, to solve Ax = b
	const int internal_nodes = (Nx-2)*(Ny-2); // Ax = b solved only for internal nodes
	const int ku = Nx-2;			  // No. of super diagonals
	const int k = (ku+1)*internal_nodes;	  // Size of banded matrix (ku+kl+1)*N
	double alpha = 2.0*(1.0/dx/dx + 1.0/dy/dy); // Coefficients of i,j
	double beta_x = -1/dx/dx;		// Coefficients of i+/-1, j
	double beta_y = -1/dy/dy;		// Coefficients of i/, j+/-1
	double norm;				// Storing 2-norm of solution difference
	int info;
	int count;
	
	// a_banded holds matrix A in banded format
	double* a_banded = new double[k];
	// rhs stores vector b for Ax = b
	double* rhs = new double[internal_nodes];
	for(int i = 0; i< internal_nodes; i++){
		rhs[i] = 0.0;
	}

	// Try and Catch if dt is too large
	double cond = Re*Lx/(Nx-1)*Ly/(Ny-1) / 4;
	try{	
		if(dt >= cond){throw std::out_of_range("");}
	}
	catch(std::out_of_range const &e){
		cout << "dt: " << dt << 
		", is too large and should be lesser than: " << cond << endl;
	      return;	      
	}

	// Generating Banded Matrix in column format for Possion Solve
	// Initialising a_banded to all zeros
	for(unsigned int i=0; i<k; i++){
		a_banded[i] = 0.0;
	}
	// Filling Bottom and Bottom +1 rows with alpha and beta_x
	a_banded[ku] = alpha;
	count = 1;
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
	
	// Starting time loop
	while (t_elapse < T){
		
		// Enforcing Left Boundary Conditions for Voticity, w
		for(unsigned int i = Nx; i < Nx*Ny; i+=Nx){
			v[i] = (s[i] - s[i+1])*2.0/dx/dx;
		}
		// Enforcing Bottom Boundary Conditions for Vorticity, w
		for(unsigned int i = 0; i < Nx; i++){
			v[i] = (s[i]-s[i+Nx])*2.0/dy/dy;
		}

		// Enforcing Right Boundary Conditions for Vorticity, w
		for(unsigned int i = Nx-1; i < Nx*Ny; i+=Nx){
			v[i] = (s[i] - s[i-1])*2.0/dx/dx;
		}

		// Enforcing Top Boundary Conditions for Vorticity, w
		for(unsigned int i = Nx*(Ny-1); i < Nx*Ny; i++){
			v[i] = (s[i]-s[i-Nx])*2.0/dy/dy - 2.0*U/dy;
		}
		// Calculation of Interior Vorticity at time t
		for(unsigned int j = 1; j<Ny-1; j++){
			for(unsigned int i = 1; i<Nx-1; i++){
				v[i+Nx*j] = -(s[i+1+Nx*j] - 2.0*s[i+Nx*j] + s[i-1+Nx*j])/dx/dx
					    -(s[i+Nx+Nx*j] - 2.0*s[i+Nx*j] + s[i-Nx+Nx*j])/dy/dy;
			}
		}

		// Calculation of Interior Vorticity at time t + dt
		for(unsigned int j = 1; j<Ny-1; j++){
			for(unsigned int i = 1; i<Nx-1; i++){
				v_new[i+Nx*j] = v[i+Nx*j] - dt*(s[i+Nx*j+Nx]-s[i+Nx*j-Nx])*(v[i+Nx*j+1]-v[i+Nx*j-1])/4.0/dx/dy + dt*(s[i+Nx*j+1]-s[i+Nx*j-1])*(v[i+Nx*j+Nx]-v[i+Nx*j-Nx])/4.0/dx/dy 
					+ dt/Re*( (v[i+Nx*j+1] - 2.0*v[i+Nx*j] + v[i+Nx*j-1])/dx/dx
						+ (v[i+Nx*j+Nx] - 2.0*v[i+Nx*j] + v[i+Nx*j-Nx])/dy/dy);
			}
		}
		
		
		// Solution of Poisson Problem to Compute Streamfunction at t + dt	
		// Mapping Global Nodes to inner Nodes
		for(int j = 0; j< Ny-2; j++){
			for(int i = 0; i<Nx-2; i++){
				rhs[i+j*(Nx-2)] = v_new[(i+1)+(j+1)*Nx];
			}
		}
		
		// Solving Using Forward Substitution
		F77NAME(dpbtrs) ('U', internal_nodes, ku, 1, a_banded, ku+1, rhs, internal_nodes, info);

		// Mapping Solution to Global Vector
		for(unsigned int j = 0; j<Ny-2; j++){
			for(unsigned int i = 0; i<Nx-2; i++){
				s_new[(i+1)+(j+1)*Nx] = rhs[i+j*(Nx-2)];

			}
		}
	// Calculating 2-norm of solution difference
	norm = 0.0;
	for(unsigned int i = 0 ; i<Nx*Ny; i++){
		norm += (s_new[i] - s[i])*(s_new[i]-s[i]);
	}
	cout << "Info: " << info << endl;
	cout << "2-Norm: " << sqrt(norm) << endl;
	// Copying streamfunction at n+1 back into n for next time iteration 
	cblas_dcopy(Nx*Ny, s_new, 1, s, 1);
	t_elapse+= dt;
	cout << "Time elapsed: " << t_elapse<< endl<< endl;
	}
}

void LidDrivenCavity::ExportSol(){
	ofstream sOut("streamfunction.txt", ios::out | ios::trunc);
	sOut.precision(5);
	sOut << setw(15) << "x"
	     << setw(15) << "y"
	     << setw(15) << "psi" << endl;
	for(int j=0; j<Ny; j++){
		for(int i=0; i<Nx; i++){
			sOut << setw(15) << i*double(Lx/(Nx-1.0))
			     << setw(15) << j*double(Ly/(Ny-1.0))
			     << setw(15) << s[i+Nx*j] << endl;
		}
		sOut << endl;
	}

	ofstream vOut("vorticity.txt", ios::out | ios::trunc);
	vOut.precision(5);
	vOut << setw(15) << "x"
	     << setw(15) << "y"
	     << setw(15) << "vorticity" << endl;
	for(int j=0; j<Ny; j++){
		for(int i=0; i<Nx; i++){
			vOut << setw(15) << i*double(Lx/(Nx-1.0))
			     << setw(15) << j*double(Ly/(Ny-1.0))
			     << setw(15) << v[i+Nx*j] << endl;
		}
		vOut << endl;
	}
	vOut.close();
}
