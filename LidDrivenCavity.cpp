#include "LidDrivenCavity.h"
#include <iostream>
#include <cstring>
#include <cblas.h>

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
	cout << "Initialising solution..." << endl;
}

void LidDrivenCavity::Integrate()
{
	double dx = Lx/(Nx-1);
	double dy = Ly/(Ny-1);
	double U = 1.0;
	double t_lapse = 0.0;

	// Starting time loop
	// while t_lapse < T{
	
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

	// Enforcing Left Boundary Conditions for Voticity, w
	for(unsigned int i = Nx; i < Nx*Ny; i+=Nx){
		v[i] = (s[i] - s[i+1])*2.0/dx/dx;
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
			v_new[i+Nx*j] = v[i+Nx*j]*dt - (s[i+Nx*j+Nx] - s[i+Nx*j-Nx])*(v[i+Nx*j+1] - v[i+Nx*j-1])/2.0/2.0/dx/dy 
					+ (s[i+Nx*j+1] - s[i+Nx*j-1])*(v[i+Nx*j+Nx] - v[i+Nx*j-Nx])/2.0/2.0/dx/dy
			       		+ ( (v[i+Nx*j+1] - 2.0* v[i+Nx*j] +v[i+Nx*j-1])/dx/dx + (v[i+Nx*j+Nx] - 2.0* v[i+Nx*j] +v[i+Nx*j-Nx])/dy/dy )/Re;
		}
	}
	
	// Solution of Poisson Problem to Compute Streamfunction at t + dt
	const int internal_nodes = (Nx-2)*(Ny-2);
	const int ku = Nx-2;
	const int k = (ku+1)*internal_nodes;
	double alpha = 2.0*(1.0/dx/dx + 1.0/dy/dy);
	double beta_x = -1/dx/dx;
	double beta_y = -1/dy/dy;

	double* a_banded = new double[k];
	cout << "k = (ku+1)*N = " << k << endl;	
	
	// Filling Banded Matrix a for conjugate gradient
	a_banded[ku] = alpha;
	for(unsigned int i=ku+ku+1; i<k; i+=(Nx-1)){
		a_banded[i] = alpha;
		a_banded[i-1] = beta_x;
	}
	for(unsigned int i=ku*(ku+1); i<k; i+=(Nx-1)){
		a_banded[i] = beta_y;
	}
	for(unsigned int i=0; i<k;i++){
		cout<< a_banded[i] << endl;
	}
	
	// Solving using conjugate gradient method
	cblas_dcopy(n

	


}
