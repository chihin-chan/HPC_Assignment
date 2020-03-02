#include "LidDrivenCavity.h"
#include <iostream>
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
}

void LidDrivenCavity::Integrate()
{
}
