#pragma once

#include <string>
using namespace std;

class LidDrivenCavity
{
public:
    LidDrivenCavity();
    ~LidDrivenCavity();

    void SetDomainSize(double xlen, double ylen);
    void SetGridSize(int nx, int ny);
    void SetTimeStep(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);
    void MatPrint(double* x, int n);

    void Initialise();
    void Integrate();
    void ExportSol();
    // Add any other public functions

private:
    double* v = nullptr;
    double* s = nullptr;
    double* v_new = nullptr;
    double* s_new = nullptr;

    double dt;
    double T;
    int    Nx;
    int    Ny;
    double Lx;
    double Ly;
    double Re;
};

