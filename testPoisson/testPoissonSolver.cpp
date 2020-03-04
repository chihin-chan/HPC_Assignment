#include <iomanip>
#include "PoissonSolver.h"
#include <iostream>
#include <string.h>
#include <math.h>

using namespace std;

int main(){

    // Creates a Poisson Solver
    int nx = 120;
    int ny = 120;
    double Lx = M_PI;
    double Ly = M_PI;
    double dx = double(Lx)/(nx-1.0);
    double dy = double(Ly)/(ny-1.0);

    // Initialises Poisson Solver
    PoissonSolver ps(nx,ny,dx,dy);
   
    // Initialise b
    double* b = new double[(nx-2)*(ny-2)];
    for(int j=0; j<(ny-2); j++){
        for(int i=0; i<(nx-2); i++){
            b[i+j*(nx-2)] = sin((i+1)*dx)*sin((j+1)*dy);
        }
    }
    for(int i=0; i<(ny-2)*(nx-2);i++){
        cout << b[i] << endl;
}
cout << endl << endl << endl;
    // Call Solver.solve
    ps.CholSolve(b);
    
    for(int i=0; i<(nx-2)*(ny-2);i++){
    }

    delete[] b;

}