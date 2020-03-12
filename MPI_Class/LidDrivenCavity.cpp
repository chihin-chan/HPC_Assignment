#include "LidDrivenCavity.h"
#include <iostream>
#include <cstring>
#include <cblas.h>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <mpi.h>

#define UP    0
#define DOWN  1
#define LEFT  2
#define RIGHT 3
#define sf	  0
#define vort  1
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
}

void LidDrivenCavity::SetGridSize(int nx, int ny)
{
	Nx = nx;
	Ny = ny;

}

void LidDrivenCavity::SetPartitionSize(int ppx, int ppy)
{
	Px = ppx;
	Py = ppy;

}
void LidDrivenCavity::SetTimeStep(double deltat)
{
	dt = deltat;
}

void LidDrivenCavity::SetFinalTime(double finalt)
{
	T = finalt;
}

void LidDrivenCavity::SetReynoldsNumber(double re)
{
	Re = re;
}

void LidDrivenCavity::GetRank(int rr)
{
	rank = rr;
}

void LidDrivenCavity::GetSize(int ss)
{
	size = ss;
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

	// Generating a communicator for 2-D Cartesian Topology
	dims[0] = Px;
	dims[1] = Py;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &mygrid);

	// Retreving coordinates of Rank
	MPI_Cart_coords(mygrid, rank, 2, coords);

    // Obtain Shifted Source and Destination ranks in both directions
    MPI_Cart_shift(mygrid, 0, 1, &nghbrs[UP], &nghbrs[DOWN]);
    MPI_Cart_shift(mygrid, 1, 1, &nghbrs[LEFT], &nghbrs[RIGHT]);
	
	// Printing Initialised Variables
	if(rank == 0 ){
		cout << "Domain Size, Lx: " << Lx <<"	Ly: " << Ly << endl;
		cout << "Grid Size, Nx: " << Nx <<"	Ny: " << Ny << endl;
		cout << "Partition Size, Px: " << Px << " Py: " << Py << endl;
		cout << "Final Time: " << T << endl;
		cout << "Time Step: " << dt << endl;
		cout << "Reynolds Number: " << Re << endl;
	}
	// Defining vorticity and streafunction matrix for subdomains
    loc_nx = ceil(double(Nx)/double(Py)) + 2;
    loc_ny = ceil(double(Ny)/double(Px)) + 2;
	if(coords[1] == Py-1){
    	if(Nx%Py != 0){
        	loc_nx = (Nx % int( ceil(double(Nx)/double(Py)) )) + 2;
    	}
    	else{
        	loc_nx = Nx/Py + 2;
    	}
	}
	if(coords[0] == Px-1){
		if(Ny%Py != 0){
		    loc_ny = (Ny % int( ceil(double(Ny)/double(Px)) )) + 2;
		}
		else{
		    loc_ny = Ny/Px + 2;
		}
	}	
	v = new double[loc_nx*loc_ny];
	s = new double[loc_nx*loc_ny];
	// Initialising Vorticity, w and Streamfunction, s to zero
	fill_n(s, loc_nx*loc_ny, 0.0);
	fill_n(v, loc_nx*loc_ny, 0.0);
	// Calculating Grid Spacing
	dx = double(Lx) / double((Nx-1.0));
	dy = double(Ly) / double((Ny-1.0));
	
}

void LidDrivenCavity::BoundaryConditions()
{
	double U = 1.0;
	///////////////
    // TOP DOMAINS
    ///////////////
    if(nghbrs[UP] == -2){
		for(int i = loc_nx*(loc_ny-1)-loc_nx+1; i<(loc_nx*loc_ny)-1-loc_nx; i++){
			v[i] = (s[i] - s[i-loc_nx])*2.0/dy/dy - 2.0*U/dy;
		}
	}

    ////////////////
    // LEFT DOMAINS
    ////////////////
    if(nghbrs[LEFT] == -2){
		for(int i=1+loc_nx; i<loc_nx*(loc_ny-1); i+=loc_nx){
			v[i] = (s[i] - s[i+1])*2.0/dx/dx;
        }
    }
    /////////////////////////////////////////
    // RIGHT DOMAINS (Potentially Deficient)
    /////////////////////////////////////////
    if(nghbrs[RIGHT] == -2){
		for(int i = (2*loc_nx) - 2; i < (loc_nx*loc_ny) - 2; i += loc_nx){
			v[i] = (s[i] - s[i-1])*2.0/dx/dx;
        }
    }
    
    ////////////////////////////////////////
    // BOTTOM DOMAINS (Potentially Deficient)
    ///////////////////////////////////////
    if(nghbrs[DOWN] == -2){
        for(int i = 1+loc_nx; i < (2*loc_nx) - 1; i++){
            v[i] = (s[i] - s[i+loc_nx])*2/dy/dy;
        }
    }
}

void LidDrivenCavity::Communicate()
{	
	double* inbuf_s;
	double* inbuf_v;
	double* outbuf_s;
	double* outbuf_v;
	
	// Sending and Receiving from DOWN Neighbours
	if(nghbrs[DOWN] != -2){
		outbuf_s = new double[loc_nx];
		outbuf_v = new double[loc_nx];
		inbuf_s = new double[loc_nx];
		inbuf_v = new double[loc_nx];
		// Packing
		int k = 0;
		for (int i=loc_nx; i<2*loc_nx; i++){
			outbuf_s[k] = s[i];
			outbuf_v[k] = v[i];
			k++;
		}
		// Sending
		cout << "I'm rank: " << rank << " and I am sending to DOWN rank: " << nghbrs[DOWN] << endl;
		MPI_Isend(outbuf_s, loc_nx, MPI_DOUBLE, nghbrs[DOWN], tag[sf], mygrid, &req);
		MPI_Isend(outbuf_v, loc_nx, MPI_DOUBLE, nghbrs[DOWN], tag[vort], mygrid, &req);
		// Receiving
		cout << "I'm rank: " << rank << " and I am recv to DOWN rank: " << nghbrs[DOWN] << endl;
		MPI_Irecv(inbuf_s, loc_nx, MPI_DOUBLE, nghbrs[DOWN], tag[sf], mygrid, &req);
		MPI_Irecv(inbuf_v, loc_nx, MPI_DOUBLE, nghbrs[DOWN], tag[sf], mygrid, &req);
		// Unpacking
		k = 0 ;
		for (int i = 0; i< loc_nx; i++){
			s[i] = inbuf_s[k];
			v[i] = inbuf_v[k];
			k++;
		}
		delete [] outbuf_s;
		delete [] outbuf_v;
		delete [] inbuf_s;
		delete [] inbuf_v;
	}
	// Sending and Receiving from UP Neighbours
	if(nghbrs[UP] != -2){
		inbuf_s = new double[loc_nx];
		inbuf_v = new double[loc_nx];
		outbuf_s = new double[loc_nx];
		outbuf_v = new double[loc_nx];
		// Receiving
		cout << "I'm rank: " << rank << " and I am recv to UP rank: " << nghbrs[UP] << endl;
		MPI_Irecv(inbuf_s, loc_nx, MPI_DOUBLE, nghbrs[UP], tag[sf], mygrid, &req);
		MPI_Irecv(inbuf_v, loc_nx, MPI_DOUBLE, nghbrs[UP], tag[vort], mygrid, &req);
		// Unpacking
		int k = 0;
		for(int i = (loc_nx)*(loc_ny-1); i<(loc_nx*loc_ny)-1; i++){
			s[i] = inbuf_s[k];
			v[i] = inbuf_v[k];
			k++;
		}
		// Packing
		k = 0;
		for (int i = loc_nx*(loc_ny-1)-loc_nx; i<(loc_nx)*(loc_ny-1); i++){
			outbuf_s[k] = s[i];
			outbuf_v[k] = v[i];
			k++;
		}
		// Sending
		cout << "I'm rank: " << rank << " and I am sending to UP rank: " << nghbrs[UP] << endl;
		MPI_Isend(outbuf_s, loc_nx, MPI_DOUBLE, nghbrs[UP], tag[sf], mygrid, &req);
		MPI_Isend(outbuf_v, loc_nx, MPI_DOUBLE, nghbrs[UP], tag[vort], mygrid, &req);
		delete [] outbuf_s;
		delete [] outbuf_v;
		delete [] inbuf_s;
		delete [] inbuf_v;
	}

	// Sending and Receiving from RIGHT Neighbours
	if(nghbrs[RIGHT] != -2){
		outbuf_s = new double[loc_ny];
		outbuf_v = new double[loc_ny];
		inbuf_s = new double[loc_ny];
		inbuf_v = new double[loc_ny];
		int k = 0;
		// Packing
		for (int i = loc_nx - 2; i < loc_nx*loc_ny-1; i+=loc_nx){
			outbuf_s[k] = s[i];
			outbuf_v[k] = v[i];
			k++;
		}
		// Sending
		cout << "I'm rank: " << rank << " and I am sending to RIGHT rank: " << nghbrs[RIGHT] << endl;
		MPI_Isend(outbuf_s, loc_ny, MPI_DOUBLE, nghbrs[RIGHT], tag[sf], mygrid, &req);
		MPI_Isend(outbuf_v, loc_ny, MPI_DOUBLE, nghbrs[RIGHT], tag[vort], mygrid, &req);
		// Receiving
		cout << "I'm rank: " << rank << " and I am recv to RIGHT rank: " << nghbrs[RIGHT] << endl;
		MPI_Irecv(inbuf_s, loc_ny, MPI_DOUBLE, nghbrs[RIGHT], tag[sf], mygrid, &req);
		MPI_Irecv(inbuf_v, loc_ny, MPI_DOUBLE, nghbrs[RIGHT], tag[vort], mygrid, &req);
		k = 0;
		// Unpacking
		for (int i = loc_nx -1; i<(loc_nx*loc_ny)-1; i+= loc_nx){
			s[i] = inbuf_s[k];
			v[i] = inbuf_v[k];
			k++;
		}
		delete [] outbuf_s;
		delete [] outbuf_v;
		delete [] inbuf_s;
		delete [] inbuf_v;
	}
	
	// Sending and Receiving from LEFT Neighbours
	if(nghbrs[LEFT] != -2){
		inbuf_s = new double[loc_ny];
		inbuf_v = new double[loc_ny];
		outbuf_s = new double[loc_ny];
		outbuf_v = new double[loc_ny];
		// Receiving
		cout << "I'm rank: " << rank << " and I am recv to LEFT rank: " << nghbrs[LEFT] << endl;
		MPI_Irecv(inbuf_s, loc_ny, MPI_DOUBLE, nghbrs[LEFT], tag[sf], mygrid, &req);
		MPI_Irecv(inbuf_v, loc_ny, MPI_DOUBLE, nghbrs[LEFT], tag[vort], mygrid, &req);
		// Unpacking
		int k = 0;			
		for(int i = 0; i<(loc_nx*loc_ny) -1; i+=loc_nx){
			s[i] = inbuf_s[k];
			v[i] = inbuf_v[k];
			k++;
		}
		// Packing
		k = 0;
		for (int i = 1; i<(loc_nx*loc_ny) - 1; i+=loc_nx){
			outbuf_s[k] = s[i];
			outbuf_v[k] = v[i];
		}
		// Sending
		cout << "I'm rank: " << rank << " and I am sending to LEFT rank: " << nghbrs[LEFT] << endl;
		MPI_Isend(outbuf_s, loc_ny, MPI_DOUBLE, nghbrs[LEFT], tag[sf], mygrid, &req);
		MPI_Isend(outbuf_v, loc_ny, MPI_DOUBLE, nghbrs[LEFT], tag[vort], mygrid, &req);
		delete [] outbuf_s;
		delete [] outbuf_v;
		delete [] inbuf_s;
		delete [] inbuf_v;
	}
	MPI_Barrier(mygrid);
}

void InteriorUpdate(){
	//////////////////
	// NOT DONE
	//////////////
	// Calculation of interior vorticity at time t
	for(int i = xstart; i<loc_nx-xend; i++){
		for(int j = jstart; j<loc_ny-jend; j++){
			v[i+nx*j] = -(s[i+loc_nx*j+1] - 2.0*s[i+loc_nx*j] + s[i+loc_nx*j-1])/dx/dx
						-(s[i+loc_nx*j+loc_nx] -2.0*s[i+loc_nx*j] + s[i+loc_nx*j-nx])/dy/dy;
		}
	}
	// Calculation of interior vorticity at time t+dt;
	for(int i = xstart; i<loc_nx-xend; i++){
		for(int j = jstart; j<loc_ny-jend; j++){
			v[i+nx*j] = dt/Re*( (v[i+loc_nx*j+1] - 2.0*v[i+loc_nx*j] + v[i+loc_nx*j-1])/dx/dx
						+ (v[i+loc_nx*j+loc_nx] - 2.0*v[i+loc_nx*j] + v[i+loc_nx*j-loc_nx])/dy/dy)
						+ dt*( (s[i+loc_nx*j+1] - s[i+loc_nx*j-1])/2.0/dx*
							   (v[i+loc_nx*j+loc_nx] - v[i+loc_nx*j-loc_nx])/2.0/dy)
						- dt*( (s[i+loc_nx*j+loc_nx] - s[i+loc_nx*j-loc_nx])/2.0/dy*
							   (v[i+loc_nx*j+1] - v[i+loc_nx*j-1])/2.0/dx) 
						+ v[i+loc_nx*j];
		}
	}
}


void LidDrivenCavity::Integrate()
{	

	int internal_nodes;
	int ku;
	cout << "Rank: " << rank << " my loc_nx: " << loc_nx << "and my loc_ny: " << loc_ny << endl;
	// Initialising properties of banded matrix A, to solve Ax = b
	if (nghbrs[RIGHT] == -2 && nghbrs[LEFT] == -2 && nghbrs[DOWN] == -2 && nghbrs[DOWN] == -2){
		internal_nodes = (loc_nx-4)*(loc_ny-4);
		ku = loc_nx-4;			  
	}
	else if(nghbrs[RIGHT] == -2 && (nghbrs[UP] == -2 || nghbrs[DOWN] == -2) && nghbrs[LEFT] == -2){
		internal_nodes = (loc_nx-4)*(loc_ny-3);
		ku = loc_nx-4;	
	}
	else if(nghbrs[UP] == -2 && nghbrs[DOWN] == -2 && (nghbrs[RIGHT] == -2 || nghbrs[LEFT] == -2)){
		internal_nodes = (loc_nx-3)*(loc_ny-4);
		ku = loc_nx-3;
	}
	else if(nghbrs[UP] == -2 && nghbrs[DOWN] == -2){
		internal_nodes = (loc_nx-2)*(loc_ny-4);
		ku = loc_nx-2;
	}
	else if(nghbrs[RIGHT] == -2 && nghbrs[LEFT] == -2){
		internal_nodes = (loc_nx-4)*(loc_ny-2);
		ku = loc_nx-4;
	}
	else if( (nghbrs[RIGHT] == -2 && nghbrs[UP] == -2) || (nghbrs[LEFT] == -2 && nghbrs[UP] == -2) || (nghbrs[RIGHT] == -2 && nghbrs[DOWN] == -2) || (nghbrs[LEFT] == -2 && nghbrs[DOWN] == -2)){
		internal_nodes = (loc_nx-3)*(loc_ny-3);
		ku = loc_nx-3;
	}
	else if( nghbrs[RIGHT] == -2 || nghbrs[LEFT] == -2){
		internal_nodes = (loc_nx-3)*(loc_ny-2);
		ku = loc_nx - 3;
	}
	else if (nghbrs[UP] == -2 || nghbrs[DOWN] == -2){
		internal_nodes = (loc_nx-2)*(loc_ny-3);
		ku = loc_nx - 2;
	}
	else{
		internal_nodes = (loc_nx-2)*(loc_ny-2);
		ku = loc_nx - 2;
	}
	
	int k = (ku+1)*internal_nodes;	  // Size of banded matrix (ku+kl+1)*N
	double alpha = 2.0*(1.0/dx/dx + 1.0/dy/dy); // Coefficients of i,j
	double beta_x = -1.0/dx/dx;		// Coefficients of i+/-1, j
	double beta_y = -1.0/dy/dy;		// Coefficients of i/, j+/-1
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
		
	double t_elapse = 0.0;
	// Starting time loop
		
		BoundaryConditions();
		Communicate();
		// Calculation of Interior Vorticity at time t


		// Calculation of Interior Vorticity at time t + dt
		
		// Solution of Poisson Problem to Compute Streamfunction at t + dt	
		// Mapping Global Nodes to inner Nodes

		// Solving Using Forward Substitution
		// F77NAME(dpbtrs) ('U', internal_nodes, ku, 1, a_banded, ku+1, rhs, internal_nodes, info);

		// Mapping Solution to Global Vector

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
	sOut.close();
}
