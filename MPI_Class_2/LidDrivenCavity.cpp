#include "LidDrivenCavity.h"
#include <iostream>
#include <cstring>
#include <cblas.h>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <vector>
#include <chrono>
#include <thread>

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
	MPI_Barrier(MPI_COMM_WORLD);
	delete[] v;
	delete[] s;
	delete[] v_new;
	delete[] rhs;
	delete [] outbuf_s_D;
	delete [] outbuf_v_D;
	delete [] inbuf_s_D;
	delete [] inbuf_v_D;

	delete [] outbuf_s_U;
	delete [] outbuf_v_U;
	delete [] inbuf_s_U;
	delete [] inbuf_v_U;

	delete [] outbuf_s_L;
	delete [] outbuf_v_L;
	delete [] inbuf_s_L;
	delete [] inbuf_v_L;

	delete [] outbuf_s_R;
	delete [] outbuf_v_R;
	delete [] inbuf_s_R;
	delete [] inbuf_v_R;


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

void LidDrivenCavity::SMatPrintRank(int r){
	if (rank == r){
		cout << "Printing Streamfunction from rank: " << r << endl;
		for(int j = loc_ny-1; j>=0;  j--){
		    for(int i = 0; i<loc_nx; i++){
		        cout << setw(12) << s[i+j*loc_nx];
		    }
		cout << endl;
		}
	}
}

void LidDrivenCavity::VMatPrintRank(int r){
	if (rank == r){
		cout << "Printing Vorticity from rank: " << r << endl;
		for(int j = loc_ny-1; j>=0;  j--){
		    for(int i = 0; i<loc_nx; i++){
		        cout << setw(12) << v[i+j*loc_nx];
		    }
		cout << endl;
		}
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
		if(Ny%Px != 0){
		    loc_ny = (Ny % int( ceil(double(Ny)/double(Px)) )) + 2;
		}
		else{
		    loc_ny = Ny/Px + 2;
		}
	}	
	v = new double[loc_nx*loc_ny];
	s = new double[loc_nx*loc_ny];
	v_new = new double[loc_nx*loc_ny];
	outbuf_s_D = new double[loc_nx];
	outbuf_v_D = new double[loc_nx];
	inbuf_s_D = new double[loc_nx];
	inbuf_v_D = new double[loc_nx];

	outbuf_s_U = new double[loc_nx];
	outbuf_v_U = new double[loc_nx];
	inbuf_s_U = new double[loc_nx];
	inbuf_v_U = new double[loc_nx];

	outbuf_s_L = new double[loc_ny];
	outbuf_v_L = new double[loc_ny];
	inbuf_s_L = new double[loc_ny];
	inbuf_v_L = new double[loc_ny];

	outbuf_s_R = new double[loc_ny];
	outbuf_v_R = new double[loc_ny];
	inbuf_s_R = new double[loc_ny];
	inbuf_v_R = new double[loc_ny];

	// Initialising Vorticity, w and Streamfunction, s to zero
	fill_n(s, loc_nx*loc_ny, 0.0);
	fill_n(v, loc_nx*loc_ny, 0.0);
	fill_n(v_new, loc_nx*loc_ny, 0.0);
	fill_n(outbuf_s_U, loc_nx, 0.0);
	fill_n(outbuf_v_U, loc_nx, 0.0);
	fill_n(inbuf_s_U, loc_nx, 0.0);
	fill_n(inbuf_v_U, loc_nx, 0.0);
	fill_n(outbuf_s_D, loc_nx, 0.0);
	fill_n(outbuf_v_D, loc_nx, 0.0);
	fill_n(inbuf_s_D, loc_nx, 0.0);
	fill_n(inbuf_v_D, loc_nx, 0.0);
	fill_n(outbuf_s_R, loc_ny, 0.0);
	fill_n(outbuf_v_R, loc_ny, 0.0);
	fill_n(inbuf_s_R, loc_ny, 0.0);
	fill_n(inbuf_v_R, loc_ny, 0.0);
	fill_n(outbuf_s_L, loc_ny, 0.0);
	fill_n(outbuf_v_L, loc_ny, 0.0);
	fill_n(inbuf_s_L, loc_ny, 0.0);
	fill_n(inbuf_v_L, loc_ny, 0.0);
	// Calculating Grid Spacing
	dx = double(Lx) / double((Nx-1.0));
	dy = double(Ly) / double((Ny-1.0));
	MPI_Barrier(MPI_COMM_WORLD);
	
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
            v[i] = (s[i] - s[i+loc_nx])*2.0/dy/dy;
        }
    }
}

void LidDrivenCavity::Communicate()
{	
	MPI_Request req1, req2, req3, req4, req5, req6, req7, req8,
				req9, req10, req11, req12, req13, req14, req15, req16;
	// Sending and Receiving from DOWN Neighbours
	if(nghbrs[DOWN] != -2){
		// Packing
		int k = 0;
		for (int i=loc_nx; i<2*loc_nx; i++){
			outbuf_s_D[k] = s[i];
			outbuf_v_D[k] = v[i];
			k++;
		}
		// Sending
		// cout << "I'm rank: " << rank << " and I am sending to DOWN rank: " << nghbrs[DOWN] << endl;
		MPI_Isend(outbuf_s_D, loc_nx, MPI_DOUBLE, nghbrs[DOWN], tag[sf], mygrid, &req1);
		MPI_Isend(outbuf_v_D, loc_nx, MPI_DOUBLE, nghbrs[DOWN], tag[vort], mygrid, &req2);
		
		// Receiving
		// cout << "I'm rank: " << rank << " and I am recv from DOWN rank: " << nghbrs[DOWN] << endl;
		MPI_Irecv(inbuf_s_D, loc_nx, MPI_DOUBLE, nghbrs[DOWN], tag[sf], mygrid,&req3);
		MPI_Irecv(inbuf_v_D, loc_nx, MPI_DOUBLE, nghbrs[DOWN], tag[vort], mygrid,&req4);

		// Wait
		 MPI_Wait(&req1, MPI_STATUS_IGNORE);
		 MPI_Wait(&req2, MPI_STATUS_IGNORE);
		 MPI_Wait(&req3, MPI_STATUS_IGNORE);
		 MPI_Wait(&req4, MPI_STATUS_IGNORE);


		// Unpacking
		k = 0;
		for (int i = 0; i< loc_nx; i++){
			s[i] = inbuf_s_D[k];
			v[i] = inbuf_v_D[k];
			k++;
		}
	}

	// Sending and Receiving from UP Neighbours
	if(nghbrs[UP] != -2){
		// Packing
		int k = 0;
		for (int i = loc_nx*(loc_ny-1)-loc_nx; i<loc_nx*(loc_ny-1); i++){
			outbuf_s_U[k] = s[i];
			outbuf_v_U[k] = v[i];
			k++;
		}
		// Sending
		// cout << "I'm rank: " << rank << " and I am sending to UP rank: " << nghbrs[UP] << endl;
		MPI_Isend(outbuf_s_U, loc_nx, MPI_DOUBLE, nghbrs[UP], tag[sf], mygrid, &req5);
		MPI_Isend(outbuf_v_U, loc_nx, MPI_DOUBLE, nghbrs[UP], tag[vort], mygrid, &req6);
		// Receiving
		// cout << "I'm rank: " << rank << " and I am recv from UP rank: " << nghbrs[UP] << endl;
		MPI_Irecv(inbuf_s_U, loc_nx, MPI_DOUBLE, nghbrs[UP], tag[sf], mygrid,&req7);
		MPI_Irecv(inbuf_v_U, loc_nx, MPI_DOUBLE, nghbrs[UP], tag[vort], mygrid,&req8);
		// Wait
		MPI_Wait(&req5, MPI_STATUS_IGNORE);
		MPI_Wait(&req6, MPI_STATUS_IGNORE);
		MPI_Wait(&req7, MPI_STATUS_IGNORE);
		MPI_Wait(&req8, MPI_STATUS_IGNORE);


		// Unpacking
		k = 0;
		for(int i = (loc_nx)*(loc_ny-1); i<(loc_nx*loc_ny); i++){
			s[i] = inbuf_s_U[k];
			v[i] = inbuf_v_U[k];
			k++;
		}
	}

	// Sending and Receiving from RIGHT Neighbours
	if(nghbrs[RIGHT] != -2){
		int k = 0;
		// Packing
		for (int i = loc_nx - 2; i < loc_nx*loc_ny-1; i+=loc_nx){
			outbuf_s_R[k] = s[i];
			outbuf_v_R[k] = v[i];
			k++;
		}
		// Sending
		// cout << "I'm rank: " << rank << " and I am sending to RIGHT rank: " << nghbrs[RIGHT] << endl;
		MPI_Isend(outbuf_s_R, loc_ny, MPI_DOUBLE, nghbrs[RIGHT], tag[sf], mygrid, &req9);
		MPI_Isend(outbuf_v_R, loc_ny, MPI_DOUBLE, nghbrs[RIGHT], tag[vort], mygrid, &req10);
		// Receiving
		// cout << "I'm rank: " << rank << " and I am recv to RIGHT rank: " << nghbrs[RIGHT] << endl;
		MPI_Irecv(inbuf_s_R, loc_ny, MPI_DOUBLE, nghbrs[RIGHT], tag[sf], mygrid,&req11);
		MPI_Irecv(inbuf_v_R, loc_ny, MPI_DOUBLE, nghbrs[RIGHT], tag[vort], mygrid,&req12);
		// Wait
		MPI_Wait(&req9, MPI_STATUS_IGNORE);
		MPI_Wait(&req10, MPI_STATUS_IGNORE);
		MPI_Wait(&req11, MPI_STATUS_IGNORE);
		MPI_Wait(&req12, MPI_STATUS_IGNORE);

		k = 0;
		// Unpacking
		for (int i = loc_nx -1; i<(loc_nx*loc_ny); i+= loc_nx){
			s[i] = inbuf_s_R[k];
			v[i] = inbuf_v_R[k];
			k++;
		}
	}

	// Sending and Receiving from LEFT Neighbours
	if(nghbrs[LEFT] != -2){
		// Packing
		int k = 0;
		for (int i = 1; i<(loc_nx*loc_ny); i+=loc_nx){
			outbuf_s_L[k] = s[i];
			outbuf_v_L[k] = v[i];
			k++;
		}
		// Sending
		// cout << "I'm rank: " << rank << " and I am sending to LEFT rank: " << nghbrs[LEFT] << endl;
		MPI_Isend(outbuf_s_L, loc_ny, MPI_DOUBLE, nghbrs[LEFT], tag[sf], mygrid, &req13);
		MPI_Isend(outbuf_v_L, loc_ny, MPI_DOUBLE, nghbrs[LEFT], tag[vort], mygrid, &req14);
		// Receiving
		// cout << "I'm rank: " << rank << " and I am recv to LEFT rank: " << nghbrs[LEFT] << endl;
		MPI_Irecv(inbuf_s_L, loc_ny, MPI_DOUBLE, nghbrs[LEFT], tag[sf], mygrid,&req15);
		MPI_Irecv(inbuf_v_L, loc_ny, MPI_DOUBLE, nghbrs[LEFT], tag[vort], mygrid,&req16);
		// Wait
		MPI_Wait(&req13, MPI_STATUS_IGNORE);
		MPI_Wait(&req14, MPI_STATUS_IGNORE);
		MPI_Wait(&req15, MPI_STATUS_IGNORE);
		MPI_Wait(&req16, MPI_STATUS_IGNORE);
		// Unpacking
		k = 0;			
		for(int i = 0; i < (loc_nx*loc_ny) -1; i+=loc_nx){
			s[i] = inbuf_s_L[k];
			v[i] = inbuf_v_L[k];
			k++;
		}
	}
}

void LidDrivenCavity::InteriorUpdate(){
	
	int xstart = 1;
	int jstart = 1;
	int xend = 1;
	int jend = 1;
	if(nghbrs[UP] == -2){
		jend = 2;
	}
	if(nghbrs[DOWN] == -2){
		jstart = 2;
	}
	if(nghbrs[RIGHT] == -2){
		xend = 2;
	}
	if(nghbrs[LEFT] == -2){
		xstart = 2;
	}
	// Calculation of interior vorticity at time t
	for(int i = xstart; i<loc_nx-xend; i++){
		for(int j = jstart; j<loc_ny-jend; j++){
			v[i+loc_nx*j] = -(s[i+loc_nx*j+1] - 2.0*s[i+loc_nx*j] + s[i+loc_nx*j-1])/dx/dx
					-(s[i+loc_nx*j+loc_nx] -2.0*s[i+loc_nx*j] + s[i+loc_nx*j-loc_nx])/dy/dy;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	Communicate();


	// Calculation of interior vorticity at time t+dt;
	for(int i = xstart; i<loc_nx-xend; i++){
		for(int j = jstart; j<loc_ny-jend; j++){
			v_new[i+loc_nx*j] = v[i+loc_nx*j] 
         - dt*(s[i+loc_nx*j+loc_nx]-s[i+loc_nx*j-loc_nx])*(v[i+loc_nx*j+1]-v[i+loc_nx*j-1])/4.0/dx/dy
         + dt*(s[i+loc_nx*j+1]-s[i+loc_nx*j-1])*(v[i+loc_nx*j+loc_nx]-v[i+loc_nx*j-loc_nx])/4.0/dx/dy
         + dt*( (v[i+loc_nx*j+1] - 2.0*v[i+loc_nx*j] + v[i+j*loc_nx-1])/dx/dx 
                + (v[i+loc_nx*j+loc_nx] - 2.0*v[i+loc_nx*j] + v[i+loc_nx*j-loc_nx])/dy/dy )/Re;
		}
	}
	for(int i = xstart; i<loc_nx-xend; i++){
		for(int j = jstart; j<loc_ny-jend; j++){
			v[i+loc_nx*j] = v_new[i+loc_nx*j];
		}
	}
}
	
void LidDrivenCavity::MapRHS(){
	int x_off = 2;
	int y_off = 2;
	int x_start = 1;
	int y_start = 1;
	if (nghbrs[UP] == -2){
		y_off += 1;
	}
	if (nghbrs[DOWN] == -2){
		y_off += 1;
		y_start = 2;
	}
	if (nghbrs[RIGHT] == -2){
		x_off += 1;
	}
	if (nghbrs[LEFT] == -2){
		x_off += 1;
		x_start = 2;
	}

	for(int i = 0; i < loc_nx-x_off; i++){
		for(int j = 0; j < loc_ny-y_off; j++){
			rhs[i+j*(loc_nx-x_off)] = v[(i+x_start) + loc_nx*(j+y_start)];
            // Subtracting streamfunction on left boundary
			if (i==0){
				rhs[i+j*(loc_nx-x_off)] = rhs[i+j*(loc_nx-x_off)] 
				+ s[(i+x_start) + loc_nx*(j+y_start) - 1]/dx/dx;
				
			}
			// Subtracting streamfunction on bottom boundary
			if (j==0){
				rhs[i+j*(loc_nx-x_off)] = rhs[i+j*(loc_nx-x_off)] 
				+ s[(i+x_start) + loc_nx*(j+y_start) - loc_nx]/dy/dy;
			}
			// Subtracting streamfunction on right boundary
			if (i == loc_nx-x_off-1){
				rhs[i+j*(loc_nx-x_off)] = rhs[i+j*(loc_nx-x_off)] 
				+ s[(i+x_start) + loc_nx*(j+y_start) + 1]/dx/dx;
			}	
			// Subtracting streamfunction on top boundary
			if (j == loc_ny-y_off-1){
				rhs[i+j*(loc_nx-x_off)] = rhs[i+j*(loc_nx-x_off)] 
				+ s[(i+x_start) + loc_nx*(j+y_start) + loc_nx]/dy/dy;
			}
		}
	}
}


void LidDrivenCavity::iMapRHS(){
	int x_off = 2;
	int y_off = 2;
	int x_start = 1;
	int y_start = 1;
	if (nghbrs[UP] == -2){
		y_off += 1;
	}
	if (nghbrs[DOWN] == -2){
		y_off += 1;
		y_start = 2;
	}
	if (nghbrs[RIGHT] == -2){
		x_off += 1;
	}
	if (nghbrs[LEFT] == -2){
		x_off += 1;
		x_start = 2;
	}
	for(int i = 0; i < loc_nx-x_off; i++){
		for(int j = 0; j < loc_ny-y_off; j++){
			s[(i+x_start) + loc_nx*(j+y_start)] = rhs[i+j*(loc_nx-x_off)];
		}
	}
}

void LidDrivenCavity::Integrate()
{	

	int internal_nodes;
	int ku;
	cout << "Rank: " << rank << "	loc_nx: " << loc_nx << "	loc_ny: " << loc_ny << endl;
	// Initialising properties of banded matrix A, to solve Ax = b
	int x_off = 2;
	int y_off = 2;
	if (nghbrs[UP] == -2){
		y_off += 1;
	}
	if (nghbrs[DOWN] == -2){
		y_off += 1;
	}
	if (nghbrs[RIGHT] == -2){
		x_off += 1;
	}
	if (nghbrs[LEFT] == -2){
		x_off += 1;
	}
	internal_nodes = (loc_nx-x_off)*(loc_ny-y_off);
	ku = loc_nx-x_off;			  // No. of super diagonals
	int k = (ku+1)*internal_nodes;	  // Size of banded matrix (ku+kl+1)*N
	double alpha = 2.0*(1.0/dx/dx + 1.0/dy/dy); // Coefficients of i,j
	double beta_x = -1.0/dx/dx;		// Coefficients of i+/-1, j
	double beta_y = -1.0/dy/dy;		// Coefficients of i/, j+/-1
	double norm;				// Storing 2-norm of solution difference

	int info;
	int count = 1;
	// a_banded holds matrix A in banded format
	double* a_banded = new double[k];
	fill_n(a_banded, k, 0.0);	
	// rhs stores vector b for Ax = b
	rhs = new double[internal_nodes];
	fill_n(rhs, internal_nodes, 0.0);
	// Try and Catch if dt is too large
	double cond = Re*Lx/(Nx-1)*Ly/(Ny-1) / 4;
	try{	
		if(dt >= cond){throw std::out_of_range("");}
	}
	catch(std::out_of_range const &e){
		if(rank == 0){		
		cout << "dt: " << dt << 
		", is too large and should be lesser than: " << cond << endl;
		}
	      return;	      
	}

	// Storing Elements of a_banded in banded format
	a_banded[ku] = alpha;
	for(int i=2*(ku+1)-1; i<k; i+=(ku+1)){
		a_banded[i] = alpha;
		// Storing beta_x
		if(count % (loc_nx-x_off) != 0){
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
	
/*
	// Printing A_banded for checking
    if (rank == 0){
        for(unsigned int i = 0; i < (ku+1) ; i++){
            for(unsigned int j = 0; j < internal_nodes; j++){
                cout << a_banded[i+j*(ku+1)] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
*/
	// Caching Cholesky factorisation
	F77NAME(dpbtrf) ('u', internal_nodes, ku, a_banded, ku+1, info);

	// Waiting for all processes before starting time-loop
		
	cout << "Rank: " << rank << "  starting time loop" << endl;	
	// Starting time loop
	double t_elapse = 0.0;
	while (t_elapse < T){	
		// Imposing Boundary Conditions
		BoundaryConditions();

		// Calculation of Interior Vorticity at time t AND t+dt
		InteriorUpdate();

        	// Solution of Poisson Problem to Compute Streamfunction at t + dt	
		for (int i = 0; i < 5; i++){
			
			// Mapping inner vorticity and boundaries of streamfunction to RHS at t+dt
			MapRHS();

			// Solving
			F77NAME(dpbtrs) ('U', internal_nodes, ku, 1, a_banded, ku+1, rhs, internal_nodes, info);

			// Mapping rhs to inner streamfunction
			iMapRHS();

			// Updating streamfunction to neighbours
			MPI_Barrier(MPI_COMM_WORLD);
			Communicate();
		}
		
		if(rank==0) cout << endl << endl << "Time-step: " << t_elapse << endl << endl;
		t_elapse += dt;
		MPI_Barrier(MPI_COMM_WORLD);		
	}
	VMatPrintRank(0);
}

void LidDrivenCavity::ExportSol(){
	
	
	int start;
	int end;
	ofstream sOut("streamfunction.txt", ios::out | ios::trunc);
   	sOut.precision(5);
	// Iterating down domain rows
	for(int i=0; i<Px-1; i++){
		// Iterating through subdomain rows
		for(int leny = 0; leny<loc_ny-2; leny++){
			// Define starting position of each subdomain
			start = loc_nx*(loc_ny-1) - loc_nx + 1 - leny*loc_nx;
			end = (loc_nx*loc_ny)-1 - loc_nx - leny*loc_nx;
			// Iterating through subdomain cols
			for(int j = 0; j<Py-1; j++){
				// Checking if rank coincide
				if(coords[0] == i && coords[1] == j){
					// Iterating through cols in each subdomain
					for(int start; start<end; start++){
						sOut << setw(15) << s[start] <<"," << setw(15);
					}
				}
			}
		sOut << endl;
		}
	}
	sOut.close();
}
