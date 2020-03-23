#include "LidDrivenCavity.h"
#include "PoissonSolver.h"
#include <iostream>
#include <cstring>
#include <cblas.h>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <thread>

#include <thread>
#include <chrono>

#define UP    0
#define DOWN  1
#define LEFT  2
#define RIGHT 3
#define sf	  0
#define vort  1
#define F77NAME(x) x##_
extern "C" {

	void Cblacs_pinfo(int*, int*);
	void Cblacs_get(int, int, int*);
	void Cblacs_gridinit(int*, const char*, int, int);
	void Cblacs_gridinfo(int , int*, int*, int*, int*);
	void Cblacs_gridexit(int);
	void Cblacs_barrier(int, const char*);
	   
	void F77NAME(pdpbtrf)(const char& UPLO, const int& N, const int& BW,
						  double* A,  const int& JA,  int* DESCA,
						  double* AF, int& LAF, double* WORK, 
						  const int& LWORK, int* INFO);				             		 
		        		 
	 void F77NAME(pdpbtrs)(const char& UPLO, const int& N, const int& BW,
						  const int& NRHS, double* A,  const int& JA,  
						  int* DESCA, double* B,  const int& IB, 
						  int* DESCB, double* AF, int& LAF, 
						  double* WORK, const int& LWORK, int* INFO);     
	
	void F77NAME(dpbtrf) (const char& UPLO, const int& n, const int& kd,
		       	      const double* ab, const int& LDAB, int& info);
	void F77NAME(dpbtrs) (const char& UPLO, const int& n, const int& KD,
			      const int& nrhs, const double* AB, const int& ldab,
			      const double* b, const int& ldb, int& info);
}

using namespace std;

// Default Constructor
LidDrivenCavity::LidDrivenCavity()
{
}

// Destructor
LidDrivenCavity::~LidDrivenCavity()
{
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


	delete[] A;
	delete[] AF;
	delete[] b;
	delete[] b_rank;
	delete[] b_cart;
	delete[] b_scal;
	delete[] workf;
	delete[] works;
	delete[] desca;
	delete[] descb;
	delete[] ipiv;
	delete[] b_scal_nx;
	delete[] b_scal_ny;

}

// Function to store domain size Lx/Ly
void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
{ 
	Lx = xlen;
	Ly = ylen;
}

// Function to obtain store number of grid points Nx/Ny
void LidDrivenCavity::SetGridSize(int nx, int ny)
{
	Nx = nx;
	Ny = ny;

}

// Function to store partition sizes Px/Py
void LidDrivenCavity::SetPartitionSize(int ppx, int ppy)
{
	Px = ppx;
	Py = ppy;

}

// Function to store dt
void LidDrivenCavity::SetTimeStep(double deltat)
{
	dt = deltat;
}

// Function to store final time of simulation
void LidDrivenCavity::SetFinalTime(double finalt)
{
	T = finalt;
}

// Function that stores Reynolds number
void LidDrivenCavity::SetReynoldsNumber(double re)
{
	Re = re;
}

// Function that stores rank
void LidDrivenCavity::GetRank(int rr)
{
	rank = rr;
}

// Function that stores sizes
void LidDrivenCavity::GetSize(int ss)
{
	size = ss;
}

// Function that prints streamfunction matrix of a given rank
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

// Function that prints vorticity matrix of a given rank
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

// Function that initialise solver settings
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
    // Conditional statement to store deficient matrix in the last column if Px is not a divisor of Nx
	if(coords[1] == Py-1){
    	if(Nx%Py != 0){
        	loc_nx = (Nx % int( ceil(double(Nx)/double(Py)) )) + 2;
    	}
    	else{
        	loc_nx = Nx/Py + 2;
    	}
	}
    // Conditional statement to store deficient matrix in the last row if Py is not a divisor of Ny
	if(coords[0] == Px-1){
		if(Ny%Px != 0){
		    loc_ny = (Ny % int( ceil(double(Ny)/double(Px)) )) + 2;
		}
		else{
		    loc_ny = Ny/Px + 2;
		}
	}
    // Allocating size of matrices
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

	// Initialising vorticity(v), streamfunction(s) and send/recv buffers to zero
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

	// Speed of Moving Wall
	U = 1.0;
	
}

// Function that implement boundary conditions
void LidDrivenCavity::BoundaryConditions()
{

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

// Function that passes boundary points to nodes
void LidDrivenCavity::Communicate()
{	
	MPI_Request req1, req2, req3, req4, req5, req6, req7, req8,
				req9, req10, req11, req12, req13, req14, req15, req16;
	////////////////////////////////////////
    // Send/Receive to/from  DOWN Neighbours
    ////////////////////////////////////////
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
    
    //////////////////////////////////////
	// Send/Receive to/from UP Neighbours
    /////////////////////////////////////
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

    /////////////////////////////////////
	// Send/Recv to/from RIGHT Neighbours
    /////////////////////////////////////
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

    ////////////////////////////////////
	// Send/Recv to/from LEFT Neighbours
    ////////////////////////////////////
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

// Function that calculates interior points of eq (10) & (11)
void LidDrivenCavity::InteriorUpdate()
{
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
	// Calculation of interior vorticity at time t ------ (10)
	for(int i = xstart; i<loc_nx-xend; i++){
		for(int j = jstart; j<loc_ny-jend; j++){
			v[i+loc_nx*j] = -(s[i+loc_nx*j+1] - 2.0*s[i+loc_nx*j] + s[i+loc_nx*j-1])/dx/dx
					-(s[i+loc_nx*j+loc_nx] -2.0*s[i+loc_nx*j] + s[i+loc_nx*j-loc_nx])/dy/dy;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	Communicate();

	// Calculation of interior vorticity at time t+dt ---- (11)
	for(int i = xstart; i<loc_nx-xend; i++){
		for(int j = jstart; j<loc_ny-jend; j++){
			v_new[i+loc_nx*j] = v[i+loc_nx*j] 
         - dt*(s[i+loc_nx*j+loc_nx]-s[i+loc_nx*j-loc_nx])*(v[i+loc_nx*j+1]-v[i+loc_nx*j-1])/4.0/dx/dy
         + dt*(s[i+loc_nx*j+1]-s[i+loc_nx*j-1])*(v[i+loc_nx*j+loc_nx]-v[i+loc_nx*j-loc_nx])/4.0/dx/dy
         + dt*( (v[i+loc_nx*j+1] - 2.0*v[i+loc_nx*j] + v[i+j*loc_nx-1])/dx/dx 
                + (v[i+loc_nx*j+loc_nx] - 2.0*v[i+loc_nx*j] + v[i+loc_nx*j-loc_nx])/dy/dy )/Re;
		}
	}
    // Updating inner vorticity points
	for(int i = xstart; i<loc_nx-xend; i++){
		for(int j = jstart; j<loc_ny-jend; j++){
			v[i+loc_nx*j] = v_new[i+loc_nx*j];
		}
	}
}

// Function that extracts the inner nodes of vorticity(w) into vector b of linear problem Ax=b	
void LidDrivenCavity::MapRHS(){
	int x_start = 1;
	int x_end = 1;
	int y_start = 1;
	int y_end = 1;
	if(nghbrs[LEFT] == -2){
		x_start += 1;
	}
	if(nghbrs[RIGHT] == -2){
		x_end += 1;
	}
	if(nghbrs[UP] == -2){
		y_end += 1;
	}	
	if (nghbrs[DOWN] == -2){
		y_start += 1;
	}
	int c = 0;
	for(int j = y_start; j < loc_ny-y_end; j++){
		for(int i = x_start; i < loc_nx-x_end; i++){
			b_rank[c] = v[i + j*loc_nx];
		//	if (rank==size-1) cout << "Index: " << i+j*loc_nx << endl;
			c++;
		}
	}
}


void LidDrivenCavity::ScalInit(){
	// Defining variables for parallelisation
	N = (Nx-2)*(Ny-2);
	NB = ceil(double(N)/double(size));
	BW = (Nx-2);
	lda = BW+1;
	nrhs = 1;
	JA = 1;
	IB = 1;
	LA = lda*NB;
	LAF = (NB+2*BW)*BW;
	lworkf = BW*BW;
	lworks = nrhs*BW;
	

	// Arrays
	b_nx = ceil(double(Nx)/Py);
    b_ny = ceil(double(Ny)/Px);
	b_rank_L = b_nx*b_ny;
	A = new double[LA];
	AF = new double[LAF];
	ipiv = new int[NB];
	b = new double[NB];
	b_rank	= new double[b_rank_L];
	workf = new double[lworkf];
	works = new double[lworks];

	// Initialisation
	fill_n(A, LA, 0.0);
	fill_n(AF, LAF, 0.0);
	fill_n(b, NB, 0.0);
	fill_n(b_rank, b_rank_L, 0.0);
	
	// Initialising BLACS
	ncol = size;
	nrow = 1;
	Cblacs_pinfo(&rank, &size);
	Cblacs_get(0, 0, &ctx);
	Cblacs_gridinit(&ctx, "Row_Major", 1, size);
	Cblacs_gridinfo(ctx, &nrow, &ncol, &myrow, &mycol);

	// Descriptors for banded Matrix A
	desca = new int[7];
	desca[0] = 501;
	desca[1] = ctx;
	desca[2] = N;
	desca[3] = NB;
	desca[4] = 0;
	desca[5] = lda;
	desca[6] = 0;
		
	// Descriptors for RHS vector B
	descb = new int[7];
	descb[0] = 502;
	descb[1] = ctx;
	descb[2] = N;
	descb[3] = NB;
	descb[4] = 0;
	descb[5] = NB;
	descb[6] = 0;
    
	// RHS Mapping variables
	if(rank == 0){
		b_cart = new double[b_rank_L*size];
		fill_n(b_cart, b_rank_L*size, 0.0);
	}

	b_scalx = b_nx;
	b_scaly = b_ny;
	x_begin = 0;
    y_begin = 0;
	x_end = 0;
	y_end = 0;

	if (nghbrs[LEFT] == -2){
		x_begin = 1;
	}
	if (nghbrs[UP] == -2){
		y_end = 1;
	}
	if (nghbrs[DOWN] == -2){
		y_begin = 1;
		if (Ny%b_ny != 0){
			y_begin += b_ny - Ny%b_ny;
		}
	}
	if (nghbrs[RIGHT] == -2){
		x_end = 1;
		if (Nx%b_nx != 0){
			x_end += b_nx - Nx%b_nx;
		}
	}

	b_scalx -= (x_end+x_begin);
	b_scaly -= (y_end+y_begin);
	
	if (rank == 0){
		b_scal_nx = new int[size];
		b_scal_ny = new int[size];
		b_scal = new double[NB*size];
	}
	MPI_Barrier(MPI_COMM_WORLD);
	cout << "Im rank; " << rank<< "	scal_nx: " << b_scalx<< "	scal_ny: " << b_scaly << endl;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gather(&b_scalx, 1, MPI_INT, b_scal_nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gather(&b_scaly, 1, MPI_INT, b_scal_ny, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
}

void LidDrivenCavity::CholFact(){

	// Constructing banded matrix A
	double alpha = 2.0*(1.0/dx/dx + 1.0/dy/dy); // Coefficients of i,j
	double beta_x = -1.0/dx/dx;		// Coefficients of i+/-1, j
	double beta_y = -1.0/dy/dy;		// Coefficients of i/, j+/-1
	info = 0;

	for(int i = 0 ; i<NB; i++){
		if(rank*NB + i < N){
			// Filling diagonal elements
			A[i*lda+BW] = alpha;
			// Filling superdiagonals of L/R Neighbours
			if((i+rank*NB)%(Nx-2) == 0){
				A[i*lda + BW -1] = 0;
			}
			else{
				A[i*lda + BW -1] = beta_x;
			}
			// Finlling superdiagonals of U/D neighbours
			if(rank*lda*NB + lda*i >= lda*(lda-1)){
				A[i*lda] = beta_y;
			}
		}
	}


/*	// Prints banded matrix A
	if (rank == 5){
		for(int i = 0; i < lda; i++){
			for(int j = 0; j < NB; j++){
				cout << setw(8) << A[i+lda*j] << setw(8);
			}
			cout << endl;
		}
	}
*/
	MPI_Barrier(MPI_COMM_WORLD);
	F77NAME(pdpbtrf) ('U', N, BW, A, JA, desca, AF, LAF, workf, lworkf, &info);
	cout << "Im rank: " << rank << "	with fact info: " << info << endl;
}
	
	
void LidDrivenCavity::CartToScal(){
	if (rank == 0){
		int offset;
		int rank_off = 0;
		int count = 0;
		for(int c = 0; c < size; c+=Py){	
			for(int k = 0; k<b_scal_ny[c]; k++){
				offset = 0;
				for(int j = 0 ; j < Py; j++){
					for(int i = b_scal_nx[j+c]*(b_scal_ny[j+c]-1) + offset - k*b_scal_nx[j] + rank_off;
						i < b_scal_nx[j+c]*b_scal_ny[j+c] + offset - k*b_scal_nx[j] + rank_off;
						i++){
						b_scal[count] = b_cart[i];

						count++;
					}
					offset += b_rank_L;
				}
			}
			for(int i = 0; i<Py; i++){
				rank_off += b_rank_L;
			}
		}

	}
}

void LidDrivenCavity::ScalToCart(){
	if (rank == 0){
		int offset;
		int rank_off = 0;
		int count = 0;
		for(int c = 0; c < size; c+=Py){	
			for(int k = 0; k<b_scal_ny[c]; k++){
				offset = 0;
				for(int j = 0 ; j < Py; j++){
					for(int i = b_scal_nx[j+c]*(b_scal_ny[j+c]-1) + offset - k*b_scal_nx[j] + rank_off;
						i < b_scal_nx[j+c]*b_scal_ny[j+c] + offset - k*b_scal_nx[j] + rank_off;
						i++){
						b_cart[i] = b_scal[count];
						count++;
					}
					offset += b_rank_L;
				}
			}
			for(int i = 0; i<Py; i++){
				rank_off += b_rank_L;
			}
		}

	}
}

void LidDrivenCavity::MapSF(){
	int x_start = 1;
	int x_end = 1;
	int y_start = 1;
	int y_end = 1;
	if(nghbrs[LEFT] == -2){
		x_start += 1;
	}
	if(nghbrs[RIGHT] == -2){
		x_end += 1;
	}
	if(nghbrs[UP] == -2){
		y_end += 1;
	}	
	if (nghbrs[DOWN] == -2){
		y_start += 1;
	}
	counter = 0;
	for(int j = y_start; j < loc_ny-y_end; j++){
		for(int i = x_start; i < loc_nx-x_end; i++){
			s[i + j*loc_nx] = b_rank[counter];
			counter++;
		}
	}

}

// Integrating all the functions together
void LidDrivenCavity::Integrate()
{	
	// Initialising Scalapack
	ScalInit();
	// Executing cholesky factorisation
	MPI_Barrier(MPI_COMM_WORLD);
	CholFact();
	// Starting time loop
	double t_elapse = 0.0;
	while (t_elapse < T){	
		
		// Imposing Boundary Conditions
		BoundaryConditions();

		// Calculation of Interior Vorticity at time t AND t+dt
		InteriorUpdate();

		// Consolidating all of vorticity into b_rank vectors
		MapRHS();
		// Gathering all subdomain RHS into b_cart
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Gather(b_rank, b_rank_L, MPI_DOUBLE, b_cart, b_rank_L, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		/*  // Checks b_cart
			if (rank == 0) {
				for (int i = 0 ; i < size*b_rank_L; i++){
					cout << "Index : " << i <<"	b_cart: " <<  b_cart[i] << endl;
				}
			}	
		*/
		// Converting global b vector b_cart into scalapack from b_scal
		MPI_Barrier(MPI_COMM_WORLD);
		CartToScal();
		// Scattering into b for parallel solve
		MPI_Barrier(MPI_COMM_WORLD);	
		MPI_Scatter(b_scal, NB, MPI_DOUBLE, b, NB, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
		// Solving
		MPI_Barrier(MPI_COMM_WORLD);
		F77NAME(pdpbtrs) ('U', N, BW, nrhs, A, JA, desca,b, IB, descb, AF, LAF, works, lworks, &info);
		//cout<<"Im rank: " << rank << " and solve info; " << info << endl;

		// Gathering into b_scal
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Gather(b, NB, MPI_DOUBLE, b_scal, NB, MPI_DOUBLE, 0, MPI_COMM_WORLD);
/*		
		if(rank == 0){
			for(int i = 0 ; i < NB*size; i++){
				//cout<< "Index: " << i << "	b[i]" << b[i] << endl;
			}
		}
*/
		// Converting back into b_cart form
		MPI_Barrier(MPI_COMM_WORLD);
		ScalToCart();

		// Scattering into local b's
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Scatter(b_cart, b_rank_L, MPI_DOUBLE, b_rank, b_rank_L, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		// Storing into new streamfunction
		MPI_Barrier(MPI_COMM_WORLD);
		MapSF();


		MPI_Barrier(MPI_COMM_WORLD);
		Communicate();

		t_elapse += dt;

		MPI_Barrier(MPI_COMM_WORLD);
		if(rank==0) cout << "Time: " << t_elapse << endl;	
		//std::this_thread::sleep_for(std::chrono::seconds(3));
	}
	SMatPrintRank(0);	

}

void LidDrivenCavity::ExportSol(){
   

}
