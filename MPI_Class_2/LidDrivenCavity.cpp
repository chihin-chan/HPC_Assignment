#include "LidDrivenCavity.h"
#include "PoissonSolver.h"
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

using namespace std;

// Default Constructor
LidDrivenCavity::LidDrivenCavity()
{
}

// Destructor
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
	
	// Prints complete initialisation!
	if (rank == 0){
		cout << endl << "Initialisation Complete!" << endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
}

// Function that implement boundary conditions
void LidDrivenCavity::BoundaryConditions()
{
    // Speed of moving wall
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

// Function that converts solution of Ax=b into interior nodes of streamfunction(s)
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

// Integrating all the functions together
void LidDrivenCavity::Integrate(LidDrivenCavity &src)
{	
	// Try and Catch if dt is too large
	double cond = Re*Lx/(Nx-1)*Ly/(Ny-1) / 4;
	try{	
		if(dt >= cond){throw std::out_of_range("");}
	}
	catch(std::out_of_range const &e){
		if (rank==0){
		cout << "dt: " << dt << 
		", is too large and should be lesser than: " << cond << endl;
		}
	      return;   
	}


	// Calling Poisson Solver for cholesky factorising caching
	PoissonSolver Cholesky;
	Cholesky.CholFact(src);
	
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
			Cholesky.CholSolve(rhs);

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
	SMatPrintRank(0);
}

void LidDrivenCavity::ExportSol(){
   
    // Declaring a pointer of pointers to hold all values of stream function
    double* vort_add = nullptr;
    double* stream_add = nullptr;
    int nnx = loc_nx - 2;
    int nny = loc_ny - 2;
    int* nx = nullptr;
    int* ny = nullptr;

    double* v_int = new double[nnx*nny];
    double* s_int = new double[nnx*nny];

    for(int i = 0; i< nnx; i++){
        for(int j = 0; j<nny; j++){
            v_int[i+j*nnx] = v[(i+1) + loc_nx*(j+1)];
            s_int[i+j*nnx] = s[(i+1) + loc_nx*(j+1)];
        }
    }

     if(rank==0){
        nx = new int[size];
        ny = new int[size];
        vort_add = new double[Nx*Ny];
        stream_add = new double[Nx*Ny];
    }
   
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(&nnx, 1, MPI_INT, nx, 1, MPI_INT, 0, MPI_COMM_WORLD); 
    MPI_Gather(&nny, 1, MPI_INT, ny, 1, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(v_int, nnx*nny, MPI_DOUBLE, vort_add, nnx*nny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(s_int, nnx*nny, MPI_DOUBLE, stream_add, nnx*nny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD); 
    if (rank == 0){
       ofstream sOut("streamfunction.txt", ios::out | ios::trunc);
       ofstream vOut("vorticity.txt", ios::out | ios::trunc);
       sOut.precision(5);
       vOut.precision(5);
 
        int offset = 0;
        int rank_off = 0;
        
        // Moving stepping to the next rank row
        for(int c = 0; c<size; c+=Py){
             cout << "Index c: " << c << endl;
            // Printing next lower row in the rank to the last in ny
            for(int j = 0; j<ny[c]; j++){
                cout << "Index j; " << j << endl;
                cout << "j*nx[j]: " << j*nx[j] << endl;
                offset = 0;
                // Printing top row on the next rank to the last in Py
                for(int k = 0; k < Py; k++){
                    // Printing top row of rank = 0
                    for(int i = offset + nx[k+c]*(ny[k+c]-1) - j*nx[k] + rank_off;
                        i < offset + nx[k+c]*(ny[k+c]) - j*nx[k] + rank_off;
                        i++){
                        cout << "Index: i " << i <<",   SF: " << stream_add[i] << endl;
                        vOut << setw(12) << vort_add[i] << setw(12);
                        sOut << setw(12) << stream_add[i] << setw (12);
                    }
                offset += nx[k+c] * ny[k+c];
                cout << "offset: " << offset<< endl << endl;
                }
            vOut << endl;
            sOut << endl;
            }
            // Distance to offset to more to next rank row
            for(int i = 0; i<Py; i++){
                rank_off += nx[i]*ny[i];
            }
        }
        vOut.close();
        sOut.close();
    }
}
