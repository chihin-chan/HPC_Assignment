#include <iostream>
#include <mpi.h>
#include <math.h>
#include <algorithm>

#define UP    0
#define DOWN  1
#define LEFT  2
#define RIGHT 3
using namespace std;

void matPrint(double* a, int nx, int ny){
    for(int j = ny-1; j>=0;  j--){
        for(int i = 0; i<nx; i++){
            cout << a[i+j*nx] << "   ";
        }
    cout << endl;
    }
}


int main(int argc, char** argv){

    int Nx = 11;
    int Ny = 11;
    int Px = 2;
    int Py = 3;
    int Lx = 1;
    int Ly = 1;
    double dx = double(Lx)/(Nx-1.0);
    double dy = double(Ly)/(Ny-1.0);
    double U = 1.0;
    int size;
    int rank;
	int nghbrs[4];
    const int ndims = 2;
    const int dims[ndims] = {Px, Py};
    const int periods[ndims] = {0, 0};
	int coords[ndims];
	
	MPI_Comm mygrid;

    // Initialising MPI
    MPI_Init(&argc, &argv);
    
    // Retreiving size
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Generating a communicator for 2-D Cartesian Topology
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, 0, &mygrid);
	
	// Retreving Rank
    MPI_Comm_rank(mygrid, &rank); 

	// Retreving coordinates of Rank
	MPI_Cart_coords(mygrid, rank, 2, coords);

    // Obtain Shifted Source and Destination ranks in both directions
    MPI_Cart_shift(mygrid, 0, 1, &nghbrs[UP], &nghbrs[DOWN]);
    MPI_Cart_shift(mygrid, 1, 1, &nghbrs[LEFT], &nghbrs[RIGHT]);

    // Defining vorticity and streafunction matrix for subdomains
    int loc_nx = ceil(double(Nx)/double(Py)) + 2;
    int loc_ny = ceil(double(Ny)/double(Px)) + 2;
    int loc_sub_nx;
    int loc_sub_ny;
    if(Nx%Py != 0){    
        loc_sub_nx = (Nx % int( ceil(double(Nx)/double(Py)) )) + 2;
    }
    else{
        loc_sub_nx = Nx/Py + 2;
    }
    if(Ny%Py != 0){
        loc_sub_ny = (Ny % int( ceil(double(Ny)/double(Px)) )) + 2;
    }
    else{
        loc_sub_nx = Ny/Px + 2;
    }

    double* s = new double[loc_nx*loc_ny];
    double* v = new double[loc_nx*loc_ny];
	fill_n(s, loc_nx*loc_ny, 0.0);
	fill_n(v, loc_nx*loc_ny, 1.0);
  
	// For x Deficient Right Subdomains except bottom right
    double* s_sub_right = new double[loc_sub_nx*loc_ny];
    double* v_sub_right = new double[loc_sub_nx*loc_ny];
	fill_n(s_sub_right, loc_sub_nx*loc_ny, 0.0);
	fill_n(v_sub_right, loc_sub_nx*loc_ny, 1.0);
   
	// For y Deficient Bottom Subdomains except bottom right
    double* s_sub_down = new double[loc_nx*loc_sub_ny];
    double* v_sub_down = new double[loc_nx*loc_sub_ny];
	fill_n(s_sub_down, loc_nx*loc_sub_ny, 0.0);
	fill_n(v_sub_down, loc_nx*loc_sub_ny, 1.0);

    // For x and y deficient bottom right subdomain
    double* s_sub = new double[loc_sub_nx*loc_sub_ny];
    double* v_sub = new double[loc_sub_nx*loc_sub_ny];
	fill_n(s_sub, loc_sub_nx*loc_sub_ny, 0.0);
	fill_n(v_sub, loc_sub_nx*loc_sub_ny, 1.0);    
	
    ///////////////
    // TOP DOMAINS
    ///////////////
    if(nghbrs[UP] == -2){
        // Top Right Subdomain which is deficient in x
        if(nghbrs[RIGHT] == -2){
            for(int i = loc_sub_nx*(loc_ny-1)-loc_sub_nx+1; i<(loc_sub_nx*loc_ny)-1-loc_sub_nx;i++){
                v_sub_right[i] = (s_sub_right[i] - s_sub_right[i-loc_sub_nx])*2.0/dy/dy - 2.0*U/dy;
            }

        }
        // Rest of Top subdomains
        else{
            for(int i = loc_nx*(loc_ny-1)-loc_nx+1; i<(loc_nx*loc_ny)-1-loc_nx; i++){
                v[i] = (s[i] - s[i-loc_nx])*2.0/dy/dy - 2.0*U/dy;
            }

        }
    }

    ////////////////
    // LEFT DOMAINS
    ////////////////
    if(nghbrs[LEFT] == -2){
        // Bottom Left Subdomain which is deficient in y
        if(nghbrs[DOWN] == -2){
            // Left B.C for vorticity
            for(int i = 1+loc_nx; i<loc_nx*(loc_sub_ny-1); i+=loc_nx){
                v_sub_down[i] = (s_sub_down[i] - s_sub_down[i+1])*2.0/dx/dx;
            }
        }
        // Rest of Right subdomains
        else{
            for(int i=1+loc_nx; i<loc_nx*(loc_ny-1); i+=loc_nx){
                v[i] = (s[i] - s[i+1])*2.0/dx/dx;
            }
        }
    }
  	/////////////////////////////////////////
    // RIGHT DOMAINS (Potentially Deficient)
    /////////////////////////////////////////
    if(nghbrs[RIGHT] == -2){
        // Bottom right subdomain which is deficient in x and y
        if (nghbrs[DOWN] == -2){
            // Bottom Right Subdomain
            for(int i = (2*loc_sub_nx) - 2; i < (loc_sub_nx*loc_sub_ny)-2; i+=loc_sub_nx){
                v_sub[i] = (s_sub[i] - s_sub[i-1])*2.0/dx/dx;  
            }
        }
        // Rest of the Right subdomains
        else{  
            for(int i = (2*loc_sub_nx) - 2; i < (loc_sub_nx*loc_ny) - 2; i += loc_sub_nx){
                v_sub_right[i] = (s_sub_right[i] - s_sub_right[i-1])*2.0/dx/dx;
            }
        }
    }
 
    
    ////////////////////////////////////////
    // BOTTOM DOMAINS (Potentially Deficient)
    ///////////////////////////////////////
    if(nghbrs[DOWN] == -2){
        // Bottom right subdomain which is deficient in x and y
        if (nghbrs[RIGHT] == -2){
            for(int i = loc_sub_nx + 1; i < (2*loc_sub_nx) - 1; i++){
                v_sub[i] = (s_sub[i] - s_sub[i+loc_sub_nx])*2/dy/dy;
            }
        }
        // Rest of the bottom subdomains 
        else{
            for(int i = loc_nx + 1; i < (2*loc_nx) - 1; i++){
                v_sub_down[i] = (s_sub_down[i] - s_sub_down[i+loc_nx])*2/dy/dy;
            }
        }
    }
	

	printf("rank = %2d coords = %2d%2d neighbors(u,d,l,r) = %2d %2d %2d %2d\n\n",
        rank,coords[0],coords[1],nghbrs[UP],nghbrs[DOWN],nghbrs[LEFT],nghbrs[RIGHT]);
  
MPI_Finalize();
}

