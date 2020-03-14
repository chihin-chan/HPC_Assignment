#include <iostream>
#include <mpi.h>
#include <math.h>
#include <algorithm>

#define UP    0
#define DOWN  1
#define LEFT  2
#define RIGHT 3
#define sf	  0
#define vort  1
using namespace std;

void matPrint(double* a, int nx, int ny){
    for(int j = ny-1; j>=0;  j--){
        for(int i = 0; i<nx; i++){
            cout << a[i+j*nx] << "   ";
        }
    cout << endl;
    }
}

void vecHonPrint(double*x , int nx){
	for(int i =0; i<nx; i++){
		cout << x[i] << "	";
	}
	cout << endl;
}

void boundaryConditions(int loc_nx, int loc_ny, int loc_sub_nx, int loc_sub_ny,
						double dx, double dy, double* v, double* s,
						double* v_sub_down, double* s_sub_down, double* v_sub_right,
						double* s_sub_right, double* s_sub, double* v_sub, int* nghbrs, double U){
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
}


int main(int argc, char** argv){

    int Nx = 11;
    int Ny = 17;
    int Px = 3;
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
	int tag[2] = {0,1};
	
	MPI_Request req;
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
	

	boundaryConditions(loc_nx,loc_ny,loc_sub_nx,loc_sub_nx, dx, dy, 
						v, s, v_sub_down, s_sub_down, v_sub_right, s_sub_right, 
						s_sub, v_sub, nghbrs, U);

	///////////////////////////////////////////////
	// Packing Information and Sending To The Right
	///////////////////////////////////////////////
	if(nghbrs[RIGHT] != -2){
		double* outbuf_s;
		double* outbuf_v;
		if(nghbrs[DOWN] == -2){
			outbuf_s = new double[loc_sub_ny];
			outbuf_v = new double[loc_sub_ny];
			int k = 0;
			// Packing vorticity and streeamfunction into outbuffers
			for (int i = loc_nx - 2; i < loc_nx*loc_sub_ny; i+=loc_nx){
				outbuf_s[k] = s_sub_down[i];
				outbuf_v[k] = v_sub_down[i];
				k++;
			}
			MPI_Isend(outbuf_s, loc_sub_ny, MPI_DOUBLE, nghbrs[RIGHT], tag[sf], mygrid, &req);
			MPI_Isend(outbuf_v, loc_sub_ny, MPI_DOUBLE, nghbrs[RIGHT], tag[vort], mygrid, &req);
		}
		else{
			outbuf_s = new double[loc_ny];
			outbuf_v = new double[loc_ny];
			int k = 0;
			for (int i = loc_nx - 2; i < loc_nx*loc_ny; i+=loc_nx){
				outbuf_s[k] = s_sub_down[i];
				outbuf_v[k] = v_sub_down[i];
				k++;
			}
			MPI_Isend(outbuf_s, loc_ny, MPI_DOUBLE, nghbrs[RIGHT], tag[sf], mygrid, &req);
			MPI_Isend(outbuf_v, loc_ny, MPI_DOUBLE, nghbrs[RIGHT], tag[vort], mygrid, &req);
		}
	}

	////////////////////////////////////////////////////////////
	// Receiving Information from the left Ngbhrs and Unpacking
	////////////////////////////////////////////////////////////
	if(nghbrs[LEFT] != -2){
		double* inbuf_s;
		double* inbuf_v;
		// Recving and unpacking to bottom right deficient subdomain
		if( (nghbrs[RIGHT] == -2) && (nghbrs[DOWN] == -2) ){
			inbuf_s = new double[loc_sub_ny];
			inbuf_v = new double[loc_sub_ny];
			MPI_Irecv(inbuf_s, loc_sub_ny, MPI_DOUBLE, nghbrs[LEFT], tag[sf], mygrid, &req);
			MPI_Irecv(inbuf_v, loc_sub_ny, MPI_DOUBLE, nghbrs[LEFT], tag[vort], mygrid, &req);
			int k = 0;			
			// Unpacking inbuffers into vorticity and streamfunction
			for(int i = 0; i<loc_sub_nx*loc_sub_ny-1; i+=loc_sub_nx){
				s_sub[i] = inbuf_s[k];
				v_sub[i] = inbuf_v[k];
			}
		}
		// Receving and unpacking into right deficient subdomains
		else if((nghbrs[RIGHT] == -2) && (nghbrs[DOWN] != -2)){
			inbuf_s = new double[loc_ny];
			inbuf_v = new double[loc_ny];
			MPI_Irecv(inbuf_s, loc_ny, MPI_DOUBLE, nghbrs[LEFT], tag[sf], mygrid, &req);
			MPI_Irecv(inbuf_v, loc_ny, MPI_DOUBLE, nghbrs[LEFT], tag[vort], mygrid, &req);
			// Unpacking inbuffers into vorticity and streamfunction
			int k = 0;			
			for(int i = 0; i<loc_sub_nx*loc_ny; i+=loc_sub_nx){
				s_sub_right[i] = inbuf_s[k];
				v_sub_right[i] = inbuf_v[k];
				k++;
			}
		}
		// Receving and unpacking into bottom deficient subdomains
		else if((nghbrs[DOWN] == -2) && (nghbrs[RIGHT] != -2)){
			inbuf_s = new double[loc_sub_ny];
			inbuf_v = new double[loc_sub_ny];
			MPI_Irecv(inbuf_s, loc_sub_ny, MPI_DOUBLE, nghbrs[LEFT], tag[sf], mygrid, &req);
			MPI_Irecv(inbuf_v, loc_sub_ny, MPI_DOUBLE, nghbrs[LEFT], tag[vort], mygrid, &req);
			int k = 0;			
			// Unpacking inbuffers into vorticity and streamfunction
			for(int i = 0; i<loc_nx*loc_sub_ny-1; i+=loc_sub_nx){
				s_sub_down[i] = inbuf_s[k];
				v_sub_down[i] = inbuf_v[k];
			}
		}
		// Receving and unpacking for the rest of non-deficient domains
		else{
			inbuf_s = new double[loc_ny];
			inbuf_v = new double[loc_ny];
			MPI_Irecv(inbuf_s, loc_ny, MPI_DOUBLE, nghbrs[LEFT], tag[sf], mygrid, &req);
			MPI_Irecv(inbuf_v, loc_ny, MPI_DOUBLE, nghbrs[LEFT], tag[vort], mygrid, &req);
			// Unpacking inbuffers into vorticity and streamfunction
			int k = 0;			
			for(int i = 0; i<loc_nx*loc_ny; i+=loc_sub_nx){
				s[i] = inbuf_s[k];
				v[i] = inbuf_v[k];
				k++;
			}
		}
	}
	
	//////////////////////////////////////////////
	// Packing and sending information to the left
	//////////////////////////////////////////////
	if (nghbrs[LEFT] != -2){
		double* outbuf_s;
		double* outbuf_v;
		// Bottom right deficient subdomain
		if ((nghbrs[RIGHT] == - 2) && (nghbrs[DOWN] == -2)){
			outbuf_s = new double[loc_sub_ny];
			outbuf_v = new double[loc_sub_ny];
			int k = 0;
			for(int i = 1; i<(loc_sub_ny*loc_sub_nx) - 1; i+=loc_sub_nx){
				outbuf_s[k] = s_sub[i];
				outbuf_v[k] = v_sub[i];
				k++;
			}
			MPI_Isend(outbuf_s, loc_sub_ny, MPI_DOUBLE, nghbrs[LEFT], tag[sf], mygrid, &req);
			MPI_Isend(outbuf_v, loc_sub_ny, MPI_DOUBLE, nghbrs[LEFT], tag[vort], mygrid, &req);
		}
		// Rightmost deficient subdomains
		else if( (nghbrs[RIGHT] == -2) && (nghbrs[DOWN] != -2) ) {
			outbuf_s = new double[loc_ny];
			outbuf_v = new double[loc_ny];
			int k = 0;
			for(int i = 1; i<(loc_ny*loc_sub_nx)-1; i+=loc_sub_nx){
				outbuf_s[k] = s_sub_right[i];
				outbuf_v[k] = v_sub_right[i];
				k++;
			}
			MPI_Isend(outbuf_s, loc_ny, MPI_DOUBLE, nghbrs[LEFT], tag[sf], mygrid, &req);
			MPI_Isend(outbuf_v, loc_ny, MPI_DOUBLE, nghbrs[LEFT], tag[vort], mygrid, &req);
		}
		// Bottom deficient subdomains
		else if((nghbrs[DOWN] == -2) && (nghbrs[RIGHT] != -2)){
			outbuf_s = new double[loc_sub_ny];
			outbuf_v = new double[loc_sub_ny];
			int k = 0;
			for(int i = 1; i<(loc_sub_ny*loc_nx) - 1; i+=loc_nx){
				outbuf_s[k] = s_sub_down[i];
				outbuf_v[k] = v_sub_down[i];
			}
			MPI_Isend(outbuf_s, loc_sub_ny, MPI_DOUBLE, nghbrs[LEFT], tag[sf], mygrid, &req);
			MPI_Isend(outbuf_v, loc_sub_ny, MPI_DOUBLE, nghbrs[LEFT], tag[vort], mygrid, &req);
		}
		// Rest of the domains
		else{
			outbuf_s = new double[loc_ny];
			outbuf_v = new double[loc_ny];
			int k = 0;
			for (int i = 1; i<loc_nx*loc_ny; i+=loc_nx){
				outbuf_s[k] = s[i];
				outbuf_v[k] = v[i];
			}
			MPI_Isend(outbuf_s, loc_ny, MPI_DOUBLE, nghbrs[LEFT], tag[sf], mygrid, &req);
			MPI_Isend(outbuf_v, loc_ny, MPI_DOUBLE, nghbrs[LEFT], tag[vort], mygrid, &req);
		}
	}
	
	///////////////////////////////////////////////////////////
	// Receving information from the right nghbrs and unpacking
	///////////////////////////////////////////////////////////
	if (nghbrs[RIGHT] != -2){
		double* inbuf_s;
		double* inbuf_v;
		// Bottom deficient subdomains
		if(nghbrs[DOWN] == -2){
			inbuf_s = new double[loc_sub_ny];
			inbuf_v = new double[loc_sub_ny];
			MPI_Irecv(inbuf_s, loc_sub_ny, MPI_DOUBLE, nghbrs[RIGHT], tag[sf], mygrid, &req);
			MPI_Irecv(inbuf_v, loc_sub_ny, MPI_DOUBLE, nghbrs[RIGHT], tag[vort], mygrid, &req);		
			// Unpacking
			int k = 0;
			for (int i = loc_nx-2; i<(loc_nx*loc_sub_ny)-1; i+= loc_nx){
				s_sub_down[i] = inbuf_s[k];
				v_sub_down[i] = inbuf_v[k];
				k++;
			}
		}
		// Rest of the domains
		else{
			inbuf_s = new double[loc_ny];
			inbuf_v = new double[loc_ny];
			MPI_Irecv(inbuf_s, loc_ny, MPI_DOUBLE, nghbrs[RIGHT], tag[sf], mygrid, &req);
			MPI_Irecv(inbuf_v, loc_ny, MPI_DOUBLE, nghbrs[RIGHT], tag[vort], mygrid, &req);
			// Unpacking
			int k = 0;
			for (int i = loc_nx -2; i<(loc_nx*loc_ny)-1; i+= loc_nx){
				s[i] = inbuf_s[k];
				v[i] = inbuf_v[k];
				k++;
			}
		}
	}
	/////////////////////////////////////////////////
	// Packing and sending information to DOWN nghbrs
	/////////////////////////////////////////////////
	if (nghbrs[DOWN] != -2){
		double* outbuf_s;
		double* outbuf_v;
		// Right deficient subdomains
		if (nghbrs[RIGHT] == -2){
			outbuf_s = new double[loc_sub_nx];
			outbuf_v = new double[loc_sub_nx];
			// Packing
			int k = 0;
			for(int i = loc_sub_nx; i<2*loc_sub_nx; i++){
				outbuf_s[k] = s_sub_right[i];
				outbuf_v[k] = v_sub_right[i];
				k++;
			}
			MPI_Isend(outbuf_s, loc_sub_nx, MPI_DOUBLE, nghbrs[DOWN], tag[sf], mygrid, &req);
			MPI_Isend(outbuf_v, loc_sub_nx, MPI_DOUBLE, nghbrs[DOWN], tag[vort], mygrid, &req);
		}
		// Rest of the domains
		else{
			outbuf_s = new double[loc_nx];
			outbuf_v = new double[loc_nx];
			// Packing
			int k = 0;
			for (int i=loc_nx; i<2*loc_nx; i++){
				outbuf_s[k] = s[i];
				outbuf_v[k] = v[i];
				k++;
			}
			MPI_Isend(outbuf_s, loc_nx, MPI_DOUBLE, nghbrs[DOWN], tag[sf], mygrid, &req);
			MPI_Isend(outbuf_v, loc_nx, MPI_DOUBLE, nghbrs[DOWN], tag[vort], mygrid, &req);
		}
	}

	////////////////////////////////////////////////////
	// Receving information from UP nghbrs and unpacking
	////////////////////////////////////////////////////
	if (nghbrs[UP] != -2) {
		double* inbuf_s;
		double* inbuf_v;
		// Bottom Right deficient subdomains
		if ((nghbrs[RIGHT] == -2) && (nghbrs[DOWN] == -2)){
			inbuf_s = new double[loc_sub_nx];
			inbuf_v = new double[loc_sub_nx];
			MPI_Irecv(inbuf_s, loc_sub_nx, MPI_DOUBLE, nghbrs[UP], tag[sf], mygrid, &req);
			MPI_Irecv(inbuf_v, loc_sub_nx, MPI_DOUBLE, nghbrs[UP], tag[vort], mygrid, &req);
			int k = 0;		
			for(int i = (loc_sub_nx)*(loc_sub_ny-1); i<(loc_sub_nx*loc_sub_ny)-1; i++){
				s_sub[i] = inbuf_s[k];
				v_sub[i] = inbuf_v[k];
				k++;
			}
		}
		// Bottom Row Deficient subdomains
		else if( (nghbrs[RIGHT] != -2) && (nghbrs[DOWN] == -2)){
			inbuf_s = new double[loc_nx];
			inbuf_v = new double[loc_nx];
			MPI_Irecv(inbuf_s, loc_nx, MPI_DOUBLE, nghbrs[UP], tag[sf], mygrid, &req);
			MPI_Irecv(inbuf_v, loc_nx, MPI_DOUBLE, nghbrs[UP], tag[vort], mygrid, &req);
			int k = 0;			
				for(int i=(loc_nx)*(loc_sub_ny-1); i<(loc_nx)*(loc_sub_ny)-1; i++){
				s_sub_down[i] = inbuf_s[k];
				v_sub_down[i] = inbuf_v[k];
				k++;
			}
		}

		else if ( (nghbrs[RIGHT] == -2) && (nghbrs[DOWN] != -2)){
			inbuf_s = new double[loc_sub_nx];
			inbuf_v = new double[loc_sub_nx];
			MPI_Irecv(inbuf_s, loc_sub_nx, MPI_DOUBLE, nghbrs[UP], tag[sf], mygrid, &req);
			MPI_Irecv(inbuf_v, loc_sub_nx, MPI_DOUBLE, nghbrs[UP], tag[vort], mygrid, &req);
			int k = 0;
			for(int i = (loc_sub_nx)*(loc_ny-1); i<(loc_sub_nx*loc_ny)-1; i++){
				s_sub_right[i] = inbuf_s[k];
				v_sub_right[i] = inbuf_v[k];
				k++;
			}
		}

		else{
			inbuf_s = new double[loc_nx];
			inbuf_v = new double[loc_nx];
			MPI_Irecv(inbuf_s, loc_nx, MPI_DOUBLE, nghbrs[UP], tag[sf], mygrid, &req);
			MPI_Irecv(inbuf_v, loc_nx, MPI_DOUBLE, nghbrs[UP], tag[vort], mygrid, &req);	
			int k = 0;
			for(int i = (loc_nx)*(loc_ny-1); i<(loc_nx*loc_ny)-1; i++){
				s[i] = inbuf_s[k];
				v[i] = inbuf_v[k];
				k++;
			}
		}
	}

	///////////////////////////////////////////////
	// Packing and sending information to UP nghbrs
	///////////////////////////////////////////////
	if (nghbrs[UP] != -2){
		double* outbuf_s;
		double* outbuf_v;
		// Bottom Right Deficient Subdomain
		if((nghbrs[DOWN] == -2) && (nghbrs[RIGHT] == -2)){
			outbuf_s = new double[loc_sub_nx];
			outbuf_v = new double[loc_sub_nx];
			int k = 0;
			// Packing
			for (int i = (loc_sub_nx)*(loc_sub_ny-1); i<(loc_sub_nx*loc_sub_ny)-1; i++){
				outbuf_s[k] = s_sub[k];
				outbuf_v[k] = v_sub[k];
				k++;
			}
		cout << "Im rank, " << rank << " sending to nghbr of rank: " << nghbrs[UP] << endl;
			MPI_Isend(outbuf_s, loc_sub_nx, MPI_DOUBLE, nghbrs[UP], tag[sf], mygrid, &req);
			MPI_Isend(outbuf_v, loc_sub_nx, MPI_DOUBLE, nghbrs[UP], tag[vort], mygrid, &req);
		}
		// Bottom deficient subdomains
		else if((nghbrs[DOWN] == -2) && (nghbrs[RIGHT] != -2)){
			outbuf_s = new double[loc_nx];
			outbuf_v = new double[loc_nx];
			int k = 0;
			// Packing
			for (int i = (loc_nx)*(loc_sub_ny-1); i<loc_nx*loc_sub_ny-1; i++){
				outbuf_s[k] = s_sub_down[i];
				outbuf_v[k] = v_sub_down[i];
				k++;
			}
		cout << "Im rank, " << rank << " sending to nghbr of rank: " << nghbrs[UP] << endl;
			MPI_Isend(outbuf_s, loc_nx, MPI_DOUBLE, nghbrs[UP], tag[sf], mygrid, &req);
			MPI_Isend(outbuf_v, loc_nx, MPI_DOUBLE, nghbrs[UP], tag[vort], mygrid, &req);
		}

		// Right deficient subdomains
		else if( (nghbrs[DOWN] != -2) && (nghbrs[RIGHT] == -2)){
			outbuf_s = new double[loc_sub_nx];
			outbuf_v = new double[loc_sub_nx];
			int k = 0;
			// Packing
		cout << "Im rank, " << rank << " sending to nghbr of rank: " << nghbrs[UP] << endl;
			for (int i = (loc_sub_nx)*(loc_ny-1); i<(loc_sub_nx*loc_ny)-1;i++){
				outbuf_s[k] = s_sub_right[i];
				outbuf_v[k] = v_sub_right[i];
				k++;
			}
			MPI_Isend(outbuf_s, loc_sub_nx, MPI_DOUBLE, nghbrs[UP], tag[sf], mygrid, &req);
			MPI_Isend(outbuf_v, loc_sub_nx, MPI_DOUBLE, nghbrs[UP], tag[vort], mygrid, &req);
		}

		else{
			outbuf_s = new double[loc_nx];
			outbuf_v = new double[loc_nx];
			int k = 0;
			// Packing
		cout << "Im rank, " << rank << " sending to nghbr of rank: " << nghbrs[UP] << endl;
			for (int i = loc_nx*(loc_ny-1); i<(loc_nx*loc_ny)-1; i++){
				outbuf_s[k] = s[i];
				outbuf_v[k] = v[i];
				k++;
			}
			MPI_Isend(outbuf_s, loc_nx, MPI_DOUBLE, nghbrs[UP], tag[sf], mygrid, &req);
			MPI_Isend(outbuf_v, loc_nx, MPI_DOUBLE, nghbrs[UP], tag[vort], mygrid, &req);
		}
	}

	///////////////////////////////////////////
	// Receiving from DOWN nghbrs and unpacking
	///////////////////////////////////////////
	if (nghbrs[DOWN] != -2){
		double* inbuf_s;
		double* inbuf_v;
		if (nghbrs[RIGHT] == -2){
			inbuf_s = new double[loc_sub_nx];
			inbuf_v = new double[loc_sub_nx];
			cout <<"im Rank: "<<rank<<" and recv from rank: " << nghbrs[DOWN] << endl;			
			MPI_Irecv(inbuf_s, loc_sub_nx, MPI_DOUBLE, nghbrs[DOWN], tag[sf], mygrid, &req);
			MPI_Irecv(inbuf_v, loc_sub_nx, MPI_DOUBLE, nghbrs[DOWN], tag[vort], mygrid, &req);
			int k = 0 ;
			for (int i = 0; i<loc_sub_nx; i++){
				s_sub_right[i] = inbuf_s[k];
				v_sub_right[i] = inbuf_v[k];
				k++;
			}
		}
		else{
			inbuf_s = new double[loc_nx];
			inbuf_v = new double[loc_nx];
			MPI_Irecv(inbuf_s, loc_nx, MPI_DOUBLE, nghbrs[DOWN], tag[sf], mygrid, &req);
			MPI_Irecv(inbuf_v, loc_nx, MPI_DOUBLE, nghbrs[DOWN], tag[sf], mygrid, &req);
			int k = 0 ;
			cout <<"im Rank: "<<rank<<" and recv from rank: " << nghbrs[DOWN] << endl;			
			for (int i = 0; i< loc_nx; i++){
				s[i] = inbuf_s[k];
				v[i] = inbuf_v[k];
				k++;
			}
		}
	}
			
		
			
		
	/*printf("rank = %2d coords = %2d%2d neighbors(u,d,l,r) = %2d %2d %2d %2d\n\n",
        rank,coords[0],coords[1],nghbrs[UP],nghbrs[DOWN],nghbrs[LEFT],nghbrs[RIGHT]); */
  
MPI_Finalize();
}

