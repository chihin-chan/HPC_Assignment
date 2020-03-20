#include "LidDrivenCavity.h"
#include "PoissonSolver.h"
#include <math.h> 
#include <iostream>
#include <iomanip>
#include <mpi.h>

using namespace std;

// CONSTRUCTOR
LidDrivenCavity::LidDrivenCavity(){
}

// DESTRUCTOR
LidDrivenCavity::~LidDrivenCavity(){
  //cout << "End LidDrivenCavity" << endl;
}

// METHODS FOR THE CLASS LidDrivenCavity

	// setting domain size
void LidDrivenCavity::SetDomainSize(double xlen, double ylen){
	Lx = xlen;
	Ly = ylen;
}

	// setting grid size
void LidDrivenCavity::SetGridSize(int nx, int ny){
	Nx = nx;
	Ny = ny;
}

	// setting partitions
void LidDrivenCavity::SetPartitions(int px, int py){
	Px = px;
	Py = py;
}
        
	// setting timestep
void LidDrivenCavity::SetTimeStep(double deltat){
	dt = deltat;
}
	// setting final time
void LidDrivenCavity::SetFinalTime(double finalt){
	T = finalt;
}
	// setting reynolds number
void LidDrivenCavity::SetReynoldsNumber(double re){
	Re = re;
}

void LidDrivenCavity::Initialise(int world_rank, int world_size){
	

	// Dimentions of the partitions in x- and y-dir
	dims[0] = Py;
	dims[1] = Px;
	// Periodic boundary conditions not activated
	period[0] = 0;
	period[1] = 0;
	  	
	if (world_size == Px*Py){
	
		// Creating vectors vorticity and stream function for each processor
		nx = ceil((double)Nx/Px);	// Points in x-dir for processor
		ny = ceil((double)Ny/Py);	// Points in y-dir for processor
		Npx = nx+2;	// Points + neighbours x-dir
		Npy = ny+2;	// Points + neighbours y-dir
		Np = (Npx)*(Npy); // Total n of points in each processor
		
		// Values needed for the algorithm
		deltax = Lx/(Nx-1);
		deltay = Ly/(Ny-1);
		dlxy   = 4*deltax*deltay;
		pdltx  = pow(deltax,2);
		pdlty  = pow(deltay,2);		
		
		// Vorticity, vorticity at tn+1 and streamfunction
		v = new double [Np];
		vt= new double [Np];		
		s = new double [Np];
		
		// Initialising vorticity and stream function
		for (int i = 0; i<Np; i++){
			v[i] = 0;
			s[i] = 0;
			vt[i] = 0;
		}
	
	   	///////////////////////////////////////////////////////////////////////   TEST
	/*
		for (int i = 0; i<Np; i++){
			v[i] = i;
			s[i] = i;
		}
		for (int i = 0; i<Np; i++){
			v[i] = myrank;
			s[i] = myrank;
		}
	*/
    	///////////////////////////////////////////////////////////////////////   TEST	
		
		// vectors for send and recieve
		svsendu = new double [2*nx];
		svrecvu = new double [2*nx];
		svsendd = new double [2*nx];
		svrecvd = new double [2*nx];
		svsendl = new double [2*ny];
		svrecvl = new double [2*ny];
		svsendr = new double [2*ny];
		svrecvr = new double [2*ny];
		
		// Initialising
		for (int i = 0; i<2*nx; i++){
			svsendu[i] = 0.0;
			svrecvu[i] = 0.0;
			svsendd[i] = 0.0;
			svrecvd[i] = 0.0;
		}
		for (int i = 0; i<2*ny; i++){		
			svsendl[i] = 0.0;
			svrecvl[i] = 0.0;	
			svsendr[i] = 0.0;
			svrecvr[i] = 0.0;
		}
		
		// Make a new communicator 2-D Cartesian topology
		MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, period, reorder, &cart_comm);
		
		// Get my rank in the new communicator
		MPI_Comm_rank(cart_comm, &myrank);
		
		// Det. process coords in cartesian topology of my rank
    	MPI_Cart_coords(cart_comm, myrank, ndims, coord);
    	
    	// Obtain the shift source and destination ranks in both directions
		MPI_Cart_shift(cart_comm, 0, 1, &nbrs[DOWN], &nbrs[UP]);	// y-dir
		MPI_Cart_shift(cart_comm, 1, 1, &nbrs[LEFT], &nbrs[RIGHT]);	// x-dir
		
		
	
	} else {
		cout << " Please, insert a number of processors compatible with the partitions in x and y directions (Px and Py). " 
		<< endl;
	}
}

void LidDrivenCavity::Communication(){


 	// MPI variables needed for send, recieve and wait communications 
 	MPI_Request req1, req2, req3, req4, req5,req6, req7, req8;
		
	// Not UP
	if (nbrs[UP] != -2){  							
					
		// Preparing the vector to send
		for (int i = 0; i<nx; i++){
        	svsendu[i]    = v[fpUP+i]; 
        	svsendu[i+nx] = s[fpUP+i];
		}
		// SEND to UP 							   		
		MPI_Isend(svsendu, 2*nx, MPI_DOUBLE, nbrs[UP], tag, cart_comm, &req1);
	
		// RECIEVE from UP
        MPI_Irecv(svrecvu, 2*nx, MPI_DOUBLE, nbrs[UP], tag, cart_comm, &req2);
    	
    	// Wait until send and recieve are completed		
        MPI_Wait(&req1, MPI_STATUS_IGNORE);
        MPI_Wait(&req2, MPI_STATUS_IGNORE);
         
        // Saving the vector recieved  	        	          		
        for (int i = 0; i<nx; i++){
        	v[fpUP+Npx+i] = svrecvu[i];    
        	s[fpUP+Npx+i] = svrecvu[i+nx];
        }
	}
	
	// Not DOWN	
	if (nbrs[DOWN] != -2){  							
		
 		// Preparing the vector to send
		for (int i = 0; i<nx; i++){
        	svsendd[i]    = v[fpDOWN+i];
        	svsendd[i+nx] = s[fpDOWN+i];
		}
		// SEND to DOWN			
		MPI_Isend(svsendd, 2*nx, MPI_DOUBLE, nbrs[DOWN], tag, cart_comm,  &req3);
		
		// RECIEVE from DOWN	
    	MPI_Irecv(svrecvd, 2*nx, MPI_DOUBLE, nbrs[DOWN], tag, cart_comm,  &req4);
    	
    	// Wait until send and recieve are completed		
        MPI_Wait(&req3, MPI_STATUS_IGNORE);
        MPI_Wait(&req4, MPI_STATUS_IGNORE);
        
       // Saving the vector recieved
        for (int i = 0; i<nx; i++){
        	v[fpDOWN-Npx+i] = svrecvd[i];
        	s[fpDOWN-Npx+i] = svrecvd[i+nx];
		}
    }
    
	// Not LEFT
	if (nbrs[LEFT] != -2){  							
		
		// Preparing the vector to send
		for (int i = 0; i<ny; i++){
        	svsendl[i]    = v[fpLEFT+i*Npx];   
        	svsendl[i+ny] = s[fpLEFT+i*Npx];
		}
		// SEND to LEFT    				
    	MPI_Isend(svsendl, 2*ny, MPI_DOUBLE, nbrs[LEFT], tag, cart_comm,  &req5);
    	
		// RECIEVE from LEFT	
    	MPI_Irecv(svrecvl, 2*ny, MPI_DOUBLE, nbrs[LEFT], tag, cart_comm,  &req6);
    	
    	// Wait until send and recieve are completed		
        MPI_Wait(&req5, MPI_STATUS_IGNORE);
        MPI_Wait(&req6, MPI_STATUS_IGNORE);
        
        // Saving the vector recieved
        for (int i = 0; i<ny; i++){
        	v[fpLEFT+i*Npx-1] = svrecvl[i];    
        	s[fpLEFT+i*Npx-1] = svrecvl[i+ny];
		}
    }
    
	// Not RIGHT   	
	if (nbrs[RIGHT] != -2){  							
					
		// Preparing the vector to send
		for (int i = 0; i<ny; i++){
        	svsendr[i]    = v[fpRIGHT+i*Npx];  
        	svsendr[i+ny] = s[fpRIGHT+i*Npx];
		}
		// SEND to RIGHT    	
    	MPI_Isend(svsendr, 2*ny, MPI_DOUBLE, nbrs[RIGHT], tag, cart_comm,  &req7);
    	
		// RECIEVE from RIGHT	
    	MPI_Irecv(svrecvr, 2*ny, MPI_DOUBLE, nbrs[RIGHT], tag, cart_comm,  &req8);

    	// Wait until send and recieve are completed		
        MPI_Wait(&req7, MPI_STATUS_IGNORE);
        MPI_Wait(&req8, MPI_STATUS_IGNORE);
        
        // Saving the vector recieved
    	for (int i = 0; i<ny; i++){
        	v[fpRIGHT+i*Npx+1] = svrecvr[i];   
        	s[fpRIGHT+i*Npx+1] = svrecvr[i+ny]; 
		}
   	} 

}


void LidDrivenCavity::Integrate(LidDrivenCavity &LDC){

	rankdisplay = 7;
	
	// Variables for the algorithm
	int    k;
	double a;
	double b;
	double c;
	double d;
	int Nextray;
	int Nextrax;
	
	// Creating a variable type PoissonSolver
	PoissonSolver SF;
		
	// Variables for changes in processor points	
	firstpointx = 1;
	firstpointy = 1; 
	lastpointx = Npx-1;
	lastpointy = Npy-1; 
	
	
	// Uneven grid points (extra points in last processor x and y)
    if (nbrs[UP] == -2) {
    	Nextray = ny*Py-Ny;
    	for (int i = 0; i< Nextray; i++) lastpointy = lastpointy - 1;
    }
    if (nbrs[RIGHT] == -2) {
    	Nextrax = nx*Px-Nx;
    	for (int i = 0; i< Nextrax; i++) lastpointx = lastpointx - 1;
    }
    
    
    // Boundaries    
	// UP first and last points
	fpUP = firstpointx + Npx*(lastpointy-1); 
	lpUP = lastpointx+(lastpointy-1)*Npx;;
	// DOWN 				
	fpDOWN = firstpointx + Npx; 
	lpDOWN = lastpointx + Npx;  
	// LEFT
	fpLEFT = firstpointx + Npx; 
	lpLEFT = firstpointx + (lastpointy-1)*Npx;
	// RIGHT
	fpRIGHT = lastpointx-1 + Npx; 
	lpRIGHT = lastpointx + (lastpointy-1)*Npx;
	
	// Changes because of boundaries
    if (nbrs[UP] == -2) lastpointy = lastpointy-1; 		 // One less last row in the loop
    if (nbrs[DOWN] == -2) firstpointy = firstpointy + 1; // Second row is first row in the loop
	if (nbrs[LEFT] == -2) firstpointx = firstpointx + 1; // Second point is first in the loop 
    if (nbrs[RIGHT] == -2) lastpointx = lastpointx - 1;  // One less last point in the loop 
    
    // Initializing the Poisson solver, constructing and factorizing the A matrix
    SF.Initializing(LDC);
	SF.MatrixConstr(LDC);
		
	
	// Time loop
    for (double dlt = dt; dlt <= T; dlt = dlt+dt ){    
    	
    	if (myrank==rankdisplay) cout<< endl << "--------------------------- Time step " << dlt;
    			
		// _______________________________________________  BOUNDARY CONDITIONS
		
		// ******************************************************************* UP  
			
    	if (nbrs[UP] == -2){    	
			// BCs UP    	
			for (int i = fpUP; i < lpUP; i++){ 
   				v[i] = (s[i]-s[i-Npx])*2/pdlty - 2*U/deltay;
			}
		}	
		
		// ******************************************************************* DOWN      		
    		    		
    	if (nbrs[DOWN] == -2){    
			// BCs DOWN
			for (int i = fpDOWN; i < lpDOWN; i++){
				v[i] = (s[i]-s[i+Npx])*2/pdlty;
			}
		}	  

		// ******************************************************************* LEFT  
    	
    	if (nbrs[LEFT] == -2){   
			// BCs LEFT
			for (int i = fpLEFT; i < lpLEFT+1; i=i+Npx){
				v[i] = (s[i]-s[i+1])*2/pdltx;
			}	
		}		
    	
		// ******************************************************************* RIGHT   
		    		
    	if (nbrs[RIGHT] == -2){   
			// BCs RIGHT
			for (int i = fpRIGHT; i < lpRIGHT; i=i+Npx){
				v[i] = (s[i]-s[i-1])*2/pdltx;
			}
		}	 
		
		// ______________________________________________________________________ ALGORITHM
						
		
		// ****************************** Computing vorticity in the interior at time t
		
		for (int j = firstpointy; j < lastpointy; j++){
			for (int i = firstpointx; i < lastpointx; i++){
				k = Npx*j+i;
				v[k] = -(s[k+1]-2*s[k]+s[k-1])/pdltx-(s[k+Npx]-2*s[k]+s[k-Npx])/pdlty;
			}
		}
			
		// PRINTING MATRIX AFTER ALGORITHM I --- vorticity t+dlt
    	if (myrank == rankdisplay){	
    		//cout<< endl << "____________________________________ Time step " << dlt;	
			cout << endl << "Processor " << myrank << endl
			<< " MATRIX v after ALGORITHM I with firstpointx and y : " << firstpointx << ", " << firstpointy << " and last x and y: " << lastpointx << ", " << lastpointy << endl << endl;
			for (int j = Npy-1; j > -1; j--){
				for (int i = 0; i<Npx; i++){
					cout <<  setw(12) << v[Npx*j+i] << "	";
				} cout << endl;
			} cout << endl;
		}
			
       	// ******************************  Wait all + communication between processors
       		
   		// Blocks until all processes have called it
   		MPI_Barrier(cart_comm);	
   		
		// Send and recieve communication
		Communication();
		       		
   		// Blocks until all processes have called it
   		MPI_Barrier(cart_comm);	
				
		// ****************************** computing the interior vorticity at time t+dlt
		
		for (int j = firstpointy; j < lastpointy; j++){
			for (int i = firstpointx; i < lastpointx; i++){
    			k = Npx*j+i; 
    			a = (s[k+Npx]-s[k-Npx])*(v[k+1]-v[k-1])/dlxy;
    			b = (s[k+1]-s[k-1])*(v[k+Npx]-v[k-Npx])/dlxy;
    			c = (v[k+1]-2*v[k]+v[k-1])/pdltx;
    			d = (v[k+Npx]-2*v[k]+v[k-Npx])/pdlty;
    			vt[k] = dt*(b - a + (c + d)/Re) + v[k];
    		}
    	}
    	  			
    			
    	// ************************************* save vorticity at time t+dlt vt[] in v[]
    	
    	for (int j = firstpointy; j < lastpointy; j++){
			for (int i = firstpointx; i < lastpointx; i++){
    			k = Npx*j+i; 
    			v[k] = vt[k];
    		}
    	}
    	
    	
    	// PRINTING MATRIX AFTER ALGORITHM II --- vorticity t+dlt
    	if (myrank == rankdisplay){	
    		//cout<< endl << "____________________________________ Time step " << dlt;	
			cout << endl << "Processor " << myrank << endl
			<< " MATRIX v after ALGORITHM II: " << endl << endl;
			for (int j = Npy-1; j > -1; j--){
				for (int i = 0; i<Npx; i++){
					cout <<  setw(12) << v[Npx*j+i] << "	";
				} cout << endl;
			} cout << endl;
		}
		
		    	
    	// ******************************************************* calling POISSON SOLVER
   		
    	// calling the function that solves Poisson in a separate class
    	SF.Solver(LDC);
    	
    	// ******************************  Wait all and communicate
       		
   		// Blocks until all processes have called it
   		MPI_Barrier(cart_comm);	
   		
		// Send and recieve communication
		Communication();
		       		
   		// Blocks until all processes have called it
   		MPI_Barrier(cart_comm);	
   		
   		
   		if (myrank == rankdisplay){	
			cout << endl << "Processor " << myrank << endl
			<< " MATRIX SOLUTION s: " << endl << endl;
			for (int j = Npy-1; j > -1; j--){
				for (int i = 0; i<Npx; i++){
					cout <<  setw(12) << s[Npx*j+i] << "	";
				} cout << endl;
			} cout << endl;
		}
		
    	
    	
    	
    } // end of time loop
    
     
    // Free dynamic memory --------------------------------------------------------------
    delete[] v;
    delete[] s;
    
    delete[] svsendu;
    delete[] svsendd;
    delete[] svsendl;
    delete[] svsendr;
    delete[] svrecvu;
    delete[] svrecvd;
    delete[] svrecvl;
    delete[] svrecvr;    

} // End of integrate

     
    
    ///////////////////////////////////////////////////////////////////////	  PRINT MATRIX
    
    /*
    	// PRINTING MATRIX
    	if (myrank == 7){	
    		//cout<< endl << "____________________________________ Time step " << dlt;	
			cout << endl << "Processor " << myrank << endl
			<< " MATRIX v: " << endl << endl;
			for (int j = Npy-1; j > -1; j--){
				for (int i = 0; i<Npx; i++){
					cout <<  setw(12) << vt[Npx*j+i] << "	";
				} cout << endl;
			} cout << endl;
		}
		
		
    	// PRINTING MATRIX BEGGINGING TIME STEP
    	if (myrank == rankdisplay){	
    		//cout<< endl << "____________________________________ Time step " << dlt;	
			cout << endl << "Processor " << myrank << endl
			<< " MATRIX v after Poisson: " << endl << endl;
			for (int j = Npy-1; j > -1; j--){
				for (int i = 0; i<Npx; i++){
					cout <<  setw(12) << v[Npx*j+i] << "	";
				} cout << endl;
			} cout << endl;
		}
		
		// PRINTING MATRIX AFTER BCS
    	if (myrank == rankdisplay){	
    		//cout<< endl << "____________________________________ Time step " << dlt;	
			cout << endl << "Processor " << myrank << endl
			<< " MATRIX v after BCs: " << endl << endl;
			for (int j = Npy-1; j > -1; j--){
				for (int i = 0; i<Npx; i++){
					cout <<  setw(12) << v[Npx*j+i] << "	";
				} cout << endl;
			} cout << endl;
		}		
		
		// PRINTING MATRIX AFTER ALGORITHM 1
    	if (myrank == rankdisplay){	
    		//cout<< endl << "____________________________________ Time step " << dlt;	
			cout << endl << "Processor " << myrank << endl
			<< " MATRIX v after ALGORITHM I: " << endl << endl;
			for (int j = Npy-1; j > -1; j--){
				for (int i = 0; i<Npx; i++){
					cout <<  setw(12) << v[Npx*j+i] << "	";
				} cout << endl;
			} cout << endl;
		}
		
		// PRINTING MATRIX AFTER SEND AND RECIEVE
    	if (myrank == rankdisplay){	
    		//cout<< endl << "____________________________________ Time step " << dlt;	
			cout << endl << "Processor " << myrank << endl
			<< " MATRIX v after SEND AND RECIEVE: " << endl << endl;
			for (int j = Npy-1; j > -1; j--){
				for (int i = 0; i<Npx; i++){
					cout <<  setw(12) << v[Npx*j+i] << "	";
				} cout << endl;
			} cout << endl;
		}
		
		
		
		*/
		
    ///////////////////////////////////////////////////////////////////////	  PRINT MATRIX
    
    
    
         	
	///////////////////////////////////////////////////////////////////////   TEST
	// Checking send and recieve    	
	/*if (myrank == Re){		
			cout << endl << endl << " MATRIX v: " << endl << endl;
			for (int j = Npy-1; j > -1; j--){
				for (int i = 0; i<Npx; i++){
					//cout << "v[" << Npx*j+i << "]" << "	";
					cout << v[Npx*j+i] << "	";
				} cout << endl;
			}
			cout << endl << endl << " MATRIX s: " << endl << endl;
			for (int j = Npy-1; j > -1; j--){
				for (int i = 0; i<Npx; i++){
					//cout << "v[" << Npx*j+i << "]" << "	";
					cout << s[Npx*j+i] << "	";
				} cout << endl;
			}
    	}
    */
	///////////////////////////////////////////////////////////////////////   TEST   

    	
 	///////////////////////////////////////////////////////////////////////   TEST nbrs   	
	/*
	cout << " ******************************************************** Processor " << myrank 
	<< " with coords [" << coord[0] << "," << coord[1] << "] has " << endl <<
	"UP: " << nbrs[UP] << endl <<
	"DOWN: " << nbrs[DOWN] << endl <<
	"LEFT: " << nbrs[LEFT] << endl <<
	"RIGHT: " << nbrs[RIGHT] << endl << endl ;
	*/
	///////////////////////////////////////////////////////////////////////   TEST nrbs    	
    		


    	
    		
// ___________________________________________________________________________________________ SERIE    	
    	
    /*
    for (double dlt = dt; dlt < T; dlt += dt ){
    
    	// Filling the vorticity in the boundaries at time t
    		// ******* Bottom
    	for (int i = 0; i < Nx; i++){
    		v[i] = (s[i]-s[i+Nx])*2/pdlty;
    		}
    		// ******* Top
		for (int i = N-1; i > N-Nx+1; i--){
    		v[i] = (s[i]-s[i-Nx])*2/pdlty - 2*U/deltay;
    		}
    		// ******* Left
		for (int i = 0; i < Ny; i++){
    		v[i*Nx] = (s[i*Nx]-s[i*Nx+1])*2/pdltx;
    		// ******* Right
    		v[Nx-1+i*Nx] = (s[Nx-1+i*Nx]-s[Nx-2+i*Nx])*2/pdltx;
    		}
    	
    		
    	// Computing vorticity in the interior in processor MYRANK at time t
    	for (int i = 1; i < Ny-1; i++){
    		for (int j = 1; j < Nx-1; j++){
    				k = Nx*j+i;
    				v[k] = -(s[k+1]-2*s[k]+s[k-1])/pdltx-(s[k+Nx]-2*s[k]+s[k-Nx])/pdlty;
    			}
    		}
    	
    	// computing the interior vorticity at time t+dlt
    	for (int i = 1; i < Ny-1; i++){
    		for (int j = 1; j < Nx-1; j++){
    				k = Nx*j+i;
    				a = (s[k+1]-s[k-1])*(v[k+1]-v[k-1])/dlxy;
    				b = (s[k+Nx]-s[k-Nx])*(v[k+Nx]-v[k-Nx])/dlxy;
    				c = (v[k+1]-2*v[k]+v[k-1])/pdltx;
    				d = (v[k+Nx]-2*v[k]+v[k-Nx])/pdlty;
    				vt[k] = dt*(a - b + (c + d)/Re) + v[k];
    			}
    		}
    	
    	v = vt;
    	
    	// calling the function that solves Poisson in a separate class
    	SF.iterationsolver();
    		    		
      }
      */	

