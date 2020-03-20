/*
 *  High-performance Computing - COURSEWORK
 *
 *  Poisson Solver
 *
 *  Author: Lidia Caros Roca, CID: 01821851
 *
 */

// pre-processor
#include "PoissonSolver.h"
#include "LidDrivenCavity.h"

#include <iostream>
using namespace std;

#include <mpi.h>
#include <math.h>
#include <iomanip>

//	BLAS library
#include "cblas.h"

// 	ScaLAPACK and BLACS
#define F77NAME(x) x##_
extern "C" {         
	   
	                  		  
	/* Cblacs declarations */
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
	 
} 

PoissonSolver::PoissonSolver()
{
}

PoissonSolver::~PoissonSolver()
{
	cout	<<	"PoissonSolver destructed" << endl;
}


//*****************************   INITIALIZING  ********************************************//

void PoissonSolver::Initializing(LidDrivenCavity& LDC)
{
	// Variable definition ------------------------------------------------------------------

	// Internal points per processor
	innernx = LDC.lastpointx - LDC.firstpointx;		// Number of points in x-dir in A
	innerny = LDC.lastpointy - LDC.firstpointy;		// Number of points in y-dir in A
	innern = innernx*innerny;						// Order of the matrix A
	
	// Max number of points in processors
	N = LDC.nx*LDC.ny;
					
	// BLACS
	npe		=	LDC.Px*LDC.Py;
	mype	=	LDC.myrank;
	ncol = npe;
	nrow = 1;

	// Dimensions
	NA	 = (LDC.Nx-2)*(LDC.Ny-2);		// Problem size
	NB	 = ceil(  (double)NA/npe );		// Blocking size				
	BWU	 = (LDC.Nx-2)-1;				// Upper bandwith
	nrhs = 1; 							// Number of rhs to solve
	JA	 = 1;							// Start offset in matrix
	IB	 = 1;							// Start offset in vector
	lda	 = 1+BWU;						// Leading dimension A
	LA	 = lda*NB;						// Local banded size
	LAF	 = (NB+2*BWU)*BWU;				// Size of Auxiliary Fillin
	lworkf	= BWU*BWU;					// Size of workspace workf (facto.)
	lworks	= nrhs*BWU;					// Size of workspace works (solve)

	// Arrays
	A 	 = new double [LA];				// Symmetric banded matrix storage
	AF 	 = new double [LAF];			// Auxiliary Fillin 
	ipiv = new int	  [NB];				// Pivoting array
	b	 = new double [NB];				// In: RHS vector, Out: Solution
	workf =	new double[lworkf];			// Workspace factorization
	works =	new double[lworks];			// Workspace solving	
	
	// Initializing A and b;
	for (int i = 0; i < LA; ++i){ A[i] = 0.0;}
	for (int i = 0; i < NB; ++i){ b[i] = 0.0;}
		
	// b with cartesian distribution
	bcart = new double [N];
	for (int i = 0; i < N    ; i++){ bcart[i]    = 0;}


	// CBLACS ------------------------------------------------------------------------------		
		
	// Init BLACS world communicator
	Cblacs_pinfo(&mype, &npe);
	
	// Get the default system context (i.e. MPI_COMM_WORLD)
	Cblacs_get( 0, 0, &ctx );
	
	// Initialise a process grid of 1 rows and npe columns
	Cblacs_gridinit( &ctx, "Row_Major", 1, npe );
	
	// Get info about the grid to verify it is set up
	Cblacs_gridinfo( ctx, &nrow, &ncol, &myrow, &mycol);
	
	
	// DESCRIPTORS -------------------------------------------------------------------------	

	// Descriptor for A	
	desca = new int [7];
	desca[0] = 501;				// Type banded matrix
	desca[1] = ctx;				// BLACS context handle
	desca[2] = NA;				// Problem size
	desca[3] = NB;				// Blocking of matrix
	desca[4] = 0;				// Process row/column
	desca[5] = lda;			 	// Local leading dim
	desca[6] = 0;				// Reserved for future use
	
			
	// Descriptor for b
	descb = new int [7];
	descb[0] = 502;				// Type banded RHS
	descb[1] = ctx;				// BLACS context handle
	descb[2] = NA;				// Problem size
	descb[3] = NB;				// Blocking of matrix
	descb[4] = 0;				// Process row/column
	descb[5] = NB;				// Local leading dim
	descb[6] = 0;				// Reserved for future use
	
	// GATHER -------------------------------------------------------------------------------
	
	if (LDC.myrank == 0){	
		innernxV = new int    [npe];
		innernyV = new int    [npe];
		for (int i = 0; i < npe  ; i++){ innernxV[i] = 0; 
									     innernyV[i] = 0;}
	}		
				
	
	MPI_Gather(&innernx, 1, MPI_INT,    innernxV, 1, MPI_INT,    0, LDC.cart_comm);	
	MPI_Gather(&innerny, 1, MPI_INT,    innernyV, 1, MPI_INT,    0, LDC.cart_comm); 
	
	// root processor -----------------------------------------------------------------------
	
	if (LDC.myrank == 0){	
	
	
		// Gathering vectors
		bcartnpe  = new double [N*npe];		
		new_inNV     = new int [npe];
		new_innernyV = new int [LDC.Py];
		new_innernxV = new int [LDC.Px];
		
		// b with cartesian distribution
		bscanpe = new double [NB*npe];	
					
		// Initializing vectors
		for (int i = 0; i < N*npe; i++){ bcartnpe[i]  = 0;}		
		for (int i = 0; i < NB*npe; ++i){ bscanpe[i] = 0.0;}
		for (int i = 0; i < LDC.Px; ++i){ new_innernxV[i] = 0.0;}
		for (int i = 0; i < LDC.Py; ++i){ new_innernyV[i] = 0.0;}
						
		// summation of all previous processor points
		for (int i = 0; i < npe; i++){
			new_inNV[i] = i*N;
		}
				
		// summation of all previous in x
		new_innernxV[0] = 0;
		for (int i = 1; i < LDC.Px; i++){
			for (int j = i-1; j > -1; j--){
				new_innernxV[i] = new_innernxV[i] + innernxV[j]; 			
			}
		}
		// summation of all previous in y
		new_innernyV[0] = 0;
		for (int i = 1; i < LDC.Py; i++){
			for (int j = i-1; j > -1; j--){
				new_innernyV[i] = new_innernyV[i] + innernyV[j*LDC.Px]; 	
			}
		}
	}		
}


//*******************  A MATRIX CONSTRUCTION AND FACTORIZATION  ****************************//

void PoissonSolver::MatrixConstr(LidDrivenCavity& LDC){
	
	// Filling the matrix packed storage ---------------------------------------------------
		
	// Constant values in A
	diag = 2*(1/LDC.pdltx+1/LDC.pdlty);		// Diagonal
	offdiagi = -1/LDC.pdltx;				// Upper Diagonal i+1
	offdiagj = -1/LDC.pdlty;				// Upper Diagonal j+1
	
	int ki = LDC.myrank*NB;	
	if ((LDC.myrank*NB)%(LDC.Nx-2)==0){ 
		ki = LDC.Nx-2; 
	} else { 
		int ko = (int)(LDC.myrank*NB-1)/(LDC.Nx-2);
		ki = LDC.myrank*NB - ko*(LDC.Nx-2); 
	}
	
	for (int i = 0; i < NB; ++i) {
        if( LDC.myrank*NB+i < NA){ 
			// Filling first upper off-diagonal
			if (ki%(LDC.Nx-2)==0){
				A[i*lda + BWU-1] = 0; // Sides and 1st point
				ki = 1;
			} else {
			   A[i*lda + BWU-1] = offdiagi;
			   ki++;
			}
		    // Filling diagonal terms
		    A[i*lda + BWU] = diag; 
		}
    }       
    // Filling second upper off-diagonal;    
    int kp=0;
   	int index = LDC.Nx-2 - LDC.myrank*NB;
   	if (index<0) { index = 0; }
    for (int i = index; i < NB; ++i) {
    	if( LDC.myrank*NB+kp < NA){
		    A[i*lda] = offdiagj;  
		    kp++;
		}  
    }
    
    if (LDC.myrank == LDC.rankdisplay){	
		//cout<< endl << "____________________________________ Time step " << dlt;	
		cout << endl << "Processor " << LDC.myrank << endl
		<< " MATRIX A: " << endl << endl;
		for (int i = 0; i < lda; i++){
			for (int j = 0; j < NB; j++){
				cout <<  setw(12) <<A[lda*j+i] << "	";
			} cout << endl;
		} cout << endl;
	}
    
    
    // A factorization ----------------------------------------------------------------------

    Cblacs_barrier( ctx, "All");
        
    //  'U':  Upper triangle of A is stored
    F77NAME(pdpbtrf)('U', NA, BWU, A, JA, desca, AF, LAF, workf, lworkf, &info);
   
    // INFO
    if (info<0) {
        cout << "Error in the LU factorization. INFO: " << info << endl;
    } else if (info>0){
    	cout << " the submatrix stored on processor " << LDC.myrank 
    	<< " and factored locally was not positive definite. INFO: " << info << endl;
    } else {
        cout << "Factorization completed on processor " << LDC.myrank << endl; 
    }
    
    if (LDC.myrank == LDC.rankdisplay){	
		//cout<< endl << "____________________________________ Time step " << dlt;	
		cout << endl << "Processor " << LDC.myrank << endl
		<< " MATRIX FACTO: " << endl << endl;
		for (int i = 0; i < lda; i++){
			for (int j = 0; j < NB; j++){
				cout <<  setw(12) <<A[lda*j+i] << "	";
			} cout << endl;
		} cout << endl;
	}

}


//************************** FILLING THE RIGHT HAND SIDE VECTOR ****************************//

void PoissonSolver::FillingRHS(LidDrivenCavity& LDC){

	// Variables for indexing bcart[] and v[]
	int z = 0;
	int p = 0;
	for (int i = 0; i < N    ; i++){ bcart[i]    = 0;}
	
	// Filling the vector
	for (int j = LDC.firstpointy; j < LDC.lastpointy; j++){
		for (int i = LDC.firstpointx; i < LDC.lastpointx; i++){
			z = LDC.Npx*j+i;
			bcart[p] = LDC.v[z]; 
			p++;
		}
	}
	
	// Adding the boundaries -----------------------------------------------------------
	
	p = 0;
	// Boundary down
	if (LDC.nbrs[LDC.DOWN] == -2){
		for (int i = LDC.firstpointx; i < LDC.lastpointx; i++){
			z = LDC.Npx*LDC.firstpointy+i;
			bcart[p] = bcart[p] + LDC.s[z-LDC.Npx]/LDC.pdlty; 
			p++;
		}
	}
	p = 0;
	// Boundary left
	if (LDC.nbrs[LDC.LEFT] == -2){
		for (int j = LDC.firstpointy; j < LDC.lastpointy; j++){
			z = LDC.Npx*j+LDC.firstpointx;
			bcart[p] = bcart[p] + LDC.s[z-1]/LDC.pdltx;
			p = p + innernx;
		}
	}
	p = innernx-1;
	// Boundary right
	if (LDC.nbrs[LDC.RIGHT] == -2){
		for (int j = LDC.firstpointy; j < LDC.lastpointy; j++){
			z = LDC.Npx*j+(LDC.lastpointx-1);
			bcart[p] = bcart[p] + LDC.s[z+1]/LDC.pdltx;
			p = p + innernx;
		}
	}
	p  = innernx*(innerny-1);
	// Boundary up
	if (LDC.nbrs[LDC.UP] == -2){
		for (int i = LDC.firstpointx; i < LDC.lastpointx; i++){
			z = LDC.Npx*(LDC.lastpointy-1)+i;
			bcart[p] = bcart[p] + LDC.s[z+LDC.Npx]/LDC.pdlty; 
			p++;
		}
	}
	
}


//******************     CARTESIAN TO SCALAPACK RHS VECTOR MAPPING     *********************//

void PoissonSolver::CartesianToScalapack(LidDrivenCavity& LDC){
	
	
	// From root processor
	if (LDC.myrank == 0){
				
		for (int procy = 0;	procy < LDC.Py; procy++){			// Processors in Y-dir
			for (int procx = 0; procx < LDC.Px; procx++){ 		// Processors in X-dir
				
				// local number of points in current processor c_proc
				cproc = procx+procy*LDC.Px; 
				localny = innernyV[cproc];
				localnx = innernxV[cproc];
				
				// number of points behind the current processor
				prevny = new_innernyV[procy];
				prevnx = new_innernxV[procx];
				prevn = new_inNV[cproc];
				
				//cout << " cproc = " << cproc << " [procx,procy] = [" << procx << "," << procy 
				//<< "]" << " ---------- prevn = " << prevn << " prevnx: " << prevnx 
				//<< ", prevny: " << prevny << endl;  
				
				
				// loop inside processor				
				for (int j = 0; j < localny; j++){				// Y-dir inside proc.
					for (int i = 0; i < localnx; i++){			// X-dir inside proc.
						
						// index
						k1 = i + j*(LDC.Nx-2) + prevnx + prevny*(LDC.Nx-2);
						k2 = prevn + i + j*localnx;
						
						//cout << " k1 = " << k1 << ", k2 = " << k2 <<  endl;
						
						bscanpe[k1] = bcartnpe[k2];
						
					}
				}
			}	
		}
		
	} // end of if root

}



//******************     SCALAPACK TO CARTESIAN RHS VECTOR MAPPING     *********************//

void PoissonSolver::ScalapackToCartesian(LidDrivenCavity& LDC){

	// From root processor
	if (LDC.myrank == 0){
				
		for (int procy = 0;	procy < LDC.Py; procy++){			// Processors in Y-dir
			for (int procx = 0; procx < LDC.Px; procx++){ 		// Processors in X-dir
				
				// local number of points in current processor c_proc
				cproc = procx+procy*LDC.Px; 
				localny = innernyV[cproc];
				localnx = innernxV[cproc];
				
				// number of points behind the current processor
				prevny = new_innernyV[procy];
				prevnx = new_innernxV[procx];
				prevn = new_inNV[cproc];
				
				//cout << " cproc = " << cproc << " [procx,procy] = [" << procx << "," << procy 
				//<< "]" << " ---------- prevn = " << prevn << " prevnx: " << prevnx 
				//<< ", prevny: " << prevny << endl;  
				
				
				// loop inside processor				
				for (int j = 0; j < localny; j++){				// Y-dir inside proc.
					for (int i = 0; i < localnx; i++){			// X-dir inside proc.
						
						// index
						k1 = i + j*(LDC.Nx-2) + prevnx + prevny*(LDC.Nx-2);
						k2 = prevn + i + j*localnx;
						
						bcartnpe[k2] = bscanpe[k1];
						
					}
				}
			}	
		}
		
	}
}




//*********************************   SOLVER   ********************************************//


void PoissonSolver::Solver(LidDrivenCavity& LDC){ 

	// Filling the right hand side vector bcart ---------------------------------------------

	FillingRHS(LDC);
		
	
	//for (int i = 0; i < N; ++i){ cout << "PROCESSOR " << LDC.myrank << " solution bcart[" << i 
	//							  << "] = " << bcart[i] << endl;}
								  	
	// GATHER ------------------------------------------------------------------------------- 
		
	// Barrier before gathering
	MPI_Barrier(LDC.cart_comm);
	
	// Gather MPI functions								 root = 0	
	MPI_Gather(bcart,    N, MPI_DOUBLE, bcartnpe, N, MPI_DOUBLE, 0, LDC.cart_comm);		
	
		
	// Filling the rhs vector b in Scalapack form -------------------------------------------- 
	
	CartesianToScalapack(LDC);
	
	/*
	if (LDC.myrank == 0){
	 for (int i = 0; i < NB*npe; ++i){ cout << "PROCESSOR " << LDC.myrank << " solution bscanpe[" 
	 << i  << "] = " << bscanpe[i] << endl;}
	}
	*/
	
				
	// SCATTER -------------------------------------------------------------------------------	
	
	// Barrier before scattering
	MPI_Barrier( LDC.cart_comm );
	
	// Scatter vector b from root = 0, to the other processors	
	MPI_Scatter(&bscanpe[NB*LDC.myrank], NB, MPI_DOUBLE, b, NB, MPI_DOUBLE, 0, LDC.cart_comm);	
	
	//for (int i = 0; i < NB; ++i){ cout << "PROCESSOR " << LDC.myrank << " solution b[" << i 
	//							  << "] = " << b[i] << endl;}
			
		
	// SOLVING WITH SCALAPACK ----------------------------------------------------------------
		
	//  'U':  Upper triangle of A is stored  
	F77NAME(pdpbtrs)('U', NA, BWU, nrhs, A, JA, desca, b, IB, descb, AF, LAF, works, lworks, &info);
	
	if (info) {
	    cout << "Error solving the linear system. Info:  " << info << endl;
	} else {
		cout << "Solving completed on processor " << LDC.myrank << endl;
	}
	
	//for (int i = 0; i < NB; ++i){ cout << "PROCESSOR " << LDC.myrank << " solution b[" << i 
	//							  << "] = " << b[i] << endl;}
	
	// GATHER ------------------------------------------------------------------------------- 
			
	// Gather MPI function getting chunks of size [NB] in processors in an ascendent way
	MPI_Gather(b, NB, MPI_DOUBLE, bscanpe, NB, MPI_DOUBLE, 0, LDC.cart_comm);
	
	// Returning to cartesian grid -----------------------------------------------------------
	
	ScalapackToCartesian(LDC);
	
	// SCATTER -------------------------------------------------------------------------------
	
	// Scatter MPI function storing chunks of size [N] in processors in an ascendent way
	MPI_Scatter(&bcartnpe[N*LDC.myrank], N, MPI_DOUBLE, bcart, N, MPI_DOUBLE, 0, LDC.cart_comm);	
	
	// Storing in streamfunction -------------------------------------------------------------
		
	q = 0, z = 0;
	for (int j = LDC.firstpointy; j < LDC.lastpointy; j++){
		for (int i = LDC.firstpointx; i < LDC.lastpointx; i++){
			z = LDC.Npx*j+i;
			LDC.s[z] = bcart[q]; 
			q++;
		}
	}
	
	/*
	// Printing the solution ----------------------------------------------------------------- 
	if (LDC.myrank == LDC.rankdisplay){	
		cout << endl << "Processor " << LDC.myrank << endl
		<< " MATRIX SOLUTION s: " << endl << endl;
		for (int j = LDC.Npy-1; j > -1; j--){
			for (int i = 0; i<LDC.Npx; i++){
				cout <<  setw(12) << LDC.s[LDC.Npx*j+i] << "	";
			} cout << endl;
		} cout << endl;
	}
	*/
	
	
	// Finalising BLACS ----------------------------------------------------------------------
	//Cblacs_gridexit( ctx );
	
		
}	
	
	
	
	
	/*	////////////////////////////////////////////////////////////

	// root processor CHECK GATHER
	/*
	if (LDC.myrank == 0){ 
		for (int i = 0; i < N; i++){
			cout << "--- proc " << LDC.myrank << " bcart[" << i << "] = " << bcart[i] << endl;
		}
		
		for (int i = 0; i < N*npe; i++){
			cout << "bcartnpe[" << i << "] = " << bcartnpe[i] << endl;
		}
		for (int i = 0; i < npe; i++){		
			cout << "innernxV[" << i << "] = " << innernxV[i] << endl;
			cout << "innernyV[" << i << "] = " << innernyV[i] << endl;
		}	
	}
	*/	
	/*
	// Printing the vector ------------------------------------------------------------------  
	if (LDC.myrank == LDC.rankdisplay){   
		cout << endl << " POISSON SOLVER : Vector b " << endl << endl;    
	 	for (int i = 0; i < innern; i++){
			cout << "Rank : "  << LDC.myrank << "     b["<< i << "] = " << b[i] << endl;
		}
	} 
	
	// Printing the solution ----------------------------------------------------------------- 
	if (LDC.myrank == LDC.rankdisplay){	
		cout << endl << "Processor " << LDC.myrank << endl
		<< " MATRIX SOLUTION s: " << endl << endl;
		for (int j = LDC.Npy-1; j > -1; j--){
			for (int i = 0; i<LDC.Npx; i++){
				cout <<  setw(12) << LDC.s[LDC.Npx*j+i] << "	";
			} cout << endl;
		} cout << endl;
	}
	*/

		 // PRINTING MATRIX
	/*
	if (LDC.myrank == LDC.rankdisplay){	
		//cout<< endl << "____________________________________ Time step " << dlt;	
		cout << endl << "Processor " << LDC.myrank << endl
		<< " MATRIX v POISSON: " << endl << endl;
		for (int j = LDC.Npy-1; j > -1; j--){
			for (int i = 0; i<LDC.Npx; i++){
				cout <<  setw(12) << LDC.v[LDC.Npx*j+i] << "	";
			} cout << endl;
		} cout << endl;
	}
	
	if (LDC.myrank == LDC.rankdisplay){	
		//cout<< endl << "____________________________________ Time step " << dlt;	
		cout << endl << "Processor " << LDC.myrank << endl
		<< " MATRIX s: " << endl << endl;
		for (int j = LDC.Npy-1; j > -1; j--){
			for (int i = 0; i<LDC.Npx; i++){
				cout <<  setw(12) << LDC.s[LDC.Npx*j+i] << "	";
			} cout << endl;
		} cout << endl;
	}*/	
	
	
	 	////////////////////////////////////////////////////////////   TEST adding boundaries   	
	/*
	
	// Adding the boundaries
		p = 0;
		// Boundary down
		for (int i = LDC.firstpointx; i < LDC.lastpointx; i++){
			z = LDC.Npx*LDC.firstpointy+i;
			b[p] = b[p] + LDC.s[z-LDC.Npx]/LDC.pdlty; 
			if (LDC.myrank == LDC.rankdisplay) cout << " down: b[" << p << "] = s[" << z-LDC.Npx << "]" << endl;
			p++;
		}
		p = 0;
		// Boundary left
		for (int j = LDC.firstpointy; j < LDC.lastpointy; j++){
			z = LDC.Npx*j+LDC.firstpointx;
			b[p] = b[p] + LDC.s[z-1]/LDC.pdltx;
			if (LDC.myrank == LDC.rankdisplay) cout << " left: b[" << p << "] = s[" << z-1 << "]" << endl;		
			p = p + innernx;
		}
		p = innernx-1;
		// Boundary right
		for (int j = LDC.firstpointy; j < LDC.lastpointy; j++){
			z = LDC.Npx*j+(LDC.lastpointx-1);
			b[p] = b[p] + LDC.s[z+1]/LDC.pdltx;
			if (LDC.myrank == LDC.rankdisplay){
				cout << " right: b[" << p << "] = s[" << z+1 << "]" << endl;}
			p = p + innernx;
		}
		p  = innernx*(innerny-1);
		// Boundary up
		for (int i = LDC.firstpointx; i < LDC.lastpointx; i++){
			z = LDC.Npx*(LDC.lastpointy-1)+i;
			b[p] = b[p] + LDC.s[z+LDC.Npx]/LDC.pdlty; 
			if (LDC.myrank == LDC.rankdisplay) {
				cout << " up: b[" << p << "] = s[" << z+LDC.Npx << "]" << endl;}
			p++;
		}
	*/
	/*
		// Printing the vector 
		if (LDC.myrank == LDC.rankdisplay){   
			cout << endl << " POISSON SOLVER : Vector s " << endl << endl;    
		 	for (int i = 0; i < innern; i++){
				cout << "Rank : "  << LDC.myrank << "     s["<< i << "] = " << b[i] << endl;
			}
		} 
		*/
	///////////////////////////////////////////////////////////////////////   TEST   


