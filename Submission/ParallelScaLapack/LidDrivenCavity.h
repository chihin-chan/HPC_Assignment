#pragma once

#include <string>
#include <mpi.h>
using namespace std;

class LidDrivenCavity
{
public:
    LidDrivenCavity();
    ~LidDrivenCavity();

    void SetDomainSize(double xlen, double ylen);
    void SetGridSize(int nx, int ny);
	void SetPartitionSize(int ppx, int ppy);
    void SetTimeStep(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double Re);
	void GetRank(int rr);
	void GetSize(int ss);
	void SMatPrintRank(int r);
	void VMatPrintRank(int r);


    void Initialise();
	void BoundaryConditions();
	void InteriorUpdate();
	void Communicate();
    void Integrate();
    void ExportSol();
	void MapRHS();
	void CholFact();
	void ScalInit();
	void CartToScal();
	void ScalToCart();
	void MapSF();
    // Add any other public functions

	friend class PoissonSolver;

private:
    double* v = nullptr;
    double* s = nullptr;
    double* rhs = nullptr;
    double* v_new = nullptr;
	
	// Communication Pointers
    double* outbuf_s_D = nullptr;
    double* outbuf_v_D = nullptr;
    double* inbuf_s_D = nullptr;
    double* inbuf_v_D = nullptr;
    double* outbuf_s_U = nullptr;
    double* outbuf_v_U = nullptr;
    double* inbuf_s_U = nullptr;
    double* inbuf_v_U = nullptr;
    double* outbuf_s_R = nullptr;
    double* outbuf_v_R = nullptr;
    double* inbuf_s_R = nullptr;
    double* inbuf_v_R = nullptr;
    double* outbuf_s_L = nullptr;
    double* outbuf_v_L = nullptr;
    double* inbuf_s_L = nullptr;
    double* inbuf_v_L = nullptr;

	double U;
    double dt;
    double T;
    int    Nx;
    int    Ny;
	int    loc_nx;
	int    loc_ny;
    double Lx;
    double Ly;
    double Re;
	double dx;
	double dy;
	
	// MPI Variables
	int Px;
	int Py;
    int size;
    int rank;
	int nghbrs[4];
    int dims[2];
    int periods[2] = {0, 0};
	int coords[2];
	int tag[2] = {0,1};
	MPI_Comm mygrid;

	// Scalapack variables
	int N;
	int NB;
	int BW;
	int lda;
	int nrhs;
	int JA;
	int IB;
	int LA;
	int LAF;
	int lworkf;
	int lworks;
	int info;
	int myrow;
	int mycol;
	int ncol;
	int nrow;
	int ctx;
	
	
	// Arrays for scalapack
	double* A;	// Stores banded matrix A in blocks
	double* AF;
	double* b;
	double* b_rank;
	double* b_cart;
	double* b_scal;
	double* workf;
	double* works;
	int* desca;
	int* descb;
	int* ipiv;
	int* b_scal_nx;
	int* b_scal_ny;
	int b_nx;
	int b_ny;
	int b_rank_L;
	int b_scalx;
	int b_scaly;
	int x_begin;
	int x_end;
	int y_begin;
	int y_end;
	int counter;
	
};

