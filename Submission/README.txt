Folders:

1. Serial
    Stores the serial code ande used Matlab to post process results for various reynolds number.

2. SerialLxLy
    Stores the serial code and use Matlab to post process results for different Lx Ly.

3. Parallel
    Stores the parallelised code by using Lapack's dpbtrf/dpbtrs and solving Poisson locally in each
    subdomains independently.

4. ParallelScaLapack
    Stores the parallelised code by using ScaLapack pdpbtrf/pdpbtrs and solving the Poisson equation
    globally. However, only works for Nx, Ny < 50
