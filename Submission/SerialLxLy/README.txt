How to run the code?

Use make to compile the code

1. $ make

Run it with in serial

2. $ ./myprog --Lx 2 --Lx 1 --Nx 161 --Ny 161 --Px 1 --Py 1 --dt 0.0009 --T 30 --Re 100

Repeat step 2 with different Lx and Ly
Lx = 1		Ly = 2

Run Matlab Postprocessing routine to generate streamfunction

3. $ matlab PostProc.m

