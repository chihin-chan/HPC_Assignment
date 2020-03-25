How to run the code?

Use make to compile the code

1. $ make

Run it with mpiexec, example of 4 proccessors is used here

2. $ mpiexec -np 4 ./myprog --Lx 1 --Lx 1 --Nx 50 --Ny 50 --Px 2 --Py 2 --dt 0.0009 --T 30 --Re 100

