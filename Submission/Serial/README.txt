How to run the code?

Use make to compile the code

1. $ make

Run it with in serial

2. $ ./myprog --Lx 1 --Lx 1 --Nx 161 --Ny 161 --Px 1 --Py 1 --dt 0.0009 --T 30 --Re 100

Repeat step 2 with different dt and reynolds number
dt = 0.003      Re = 400
dt = 0.005      Re = 1000
dt = 0.003      Re = 3200

Run Matlab Postprocessing routine

3. $ matlab PostProc.m

Note:
    stream2U.m -> Converts streamfunction to velocity field U
    stream2V.m -> Converts streamfunction to velocity field V
