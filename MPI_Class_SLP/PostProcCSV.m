clear all
clc

Lx = 1;
Ly = 1;
Nx = 16;
Ny = 16;

[X, Y] = meshgrid(linspace(0,Lx,Nx), linspace(0,Ly,Ny));

V = readmatrix("vorticity.txt");
S = readmatrix("streamfunction.txt");

contourf(S, 20);
axis ij