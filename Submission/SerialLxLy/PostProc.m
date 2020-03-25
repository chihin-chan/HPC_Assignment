clc
clear all
close all

Nx = 161;
Ny = 161;

A = importdata('streamfunction_Lx1Ly2.txt');
x_1 = reshape(A.data(:,1), [Nx,Ny]);
y_1 = reshape(A.data(:,2), [Nx,Ny]);
s_1 = reshape(A.data(:,3), [Nx,Ny]);

B = importdata('streamfunction_Lx2Ly1.txt');
x_2 = reshape(B.data(:,1), [Nx,Ny]);
y_2 = reshape(B.data(:,2), [Nx,Ny]);
s_2 = reshape(B.data(:,3), [Nx, Ny]);

figure;
contourf(x_1, y_1, s_1);
title('Lx = 1, Ly = 2, Re = 100');
xlabel('x');
ylabel('y');
colorbar;
pbaspect([max(max(x_1)) max(max(y_1)) 1]);

figure;
contourf(x_2, y_2, s_2);
title('Lx = 2, Ly = 1, Re = 100');
xlabel('x');
ylabel('y');
colorbar;
pbaspect([max(max(x_2)) max(max(y_2)) 1]);

