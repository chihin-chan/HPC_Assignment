clc
clear all
close all

Nx = 81;
Ny = 81;

A = importdata('streamfunction.txt');
x = reshape(A.data(:,1), [Nx,Ny]);
y = reshape(A.data(:,2), [Nx,Ny]);
s = reshape(A.data(:,3), [Nx,Ny]);

contourf(x,y,s,20);
pbaspect([max(max(x)) max(max(y)) 1]);



