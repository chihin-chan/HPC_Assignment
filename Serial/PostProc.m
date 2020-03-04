clc
clear all
close all

Nx = 161;
Ny = 161;

A = importdata('streamfunction.txt');
x = reshape(A.data(:,1), [Nx,Ny]);
y = reshape(A.data(:,2), [Nx,Ny]);
s = reshape(A.data(:,3), [Nx,Ny]);

contourf(x,y,s,100);



