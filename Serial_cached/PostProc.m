clc
clear all
close all

A = importdata('streamfunction.txt');
x = reshape(A.data(:,1), [161,161]);
y = reshape(A.data(:,2), [161,161]);
s = reshape(A.data(:,3), [161,161]);

contourf(x,y,s,10);



