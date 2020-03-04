clf;clear;clc;
close all

% Problem Paramters
nu=0.01/10;                                            % Kinematic viscosity
ul=1;                                               % Lid velocity
% Grid Paramters
n=161;                                                % imax=jmax=n
h=1/(n-1);                                          % Grid step size
% Time March Paramters
mts = 1e5;                                            % Maximum number of time steps
dt=8e-4/10;                                             % Time step size
cts=0;                                              % Current time step
% SOR Paramters
b=1.5;                                              % Over-relaxation factor
mk=100;                                            % Maximum number of SOR iterations
% Declare and Initialize
[v,s]=deal(zeros(n,n));                             % Vorticity and stream function
% Time March Loop
while cts<mts
    % Compute the boundary nodes vorticity
    for i=2:n-1
        v(i,1)=2*(s(i,1)-s(i,2))/(h*h);             % vorticity on bottom wall
        v(i,n)=2*(s(i,n)-s(i,n-1))/(h*h)-2*ul/h;    % vorticity on top wall
        v(1,i)=2*(s(1,i)-s(2,i))/(h*h);             % vorticity on left wall
        v(n,i)=2*(s(1,i)-s(n-1,i))/(h*h);           % vorticity on right wall
    end
    % Compute the interior nodes vorticity
    for i=2:n-1
        for j=2:n-1
            v(i,j)=(4*s(i,j)-s(i+1,j)-s(i-1,j)-s(i,j+1)-s(i,j-1))/(h*h);
        end
    end
    % Compute the interior nodes new time step vorticity
    for i=2:n-1
        for j=2:n-1
            v(i,j)=v(i,j)+dt*((((s(i+1,j)-s(i-1,j))*(v(i,j+1)-v(i,j-1)))/(4*h*h))...
                -(((v(i+1,j)-v(i-1,j))*(s(i,j+1)-s(i,j-1)))/(4*h*h))...
                +(v(i+1,j)+v(i-1,j)+v(i,j+1)+v(i,j-1)-4*v(i,j))*(nu/(h*h)));
        end
    end
    % Compute the interior nodes new time step stream function
    ck=0;                                           % Current iteration number
    while ck<mk
        for i=2:n-1
            for j=2:n-1
                s(i,j)=(1-b)*s(i,j)+(s(i+1,j)+s(i-1,j)+s(i,j+1)+s(i,j-1)+v(i,j)*(h*h))*(b/4);
            end
        end
        ck=ck+1;
  
    
    end
    cts=cts+1
    % Increment current time
    % Prints every 1000 time-steps
    if(mod(cts,1000) == 0)
        cts
    end
end
% Time March Plot
    subplot(121), contour(rot90(fliplr(v))), axis('square');            % Vorticity
    subplot(122), contour(rot90(fliplr(s))), axis('square');            % Stream function

