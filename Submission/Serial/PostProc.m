clc
clear all
close all

% Initialising Grid Data
Nx = 161;
Ny = 161;


%% Post Processing Rountines for Re = 100

% Importing Streamfunction for Re = 100
S_100 = importdata('streamfunction_100.txt');
x = reshape(S_100.data(:,1), [Nx,Ny]);
y = reshape(S_100.data(:,2), [Nx,Ny]);
s_100 = reshape(S_100.data(:,3), [Nx,Ny]);

% Importing Vorticity for Re = 100
VORT_100 = importdata('vorticity_100.txt');
vort_100 = reshape(VORT_100.data(:,3), [Nx,Ny]);

% Generating velocity plot for Re = 100
U_100 = stream2U(s_100,Nx,Ny);
V_100 = stream2V(s_100,Nx,Ny);

% Finding minimum streamfunciton for Re = 100
minimum = min(min(s_100));
[x_min, y_min] = find(s_100==minimum);
disp("Minimum of SF at Re = 100: "+ minimum + " at x: "+ x(x_min,y_min) + "  y:" + y(x_min,y_min));

% Importing Streamfunction for Re = 400
S_400 = importdata('streamfunction_400.txt');
s_400 = reshape(S_400.data(:,3), [Nx,Ny]);


%% Post Processing Rountines for Re = 400

% Generating velocity plot for Re = 400
U_400 = stream2U(s_400,Nx,Ny);
V_400 = stream2V(s_400,Nx,Ny);

% Finding minimum streamfunciton for Re = 400
minimum = min(min(s_400));
[x_min, y_min] = find(s_400==minimum);
disp("Minimum of SF at Re = 400:" +minimum+ " at x: "+ x(x_min,y_min) + "  y:" + y(x_min,y_min));

%% Post Processing Rountines for Re = 1000

% Importing Streamfunction for Re = 1000
S_1000 = importdata('streamfunction_1000.txt');
s_1000 = reshape(S_1000.data(:,3), [Nx,Ny]);

% Generating velocity plot for Re = 1000
U_1000 = stream2U(s_1000,Nx,Ny);
V_1000 = stream2V(s_1000,Nx,Ny);

% Finding minimum streamfunciton for Re = 1000
minimum = min(min(s_1000));
[x_min, y_min] = find(s_1000==minimum);
disp("Minimum of SF at Re = 1000: "+minimum+ " at x: "+ x(x_min,y_min) + "  y:" + y(x_min,y_min));

%% Post Processing Rountines for Re = 3200

% Importing Streamfunction for Re = 3200
S_3200 = importdata('streamfunction_3200.txt');
s_3200 = reshape(S_3200.data(:,3), [Nx,Ny]);

% Generating velocity plot for Re = 3200
U_3200 = stream2U(s_3200,Nx,Ny);
V_3200 = stream2V(s_3200,Nx,Ny);

% Finding minimum streamfunciton for Re = 3200
minimum = min(min(s_3200));
[x_min, y_min] = find(s_3200==minimum);
disp("Minimum of SF at Re = 3200: "+minimum+ " at x: "+ x(x_min,y_min) + "  y:" + y(x_min,y_min));


%% Plotting Rountines 

% Plotting Streamfunction contour at Re = 100
figure;
contourf(x,y,s_100);
title("Streamfunction at Re = 100");
xlabel("x");
ylabel("y");
colorbar;
pbaspect([max(max(x)) max(max(y)) 1]);

% Plotting Vorticity contour at Re = 100
figure;
contourf(x,y,vort_100, 400, "LineColor", "none");
caxis([-4 4]);
title("Vorticity at Re = 100");
xlabel("x");
ylabel("y");
colorbar;
pbaspect([max(max(x)) max(max(y)) 1]);

% Plotting U(y) at x = 0.5
yy = linspace(0,1,Ny);
x5 = (Nx+1)/2;
figure;
plot(yy, U_100(x5,:), yy, U_400(x5,:), yy, U_1000(x5,:), yy, U_3200(x5,:));
legend("Re = 100", "Re = 400", "Re = 1000", "Re = 3200");
xlabel('y');
ylabel('U(y)');
title('U(y) along x = 0.5');
grid on;

% Plotting V(x) at y = 0.5
xx = linspace(0,1,Nx);
y5 = (Ny+1)/2;
figure;
plot(xx, V_100(:,y5), xx, V_400(:, y5), xx, V_1000(:, y5), xx, V_3200(:, y5));
legend("Re = 100", "Re = 400", "Re = 1000", "Re = 3200");
xlabel('x');
ylabel('V(x)');
title('V(x) along y = 0.5');
grid on;

