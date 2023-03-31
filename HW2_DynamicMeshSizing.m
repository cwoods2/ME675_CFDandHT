%ME 675 - Homework 2, Problem 3
%Camden Woods - Spring 2023
clear, clc, close all

%% Part (b) - Dynamic Meshing
alpha = 1.05;   %Stretching factor
n = 128;        %Number of points
d = [0, 1];     %Domain boundaries

%Size of the first term ∆y1
dy1 = (1-alpha)/(1-alpha^n);

y1 = mesh_generation(dy1, alpha, d, n);

%Plotting points
figure(1)
scatter(y1, linspace(0,0,n))
xlabel('Grid Spacing')

%% Partc (c) - Boundary Layer
d99 = 0.176;
p = 0.7;
m = round(n*p);

%Setup for Newton-Raphson
x = []; x(1) = 1.05;   %Setting up solution vector and populating with first guess
err = 1; tol = 1e-8; ctr = 1;

while err > tol

    f = d99*(1-x(ctr)^n) - (1 - x(ctr)^m);
    fp = -d99*n*x(ctr)^(n-1) + m*x(ctr)^(m-1);
    
    x(ctr + 1) = x(ctr) - f/fp;
    
    ctr = ctr + 1;
    err_last = err;
    err = abs(x(ctr) - x(ctr - 1));
    
    %Exit conditions in case solution isn't converging (or converging slowly)
    if ctr > 1000
        fprintf("Number of iterations growing large. Try a better guess")
        break;
    elseif err_last < err
        fprintf("Solution seems to be nonconvergent. Try a better guess\n")
        fprintf("Current Error: %6.5f, Last Error: %6.5f", err, err_last)
        break;
    end   
end

%Size of the first term ∆y1
first_term = (1-x(end))/(1-x(end)^n);

y2 = mesh_generation(first_term, x(end), d, n);

% pnt_ctr = 0;
% iter = 2;
% sum = 0;
% while sum <= d99
%     pnt_ctr = pnt_ctr + 1;
%     sum = sum + (y2(iter) - y2(iter-1));
%     iter = iter + 1;
% 
% end

% Plotting points
figure(3)
hold on
scatter(y1, linspace(0.1,0.1,n))
scatter(y2, linspace(-0.1,-0.1,n), '+')
xline(d99)
hold off
legend('Standard Mesh', 'd99 Mesh')
xlabel('Grid Spacing')

function y = mesh_generation(dy1, alpha, d, n)
%This function will generate a non-uniform mesh
%INPUTS:
%   dy1: Size of the first subunit
%   alpha: scaling parameter
%   d: start and end points of the mesh
%   n: number of points
    %Generating mesh
    y = []; y(1) = d(1);
    dy = []; dy(1) = dy1;

    for i = 2:n
        y(i) = alpha*dy(i-1) + y(i-1);
        dy(i) = dy(i-1)*alpha;
    end
end
