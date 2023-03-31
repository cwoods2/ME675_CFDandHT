%ME 675 - Homework 3
% Camden Woods - Spring 2023
clear, clc, close all

%% Problem Setup
%Given params
alpha = 1e-3;

%Create uniform mesh
nx=50;                         % Number of grid points in x
ny=nx;                         % Number of grid points in y
x=linspace(0,1,nx+1);          % x-Grid (location of cell faces)
y=linspace(0,1,ny+1);          % y-Grid (location of cell faces)
xm=x(1:end-1)+(x(2)-x(1))/2;   % x-Grid (location of cell centers)
ym=y(1:end-1)+(y(2)-y(1))/2;   % y-Grid (location of cell centers)
dx=xm(2)-xm(1);                % Grid spacing in x
dy=ym(2)-ym(1);                % Grid spacing in y

% Store velocities at faces
u=zeros(nx+1,ny);
v=zeros(nx,ny+1);
for i=1:nx+1
    for j=1:ny
        u(i,j)=1/10-(sin(pi*x(i)))^2*(sin(pi*(ym(j)-0.05))*...
            cos(pi*(ym(j)-0.05))-sin(pi*(ym(j)+0.05))*cos(pi*(ym(j)+0.05)));
    end
end
for i=1:nx
    for j=1:ny+1
        v(i,j)=sin(pi*xm(i))*cos(pi*xm(i))*((sin(pi*(y(j)-0.05)))^2-...
            (sin(pi*(y(j)+0.05)))^2);
    end
end
u_old = u;
v_old = v;


% Generate RHS
b=zeros(nx*ny,1);
for i=1:nx
    for j=1:ny
        % Lexicographic ordering
        n=i+nx*(j-1);
        % Source (store at cell center)
        b(n)=exp(-((xm(i)-0.75)^2+(ym(j)-0.5)^2)/(0.05^2))-exp(-((xm(i)-0.25)^2+(ym(j)-0.5)^2)/(0.05^2));
    end
end

% %% Plot source and velocity fields
% figure(1)
% bb=transpose(reshape(b,[nx ny]));
% contour(bb)
% imagesc(bb)
% set(gca,'ydir','normal')
% title('b plot')
% colorbar

% figure(2)
% uu=transpose(reshape(u,[nx+1 ny]));
% contour(uu)
% imagesc(uu)
% set(gca,'ydir','normal')
% title('u plot')
% 
% figure(3)
% vv=transpose(reshape(v,[nx ny+1]));
% contour(vv)
% imagesc(vv)
% set(gca,'ydir','normal')
% title('v plot')

%% Part A - First Order Upwinding
A_PtA = MeshGen_Upwind(nx, ny, u, v, dx, dy, alpha);
phi_A = mldivide(A_PtA, b);

%% Part B - Second Order Linear Reconstruction
A_PtB = MeshGen_2ndO(nx, ny, u, v, dx, dy, alpha);
phi_B = mldivide(A_PtB,b);

%% Part C - QUICK Scheme
A_PtC = MeshGen_QUICK(nx, ny, u, v, dx, dy, alpha);
phi_C = mldivide(A_PtC,b);

%% Plotting Results
% Plot Answer to A
figure(4)
AA=transpose(reshape(phi_A,[nx ny]));
AA=AA-AA(nx/2,1);
hold on
contourf(xm,ym,AA)
hs=streamslice(xm,ym,u(1:nx,1:ny)',v(1:nx,1:ny)');
set(hs,'color','w','linewidth',1)
axis square
title('$\phi$: Upwind Method','interpreter','latex','fontsize',12)
set(gcf, 'Color',[1,1,1]);
hold off
colorbar

% Plot Answer to B
figure(5)
BB=transpose(reshape(phi_B,[nx ny]));
BB=BB- BB(nx/2,1);
hold on
contourf(xm,ym,BB)
hs=streamslice(xm,ym,u(1:nx,1:ny)',v(1:nx,1:ny)');
set(hs,'color','w','linewidth',1)
axis square
title('$\phi$: 2nd Order','interpreter','latex','fontsize',12)
set(gcf, 'Color',[1,1,1]);
hold off
colorbar

% Plot Answer to C
figure(6)
CC=transpose(reshape(phi_C,[nx ny]));
CC=CC-CC(nx/2,1);
hold on
contourf(xm,ym,CC)
hs=streamslice(xm,ym,u(1:nx,1:ny)',v(1:nx,1:ny)');
set(hs,'color','w','linewidth',1)
axis square
title('$\phi$: QUICK Method','interpreter','latex','fontsize',12)
set(gcf, 'Color',[1,1,1]);
colorbar
hold off

% Plot Answer to A
figure(7)
tiledlayout(1,3)
nexttile
AA=transpose(reshape(phi_A,[nx ny]));
AA=AA-AA(nx/2,1);
hold on
contourf(xm,ym,AA)
hs=streamslice(xm,ym,u(1:nx,1:ny)',v(1:nx,1:ny)');
set(hs,'color','w','linewidth',1)
axis square
title('$\phi$: Upwind Method','interpreter','latex','fontsize',12)
set(gcf, 'Color',[1,1,1]);
colorbar
hold off
nexttile
BB=transpose(reshape(phi_B,[nx ny]));
BB=BB- BB(nx/2,1);
hold on
contourf(xm,ym,BB)
hs=streamslice(xm,ym,u(1:nx,1:ny)',v(1:nx,1:ny)');
set(hs,'color','w','linewidth',1)
axis square
title('$\phi$: 2nd Order','interpreter','latex','fontsize',12)
set(gcf, 'Color',[1,1,1]);
colorbar
hold off
nexttile
CC=transpose(reshape(phi_C,[nx ny]));
CC=CC-CC(nx/2,1);
hold on
contourf(xm,ym,CC)
hs=streamslice(xm,ym,u(1:nx,1:ny)',v(1:nx,1:ny)');
set(hs,'color','w','linewidth',1)
axis square
title('$\phi$: QUICK Method','interpreter','latex','fontsize',12)
set(gcf, 'Color',[1,1,1]);
colorbar
hold off

% Centerline plots
figure(8)
tiledlayout(1,3)
nexttile
plot(AA(:,size(AA,2)/2), ym, 'Linewidth', 2)
title('Centerline Plot for Upwind')
xlabel('$\phi$', 'interpreter', 'latex')
ylabel('y')
nexttile
plot(BB(:,size(BB,2)/2), ym, 'Linewidth', 2)
title('Centerline Plot for 2nd Order Reconstruction')
xlabel('$\phi$', 'interpreter', 'latex')
ylabel('y')
nexttile
plot(CC(:,size(CC,2)/2), ym, 'Linewidth', 2)
title('Centerline Plot for QUICK')
xlabel('$\phi$', 'interpreter', 'latex')
ylabel('y')

%% Functions
function output = MeshGen_QUICK(nx, ny, u, v, dx, dy, alpha)
% A function to generate the diffusion operator of the linear system A*p = w
% INPUTS
%   nx: number of FV cells in x-direction
%   ny: number of FV cells in y-direction
%   u: x-component of velocity field
%   v: y-component of velocity field
%   alpha: diffusion constant

%OUTPUTS
%   A: Coefficient matrix for linear system

    % Generate diffusion operator
    tot = nx*ny;
    output = zeros(tot,tot);
    u_old = u;
    v_old = v;
    
    u_ctr = 1;
    v_ctr = 1;      %Counter to track placement in a row
    q = 1;          %Indexing term to read correct cell face values
    r = 1;
    
    %---- Generating index matrix ----
    indMat = zeros(nx+4, ny+4);
    indMat(3, [3:nx+2]) = [1:nx];
    
    %Generating initial indices
    for i = 4:ny+2  
        indMat(i, [3:nx+2]) = indMat(i-1, [3:nx+2]) + nx;
    end
    
    %Filling in periodic indices
        %Top and Bottom
        indMat([1:2], [3:nx+2]) = indMat([end-3:end-2], [3:nx+2]);
        indMat([end-1:end], [3:nx+2]) = indMat([3:4], [3:nx+2]);
        
        %Left and Right
        indMat([3:nx+2], [1:2]) = indMat([3:nx+2], [end-3:end-2]);
        indMat([3:nx+2], [end-1:end]) = indMat([3:nx+2], [3:4]);
        
    %Creating respective index matrices and decomposing
    Nindex = indMat([2:nx+1], [3:nx+2]);
    NNindex = indMat([1:nx], [3:nx+2]);
    Sindex = indMat([4:nx+3], [3:nx+2]);
    SSindex = indMat([5:nx+4], [3:nx+2]);
    Windex = indMat([3:nx+2], [2:nx+1]);
    WWindex = indMat([3:nx+2], [1:nx]);
    Eindex = indMat([3:nx+2], [4:nx+3]);
    EEindex = indMat([3:nx+2], [5:nx+4]);
    
    Nindex = reshape(transpose(Nindex), 1, []); Sindex = reshape(transpose(Sindex), 1, []);
    Eindex = reshape(transpose(Eindex), 1, []); Windex = reshape(transpose(Windex), 1, []);
    NNindex = reshape(transpose(NNindex), 1, []); SSindex = reshape(transpose(SSindex), 1, []);
    EEindex = reshape(transpose(EEindex), 1, []); WWindex = reshape(transpose(WWindex), 1, []);
    
    u_ctr = 1;
    v_ctr = 1;      %Counter to track placement in a row
    q = 1;          %Indexing term to read correct cell face values
    r = 1;
    
    %---- Generating coeff matrix, on a per-row basis ----
    for i = 1:tot
    % Tracks which row is being filled

        % Need to manually change the values of v due to how they were
        %   generated in the driver. 
        v = v_old(v_ctr, q:q+1);
        u = u_old(r:r+1, u_ctr);
        
        if mod(i, nx) == 0  
            v_ctr = 1;  %Resets u_ctr, v_ctr to 1 so we start at the top after we reach the end of a row of cells
            q = q + 1;  %at the same time, increments counter by one so we shift over which values we read from v_old, u_old
            
            r = 1;
            u_ctr = u_ctr + 1;
        else
            v_ctr = v_ctr + 1;
            r = r + 1;
        end
        
        %%%%%%%%%%
        % Put a break point below here to debug on a per-loop basis
        %%%%%%%%%%
    
        % Initializing per-loop terms
        EE = 0; WW = 0; NN = 0; SS = 0;
        temp = zeros(1,tot);
        
        P = 0; E = 0; W = 0; N = 0; S = 0;
        % Populating matrix coeffs with diffusion terms
        P = P + (2*alpha)*((1/(dx^2)) + (1/(dy^2)));
        E = E - alpha*(1/(dx^2));
        W = W - alpha*(1/(dx^2));
        N = N - alpha*(1/(dy^2));
        S = S - alpha*(1/(dy^2));
        
        % Populating matrix coeffs with remaining QUICK terms
        %  (+) ->          N
        %  |        i, j _____ i+1
        %  v            |     |
        %      WW    W  |  P  |  E    EE
        %               |_____|
        %            j+1
        %                  S
        
        if u(2) > 0       %u(i+1/2) > 0
            %Executes if flow is going P -> E
            P = P + (5/6)*(u(2)/dx);
            E = E + (2/6)*(u(2)/dx);
            W = W + (-1/6)*(u(2)/dx);

        else 
            %Executes if flow is going P <- E
            P = P + (2/6)*(u(2)/dx);
            E = E + (5/6)*(u(2)/dx);
            EE = EE + (-1/6)*(u(2)/dx);
            
        end
        
        if u(1) > 0         %u(i-1/2) > 0 
            %Executes if flow is going W -> P
            P = P - (2/6)*(u(1)/dx);
            W = W - (5/6)*(u(1)/dx);
            WW = WW - (-1/6)*(u(1)/dx);

        else
            %Executes if flow is going W <- P
            P = P - (5/6)*(u(1)/dx);
            W = W - (2/6)*(u(1)/dx);
            E = E - (-1/6)*(u(1)/dx);

        end
        
        if v(2) > 0       %v(j+1/2) > 0
            %Executes if flow is going P v S
            P = P + (5/6)*(v(2)/dy);
            N = N + (-1/6)*(v(2)/dy);
            S = S + (2/6)*(v(2)/dy);

        else
            %Executes if flow is going S ^ P
            P = P + (2/6)*(v(2)/dy);
            S = S + (5/6)*(v(2)/dy);
            SS = SS + (-1/6)*(v(2)/dy);
    
        end
        
        if v(1) > 0         %v(j-1/2) > 0 
            %Executes if flow is going N ^ P
            P = P - (2/6)*(v(1)/dy);
            N = N - (5/6)*(v(1)/dy);
            NN = NN - (-1/6)*(v(1)/dy);
            
        else
            %Executes if flow is going P v N
            P = P - (5/6)*(v(1)/dy);
            S = S - (-1/6)*(v(1)/dy);
            N = N - (2/6)*(v(1)/dy);

        end

        %---- Placing terms ----
        temp(i) = P;
        temp(Nindex(i)) = N;
        temp(Eindex(i)) = E;
        temp(Windex(i)) = W;
        temp(Sindex(i)) = S;
        temp(NNindex(i)) = NN;
        temp(EEindex(i)) = EE;
        temp(WWindex(i)) = WW;
        temp(SSindex(i)) = SS;
        
        output(i,:) = temp;

    end
end

function output = MeshGen_Upwind(nx, ny, u, v, dx, dy, alpha)
% A function to generate the diffusion operator of the linear system A*p = w
% INPUTS
%   nx: number of FV cells in x-direction
%   ny: number of FV cells in y-direction
%   u: x-component of velocity field
%   v: y-component of velocity field
%   alpha: diffusion constant

%OUTPUTS
%   A: Coefficient matrix for linear system

    % Generate diffusion operator
    tot = nx*ny;
    output = zeros(tot,tot);
    u_old = u;
    v_old = v;
    
    u_ctr = 1;
    v_ctr = 1;      %Counter to track placement in a row
    q = 1;          %Indexing term to read correct cell face values
    r = 1;
    lc_ctr = 1;
    
    % Generating coeff matrix, on a per-row basis
    for i = 1:tot
    % Tracks which row is being filled
    
        % Need to manually change the values of v due to how they were
        %   generated in the driver. 
        v = v_old(v_ctr, q:q+1);
        u = u_old(r:r+1, u_ctr);
        
        if mod(i, nx) == 0  
            v_ctr = 1;  %Resets u_ctr, v_ctr to 1 so we start at the top after we reach the end of a row of cells
            q = q + 1;  %at the same time, increments counter by one so we shift over which values we read from v_old, u_old
            
            r = 1;
            u_ctr = u_ctr + 1;
        else
            v_ctr = v_ctr + 1;
            r = r + 1;
        end
        
        P = 0; E = 0; W = 0; N = 0; S = 0;
        % Populating matrix coeffs with diffusion terms
        P = P + (2*alpha)*((1/(dx^2)) + (1/(dy^2)));
        E = E - alpha*(1/(dx^2));
        W = W - alpha*(1/(dx^2));
        N = N - alpha*(1/(dy^2));
        S = S - alpha*(1/(dy^2));
        
        % Populating matrix coeffs with remaining Upwind terms
        %  (+) ->     N
        %  |      _|_____|_
        %  v       |     |
        %       W  |  P  |  E
        %         _|_____|_
        %          |     |
        %             S
        
        if u(2) > 0         %u(i+1/2) > 0
            %Executes if flow is going P -> E, Upwind is P
            P = P + (u(2)/dx);
        else 
            %Executes if flow is going P <- E, Upwind is E
            E = E + (u(2)/dx);
        end
        
        if u(1) > 0         %u(i-1/2) > 0 
            %Executes if flow is going W -> P, Upwind is W
            W = W + (-u(1)/dx);
        else
            %Executes if flow is going W <- P, Upwind is P
            P = P + (-u(1)/dx);
        end
        
        if v(2) > 0         %v(i+1/2) > 0
            %Executes if flow is going P v S, Upwind is P
            P = P + (v(2)/dy);
        else
            %Executes if flow is going P ^ S, Upwind is S
            S = S + (v(2)/dy);     
        end
        
        if v(1) > 0         %v(i-1/2) > 0 
            %Executes if flow is going N v P, Upwind is N
            N = N + (-v(1)/dy);
        else
            %Executes if flow is going N ^ P, Upwind is P
            P = P + (-v(1)/dy);
        end
        
        % Collecting data and calling coeff placement algo
        coeffs = [P, E, W, N, S];
        dat = [i, tot, nx, ny, lc_ctr];
        [output(i,:), lc_ctr] = periodicMesh(coeffs, dat, lc_ctr);

    end
end

function output = MeshGen_2ndO(nx, ny, u, v, dx, dy, alpha)
% A function to generate the diffusion operator of the linear system A*p = w
% INPUTS
%   nx: number of FV cells in x-direction
%   ny: number of FV cells in y-direction
%   u: x-component of velocity field
%   v: y-component of velocity field
%   alpha: diffusion constant

%OUTPUTS
%   A: Coefficient matrix for linear system

    % Generate diffusion operator
    tot = nx*ny;
    output = zeros(tot,tot);
    u_old = u;
    v_old = v;
    u = u(:);
    
    u_ctr = 1;
    v_ctr = 1;      %Counter to track placement in a row
    q = 1;          %Indexing term to read correct cell face values
    r = 1;
    lc_ctr = 1;
    
    % Generating coeff matrix, on a per-row basis
    for i = 1:tot
    % Tracks which row is being filled

        % Need to manually change the values of v due to how they were
        %   generated in the driver. 
        v = v_old(v_ctr, q:q+1);
        u = u_old(r:r+1, u_ctr);
        
        if mod(i, nx) == 0  
            v_ctr = 1;  %Resets u_ctr, v_ctr to 1 so we start at the top after we reach the end of a row of cells
            q = q + 1;  %at the same time, increments counter by one so we shift over which values we read from v_old, u_old
            
            r = 1;
            u_ctr = u_ctr + 1;
        else
            v_ctr = v_ctr + 1;
            r = r + 1;
        end
        
        P = 0; E = 0; W = 0; N = 0; S = 0;
        % Populating matrix coeffs with diffusion terms
        P = P + (2*alpha)*((1/(dx^2)) + (1/(dy^2)));
        E = E - alpha*(1/(dx^2));
        W = W - alpha*(1/(dx^2));
        N = N - alpha*(1/(dy^2));
        S = S - alpha*(1/(dy^2));
        
        % Generating matrix coeffs
        P = P + ((u(2) - u(1))/(2*dx)) + ((v(2) - v(1))/(2*dy));
        E = E + (u(2)/(2*dx));
        W = W - (u(1)/(2*dx));
        N = N - (v(1)/(2*dy));
        S = S + (v(2)/(2*dy));

        % Collecting data and calling coeff placement algo
        coeffs = [P, E, W, N, S];
        dat = [i, tot, nx, ny, lc_ctr];
        [output(i,:), lc_ctr] = periodicMesh(coeffs, dat, lc_ctr);
        
    end    
end

function [temp, lc_ctr] = periodicMesh(coeffs, dat, lc_ctr)
%This function populates FV coefficients based on 2D periodic boundary conditions
% INPUTS
%   coeffs: vector of coefficient data for the stencil
%   dat: necessary loop data for each rock of "A" matrix being populated
% OUTPUTS
%   temp: populated row vector of "i"th row of "A" coeff matrix

        %Unwrapping Data
        P = coeffs(1); E = coeffs(2); W = coeffs(3); N = coeffs(4); S = coeffs(5);
        i = dat(1); n = dat(2); nx = dat(3); ny = dat(4); lc_ctr = dat(5);

        %Setting up temporary rows for computation
        temp = zeros(1,n);
        
        %Populating temp row
        temp(i) = P; 
        
        
        switch true
            case i == 1                 %Upper left corner
                temp(i+1) = E; temp(nx) = W; temp(i+nx) = S; temp(n-(nx-1)) = N;
                lc_ctr = lc_ctr + nx;
            
            case i == n                 %Lower right corner
                temp(n-(nx-1)) = E; temp(i-1) = W; temp(nx) = S; temp(n-nx) = N;
                
            case i == n-(nx-1)          %Lower left corner
                temp(i+1) = E; temp(n) = W; temp(1) = S; temp(i-nx) = N;
                
            case i == nx                %Upper right corner
                temp(1) = E; temp(i-1) = W; temp(i+nx) = S; temp(n) = N;
                
            case i > 1 && i < nx        %Top row BC
                temp(i+1) = E; temp(i-1) = W; temp(i+nx) = S; temp(n-nx+i) = N;
                
            case i > (n-nx+1) && i < n  %Bottom row BC
                temp(i+1) = E; temp(i-1) = W; temp(i-(n-nx)) = S; temp(i-nx) = N;
                
            case i == lc_ctr            %Left column BC
                temp(i+1) = E; temp(i+(nx-1)) = W; temp(i+nx) = S; temp(i-nx) = N;
                lc_ctr = lc_ctr + nx;
                
            case mod(i, nx) == 0        %Right Column BC
                temp(i-(nx-1)) = E; temp(i-1) = W; temp(i+nx) = S; temp(i-nx) = N;
            
            otherwise                   %Interior nodes
                temp(i+1) = E; temp(i-1) = W; temp(i+nx) = S; temp(i-nx) = N;
               
        end
end

