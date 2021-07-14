%%
%Advanced Numerical Analysis Homework 1 - Keerthana C J
%%
%Exercise 1 & 2
k = 1; %for the plots
for nx = 20:20:80 %nx ==> number of divisions in each side of the square domain
    %discretization using FDM 
    A = delsq(numgrid('S',nx+2));        
    n = size(A,1); 
    b = (1/(nx+1)^2) * ones(n,1);
    x0 = zeros(n,1);
    
    tol = 10^(-8); %tolerance
    nmax = 10^5; %maximum number of iterations
    
    [xJ, it_J, err_J] = jor(A,b,x0,nmax,tol,1); %Jacobi (omega = 1)
    [xGS, it_GS, err_GS] = sor(A,b,x0,nmax,tol,1); % Guass Seidel (omega = 1)
    
    %plots
    itJ = 1:it_J;
    itS = 1:it_GS;
    subplot(2,2,k)
    plot(itJ,err_J,'b',itS,err_GS,'r'), xlabel('Iterations'), ylabel('Error')
    title(nx)
    legend('Jacobi','Gauss Siedel')
    k = k+1;
end
clear;

%Exercise 3
% FDM discretization of the domain
nx = 30;
A = delsq(numgrid('S',nx+2));         
n = size(A,1);
b = (1/(nx+1)^2) * ones(n,1);
x0 = zeros(n,1);
tol = 10^(-8);
nmax = 10^5;

D = diag(diag(A));
I = eye(nx^2,nx^2);
L = tril(A,-1);
U = triu(A,1);
H_J = I - D\A; %Jacobi Iteration Matrix
H_GS = -(D+L)\U; %GS iteration Matrix

%comparing the different methods
[xJ, it_J, err_J] = jor(A,b,x0,nmax,tol,1); % Jacobi (omega = 1)
[xGS, it_GS, err_GS] = sor(A,b,x0,nmax,tol,1); % Guass Siedel (omega = 1)

h = 1/(nx+1); % size of each interval 
rho_HJ = 2*(cos(h*pi/2))^2 - 1; % spectral radius of the Jacobi iteration matrix
omega_opt = 2/(1+sqrt(1-(rho_HJ)^2)); % optimum value for Omega

H_SOR_opt = (D+omega_opt*L)\((1-omega_opt)*D-omega_opt*U); %SOR optimum iteration matrix

%Relaxation Methods with optimum value of Omega, 
%omega for SOR in the interval of (0,2)
%omega for JOR in the interval of (0,1] as Jacobi is convergent
[xJOR, it_JOR, err_JOR] = jor(A,b,x0,nmax,tol,0.5); % Jacobi over-relaxation
[xS, it_S, err_S] = sor(A,b,x0,nmax,tol,omega_opt); %Siedel Over-relaxation

rho_HGS = eigs(H_GS,1); % spectral radius of GS iteration matrix
rho_HSOR = eigs(H_SOR_opt,1); % spectral raidus of SOR iteration matrix

p = -log10(tol);

%Estimated number of iterations
k_J = p/(-log10(rho_HJ));  % Jacobi
k_GS = p/(-log10(rho_HGS)); % GS
k_SOR = p/(-log10(rho_HSOR)); %SOR

itJ = 1:it_J;
itJOR = 1:it_JOR;
itGS = 1:it_GS;
itS = 1:it_S;
figure()
plot(itJ,err_J,'b',itGS,err_GS,'r',itJOR,err_JOR,'y',itS,err_S,'g'), xlabel('Iterations'), ylabel('Error')
title(nx)
legend('Jacobi','GS','JOR','SOR')
it_J
it_GS
it_JOR
it_S
k_J
k_GS
k_SOR
   