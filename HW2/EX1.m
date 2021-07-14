%%
%Adavnced Numerical Analysis - HW 2 - Keerthana C J
%%
%Exercise-1
A=delsq(numgrid('S',102));
L=ichol(A);
n=size(A,1);
b = A * ones(n,1);
tol = 10^(-8);
maxit=500;
M = L*L'; %Preconditioner is LL^T

[x,flag2,relres, iter2,resvec2]=pcg(A,b,tol,maxit,L,L'); %matlab function
[x,resvec,iter] = mypcg(A,b,tol,maxit,M); %implementation of pcg with fuction mypcg

%plots
figure(1)
semilogy(0:length(resvec)-1,resvec/norm(b),'-o')
hold on
semilogy(0:length(resvec2)-1,resvec2/norm(b),'-o')
yline(tol,'r--');
legend('mypcg','matlab')
xlabel('Iteration number')
ylabel('Residual')
clear;