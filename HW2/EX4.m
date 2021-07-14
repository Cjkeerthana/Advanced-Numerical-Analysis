%%
%Adavnced Numerical Analysis - HW 2 - Keerthana C J
%%
%Exercise 4
A = gallery('wathen',100,100);
n=size(A,1);
b=rand(n,1);
tol = 10^(-8);
maxit=500;

[x,flag1,relres, iter1,resvec1]=pcg(A,b,tol,maxit);  %No preconditioner

L=ichol(A);
[x,flag2,relres, iter2,resvec2]=pcg(A,b,tol,maxit,L,L');%IC(0) preconditioner

opts.type='ict';
opts.droptol=10^-2;
L=ichol(A,opts);
[x,flag3,relres, iter3,resvec3]=pcg(A,b,tol,maxit,L,L');%ICT with drop tol 10^-2

opts.type='ict';
opts.droptol=10^-3;
L=ichol(A,opts);
[x,flag4,relres, iter4,resvec4]=pcg(A,b,tol,maxit,L,L');%ICT with drop tol 10^-3

L = diag(diag(A));
[x,flag5,relres, iter5,resvec5]=pcg(A,b,tol,maxit,L,[]);%Jacobi preconditioner

%plots
figure(2)
semilogy(0:length(resvec1)-1,resvec1/norm(b),'-o')
hold on
semilogy(0:length(resvec2)-1,resvec2/norm(b),'-o')
semilogy(0:length(resvec3)-1,resvec3/norm(b),'-o')
semilogy(0:length(resvec4)-1,resvec4/norm(b),'-o')
semilogy(0:length(resvec5)-1,resvec5/norm(b),'-o')
yline(tol,'r--');
legend('nopreconditioner','IC0','IC with 10^-2','IC with 10^-3','Jacobi')
xlabel('Iteration number')
ylabel('Residual')
clear;