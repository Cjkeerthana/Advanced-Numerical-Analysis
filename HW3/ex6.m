%%
%Exercise 6 - A saddle point problem has been solved using the MINRES
%function. 
%%
%Input Matrices of the Saddle point problem
load 'Q.dat';
load 'A.dat';
Q = spconvert(Q);
Q = Q+Q'-diag(diag(Q));
A = spconvert(A);
[m,n] = size(A);

x1 = 10^(-6)*rand(n,1);
x2 = rand(m,1);
exact = [x1;x2];
b1 = Q*x1+A'*x2;
b2 = A*x1;
b = [b1;b2];

maxit = 550;
tol = 10^(-12);

%Block Preconditioners
E = diag(diag(Q));
S = A*(E\A');

%Cholesky factorization of the Schur Complement
p = amd(S);
S0 = S(p,p);
U = chol(S0);

%Intial guess vector is taken as M-1.b
x01 = E\b1;
x02 = applyschur(b2,p,U);
x0 = [x01;x02];

%MINRES method
tic
[x,flag,relres,iter,resvec] = minres(@(x)[Q*x(1:n)+A'*x(n+1:n+m);A*x(1:n)],b,tol,maxit,@(x)[E\x(1:n);applyschur(x(n+1:n+m),p,U)],[],x0);
time = toc;

fprintf("The time taken for MINRES is %f s \n", time);
fprintf("No. of iterations = %d \n",iter);

figure(1)
semilogy(0:length(resvec)-1,resvec/norm(b),'-o')
hold on
yline(tol,'r--');
xlabel('Iteration number')
ylabel('Relative Residual')
title('MINRES - Saddle point problem')
legend('minres','tol')
clear;

%Function to solve a linear system with Schur Complement
function w = applyschur(v,p,U)
n = size(U,1);
z = v(p);
u = U\(U'\z);
pinv(p) = 1:n;
w = u(pinv);
end
