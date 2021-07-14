%%
%Adavnced Numerical Analysis - HW 2 - Keerthana C J
%%
A = load('apache1.mtx','-ascii'); % Apache matrix
%A = load('bcsstk13.mtx','-ascii'); % A structural sparse matrix;
A = spconvert(A);
A = A+A'-diag(diag(A));
n= size(A,1);
b = A*ones(n,1);
L = ichol(A); % not to be used for the structural spare matrix
tol = 10^(-8);
maxit=2000;

[x,flag1,relres, iter1,resvec1]=pcg(A,b,tol,maxit);
[x,flag,relres, iter,resvec]=pcg(A,b,tol,maxit, L, L');

figure(1)
semilogy(0:length(resvec1)-1,resvec1/norm(b),'-o')
hold on
semilogy(0:length(resvec)-1,resvec/norm(b),'-o')
yline(tol,'r--');
legend('no preconditioner','IC(0)')
xlabel('Iteration number')
ylabel('Residual')
title('Apache matrix without scaling')

[i,j,a]=find(A);
maxval = max(abs(a)); % maximum entry of A 

%scaling the apache system by dividing with the maximum value of A
A1 = A./maxval;
b1=b./maxval;
L1 = ichol(A1);

% for the structural sparse matrix
%alpha = max(sum(abs(A),2)./diag(A))-2;
%alpha = 0.1;
%L = ichol(A, struct('type','ict','droptol',1e-3,'diagcomp',alpha));

[x,flag2,relres, iter2,resvec2]=pcg(A1,b1,tol,maxit, L1 ,L1');
[x,flag3,relres, iter3,resvec3]=pcg(A1,b1,tol,maxit);

%plots
figure(2)
semilogy(0:length(resvec3)-1,resvec3/norm(b1),'-o')
hold on
semilogy(0:length(resvec2)-1,resvec2/norm(b1),'-o')
yline(tol,'r--');
legend('no preconditioner', 'IC(0)')
xlabel('Iteration number')
ylabel('Residual')
title('Apache matrix after scaling')
%title('A strcutural spare matrix')
clear;