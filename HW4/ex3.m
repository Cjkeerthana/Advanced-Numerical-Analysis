%%
%Exercise 2:
%Solving a non-linear solution involving a Laplacian matrix
%%
nx = 48;
A = delsq(numgrid('S',nx+2));
n = size(A,1);
x0 = 1000* sin([1:n]');
tol = 1e-12;
itmax = 100;
ffun = @(x)(A*x-0.1*sin(x)-5); %The non-linear system of eqns F(x)

%For solving the each non-linear iteration
%Using the \ operator
tic
[xstar1,iter1,resvec1,stepdiff1] = newton(x0,@(x)ffun(x),@(x)jfun(x,ffun),tol,itmax,1);
time1 = toc;

%Using GMRES
tic
[xstar0,iter0,resvec0,stepdiff0] = newton(x0,@(x)ffun(x),@(x)jfun(x,ffun),tol,itmax,0);
time0 = toc;

tic
[xstar2,iter2,resvec2,stepdiff2] = newton(x0,@(x)ffun(x),@(x)jfun(x,ffun),tol,itmax,-1);
time2 = toc;

fprintf("\n Time for solution using direct method = %f",time1);
fprintf("\n Time for solution using GMRES method = %f",time0);
fprintf("\n Time for solution using LU factorization = %f\n",time2);

figure(1)
semilogy(0:length(resvec1)-1,resvec1,'-o')
hold on
semilogy(0:length(resvec2)-1,resvec2,'-o')
semilogy(0:length(resvec0)-1,resvec0,'-o')
yline(tol,'r--');
legend('MatlabOperator','LU','GMRES')
xlabel('Iteration number')
ylabel('Non-linear Residual')
title('Newtons Method');
clear;
function y = jfun(x,ffun)
%Function for calculating the Jacobian using a numerical scheme
n = size(x,1);
y = zeros(n);
h = diag(1e-6*ones(size(x)));
for i = 1:n
    y(:,i) = (ffun(x+h(:,i))-ffun(x-h(:,i)))/(2*h(i,i));
end
end