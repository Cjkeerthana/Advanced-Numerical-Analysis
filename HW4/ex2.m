%%
%Exercise 1:
%Intersection two quadratic curves in 2D is solved using the Newton's method 
%%
x0 = [0;1];
ffun = @(x)[x(1)*x(1)+x(2)*x(2)-4; x(1)*x(2)-1]; %The non-linear system of eqns F(x)
itmax = 50;
tol = 1e-8;
times = [];

%For solving the each non-linear iteration
%Using the \ operator
tic
[xstar1,iter1,resvec1,stepdiff1] = newton(x0,@(x)ffun(x),@(x)jfun(x,ffun),tol,itmax,1); 
times = [times;toc];

%Using GMRES
tic
[xstar2,iter2,resvec2,stepdiff2] = newton(x0,@(x)ffun(x),@(x)jfun(x,ffun),tol,itmax,0); 
times = [times;toc];

%Using LU Factorization
tic
[xstar3,iter3,resvec3,stepdiff3] = newton(x0,@(x)ffun(x),@(x)jfun(x,ffun),tol,itmax,-1);
times = [times;toc];

figure(1)
semilogy(0:length(resvec1)-1,resvec1,'-o')
hold on
semilogy(0:length(resvec2)-1,resvec2,'-o')
semilogy(0:length(resvec3)-1,resvec3,'-o')
yline(tol,'r--');
legend('MatlabOperator','GMRES','LU')
xlabel('Iteration number')
ylabel('Non-linear Residual')
title('Newtons Method');
Iter_No = ['1';'2';'3';'4'];
x1 = xstar1(1,2:5)';
x2 = xstar1(2,2:5)';
%Table displaying the values at each iteration
T = table(Iter_No,x1,x2,resvec1,stepdiff1);
display(T);
fprintf("Time for solution using newton's method is %f",times(1));
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