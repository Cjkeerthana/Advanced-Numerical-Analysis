%%
%Exercise 5:
%Intersection two quadratic curves in 2D is solved using the Broyden's method 
%%
x0 = [0;1];
ffun = @(x)[x(1)*x(1)+x(2)*x(2)-4; x(1)*x(2)-1];
itmax = 50;
tol = 1e-8;

tic
[xstar,iter,resvec,stepdiff] = quasinewton(x0,@(x)ffun(x),@(x)jfun(x,ffun),tol,itmax);
time = toc;
fprintf("Time for solution using Broyden's method is %f\n",time);
clear;
function y = jfun(x,ffun)
n = size(x,1);
y = zeros(n);
h = diag(1e-6*ones(size(x)));
for i = 1:n
    y(:,i) = (ffun(x+h(:,i))-ffun(x-h(:,i)))/(2*h(i,i));
end
end