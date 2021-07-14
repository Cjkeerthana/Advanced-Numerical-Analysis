%%
%Exercise 2:
%Finding an eigen vector pair(x,lambda) given their initial guesses using
%Netwon's method with the constrain ||x||=1
%%
nx = 50;
A = delsq(numgrid('S',nx+2));
n = size(A,1);
[u0,lamda0] = eigs(A,1,'sm'); %The real eigen value
%Perturbing the eigen value and vector
u = u0 + ones(n,1)*1e-2; 
lamda = lamda0 + 1e-2;
x0 = [u;lamda];

%The non-linear function vector F(x) = (Ax-lambda*x) + x'x -1
ffun = @(x)([(A*x(1:size(x,1)-1)-x(size(x,1))*x(1:size(x,1)-1));x(1:size(x,1)-1)'*x(1:size(x,1)-1)-1]);
tol = 1e-12;
itmax=20;

tic
[xstar,iter,resvec,stepdiff] = newton(x0,@(x)ffun(x),@(x)jfun(x,ffun),tol,itmax,-1);
time = toc;

fprintf("The Calculated eigen value is %f\n", xstar(n+1));
fprintf("The precision is %f\n",abs(xstar(n+1)-lamda));

figure(1)
semilogy(0:length(resvec)-1,resvec,'-o')
hold on
yline(tol,'r--');
legend('GMRES')
xlabel('Iteration number')
ylabel('Non-linear Residual')
title('Newtons Method');
clear;
function y = jfun(x,ffun)
n = size(x,1);
y = zeros(n);
h = diag(1e-6*ones(size(x)));
for i = 1:n
    y(:,i) = (ffun(x+h(:,i))-ffun(x-h(:,i)))/(2*h(i,i));
end
end