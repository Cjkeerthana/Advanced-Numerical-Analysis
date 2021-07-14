function [xstar,iter,resvec,stepdiff] = quasinewton(x0,F,Jac,tol,itmax)
% [xstar,iter,resvec,stepdiff] = quasinewton(x0,F,Jac,tol,itmax)
% OUTPUT parameters: 
% xstar: solution vector for all iteration steps
% iter: number of outer iterations
% resvec: vector containing the non-linear residual norm : ||F(x_k)||
% stepdiff: norm of each step: ||alpha*s||
% INPUT parameters:  
% x0: initial guess 
% F: The non-linear system of equations F(x) given as a function handle
% returning a vector
% Jac: The Jacobian of the system F'(x) which returns a matrix of n by n
%      The initial B0 matrix = F'(x)
% tol: tolerance for the exit test (||F(x_k)||/||F(x0)||)
% itmax: maximum number of outer iterations

%LU Factorization of the initial B0 matrix
[L,U,P] = lu(Jac(x0));
y = L\(P*(-F(x0)));
s0 = U\y; 
err = norm(F(x0));
iter = 0;
s = [s0];
k = 0;
xstar = [x0];
resvec = [err];
stepdiff = [norm(s0)];
x = x0;
while err > tol && iter < itmax
    x = x + s(:,k+1);
    y = L\(P*(-F(x)));
    z = U\y;
    k = k+1;
    if k > 1
        for j = 1:(k-1)
            temp = s(:,j)'*z;
            z = z + s(:,j+1)*temp/(norm(s(:,j)))^2;
        end    
    end
    temp = s(:,k)'*z/norm(s(:,k))^2;
    sk_1 = z/(1-temp);
    s = [s,sk_1];
    err = norm(F(x))/norm(F(x0));
    xstar = [xstar,x];
    resvec = [resvec;norm(F(x))];
    stepdiff = [stepdiff;norm(s(:,k+1))];
    iter = iter + 1;
end
end