%Generalization of Guass Seidel ==> Seidel Over relaxation Method
%Args ==> A = input matrix; b = right hand vector; x0 = intial guess vector;
          %nmax = max number of intial iterations; tol=tolerance
          %omega=relaxation parameter
function[x,iter, err_S] = sor(A,b,x0,nmax,tol,omega)
[n,m] = size(A);
if n ~= m, error('Only Square Systems'); end %checking the matrix dimensions
iter = 0;
r = b - A*x0; %residual vector
%r0 = norm(r);
err = norm(r); %intial error is taken as the norm of the residual
x = x0; 
xold = x0;
err_S = [];

D = diag(diag(A));
L = tril(A,-1);
%E = -L;

while err > tol && iter < nmax
    iter = iter + 1;
    temp = (1/omega) * D + L; %the vector form of the SOR method is given as 
    x = xold + temp\r; % x = xold + inverse((1/omega)*D-E)*r
    
    err = norm(x-xold); %error estimation
    
    xold = x; %update xold
    r = b-A*x; %update residual
    %err = norm(r)/r0; %another measure of error
    err_S = [err_S; err]; %error at the end of each iteration is stored
end
return
