%Generalization of Jacobi ==> Jacobi Over relaxation Method
%Args ==> A = input matrix; b = right hand vector; x0 = intial guess vector;
          %nmax = max number of intial iterations; tol=tolerance
          %omega=relaxation parameter
function[x,iter, err_J] = jor(A,b,x0,nmax,tol,omega)
[n,m] = size(A);
if n ~= m, error('Only Square Systems'); end
iter = 0;
r = b - A*x0; %residual vector
%r0 = norm(r);
err = norm(r); %initial error is taken as the norm of the residual
x = x0; 
err_J = [];
D = diag(diag(A)); 
while err > tol && iter < nmax
    iter = iter + 1;
    temp = D\r; % the vector form of the Jacobi Over relaxation method is given as 
    xnew = x + omega*temp; % xnew = x + omega*inverse(D)*r
    
    err = norm(xnew-x); %error is measured as the distance between 
                        %solutions of the successive iterations
    x = xnew; %update x
    r = b-A*x; %update residual
    %err = norm(r)/r0; %another measure of error
    err_J = [err_J;err]; %error at the end of each iteration is stored
end
return
