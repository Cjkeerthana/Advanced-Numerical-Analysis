function[x,iter,resvec,flag]=mygmres(A,b,x0,maxit,tol)
% GMRES - Krylov Subspace method
% [x,iter,resvec,flag]=mygmres(A,b,x0,maxit,tol)
% OUTPUT parameters: 
% x: solution vector
% iter: number of iterations
% resvec: vector containing the residual norm : |beta*Q(1,k+1)|/||b||
% flag: values = 1->gmres converged within the max iterations; 0->Lucky breakdown of
% gmres; -1 -> gmres did not converge within the max no of iterations
% INPUT parameters: 
% A,b: coefficient matrix and right hand side
% x0: initial guess 
% tol: tolerance for the exit test 
% maxit: maximum number of iterations

r0 = b-A*x0;
iter=0;
rho = norm(r0);
beta = rho;
v1 = r0/beta;
[n,n] = size(A);
m = maxit;
V = zeros(n,m+1);
H = zeros(m+1,m);
V(:,1)=v1;
resvec=[];
while rho > tol*norm(b) && iter<maxit 
    iter = iter+1;
    V(:,iter+1)=A*V(:,iter);
    for j = 1:iter
        H(j,iter) = (V(:,iter+1)')*V(:,j);   %%the modified gram -schmidt to generate orthonormal basis
        V(:,iter+1) = V(:,iter+1)-H(j,iter)*V(:,j);
    end
    H(iter+1,iter) = norm(V(:,iter+1));
    if (H(iter+1,iter) > 10^(-15))
        V(:,iter+1) = V(:,iter+1)/H(iter+1,iter);
        [Q,R]=qr(H); % QR factorization of the Hessenberg matrix
        rho = abs(Q(1,iter+1)*beta);
        resvec = [resvec;rho];
        flag = 1;
    else
        flag = 0; %lucky breakdown of gmres i.e. exact soultion has been found
        break;
    end
end
if flag == 1
    [Q,R] = qr(H(1:iter+1,1:iter)); %QR factorization of the final Hessenberg matrix
    Rhat = R(1:iter,1:iter);
    betae1 = beta*Q(1,1:iter)';
    y = Rhat\betae1;
    x = x0+V(1:n,1:iter)*y; %solution vector
elseif flag == 0
    betae1 = beta*eye(iter,1);
    y = H(1:iter)\betae1;
    x = x0+V(1:n,1:iter)*y;
end
if iter >= maxit %GMRES did not converge within the maximum number of iterations provided.
    flag = -1;
end
return