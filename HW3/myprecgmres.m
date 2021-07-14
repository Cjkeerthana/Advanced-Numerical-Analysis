function[x,iter,resvec,flag]=myprecgmres(A,b,tol,maxit,x0,ptype,L,U)
% GMRES with preconditioning - Krylov Subspace method
% [x,iter,resvec,flag]=mygmres(A,b,tol,maxit,x0,ptype,L,U)
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
% ptype: Type of precondtioning
% L,U: ILU preconditionder matrices

U_inv = inv(U); %inverting the L and U matrices as they are triangular
L_inv = inv(L);

if ptype == 'L'
    r0 = U_inv*(L_inv*(b-A*x0));
    check_criteria = norm(U_inv*(L_inv*b));
elseif ptype == 'R'
    r0 = b-A*x0;
    check_criteria = norm(b);
elseif ptype == 'S'
    r0 = L_inv*(b-A*x0);
    check_criteria = norm(L_inv*(b));
end
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

while rho > tol*check_criteria && iter<maxit 
    iter = iter+1;
    if ptype == 'L'
        z = A*V(:,iter);
        w = L\z;
        V(:,iter+1) = U\w;
    elseif ptype == 'R'
        z = L\V(:,iter);
        w = U\z;
        V(:,iter+1) = A*w;
    elseif ptype == 'S'
        z = U\V(:,iter);
        w = A*z;
        V(:,iter+1) = L\w;
    end
    for j = 1:iter
        H(j,iter) = (V(:,iter+1)')*V(:,j);
        V(:,iter+1) = V(:,iter+1)-H(j,iter)*V(:,j);
    end
    H(iter+1,iter) = norm(V(:,iter+1));
    if (H(iter+1,iter) > 10^(-15))
        V(:,iter+1) = V(:,iter+1)/H(iter+1,iter);
        [Q,R]=qr(H);
        rho = abs(Q(1,iter+1)*beta);
        resvec = [resvec;rho];
        flag = 1;
    else
        flag = 0;
        break;
    end
end
if flag == 1
    [Q,R] = qr(H(1:iter+1,1:iter));
    Rhat = R(1:iter,1:iter);
    betae1 = beta*Q(1,1:iter)';
    y = Rhat\betae1;
    if ptype == 'L'
        x = x0+V(1:n,1:iter)*y;
    elseif ptype == 'R'
        x = x0 + U_inv*(L_inv*(V(1:n,1:iter)*y));
    elseif ptype == 'S'
        x = x0 + U_inv*(V(1:n,1:iter)*y);
    end
elseif flag == 0
    betae1 = beta*eye(iter,1);
    y = H(1:iter)\betae1;
    if ptype == 'L'
        x = x0+V(1:n,1:iter)*y;
    elseif ptype == 'R'
        x = x0 + U_inv*(L_inv*(V(1:n,1:iter)*y));
    elseif ptype == 'S'
        x = x0 + U_inv*(V(1:n,1:iter)*y);
    end
end
if iter == maxit
    flag = -1;
end
return