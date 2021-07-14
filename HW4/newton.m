function [xstar,iter,resvec,stepdiff] = newton(x0,F,Jac,tol,itmax,lsol)
% Newton's method for solving non-linear equations
% [xstar,iter,resvec,stepdiff] = newton(x0,F,Jac,tol,itmax,lsol,varargin)
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
% tol: tolerance for the exit test (||F(x_k)||/||F(x0)||)
% itmax: maximum number of outer iterations
% lsol: For solving each non-linear iteration
%       1: Use matlab operator \ ; 0: GMRES ; 1: LU Factorization
    err = norm(F(x0));
    x = x0;
    resvec = [];
    stepdiff = [];
    xstar = [x0];
    iter = 0;
    tau = 0.5;
    while err > tol*norm(F(x0)) && iter < itmax
        if lsol == 1
            s = -Jac(x)\F(x);
        elseif lsol == 0
            J = -Jac(x);
            if(size(J,1)<10)  %If the matrix dimensions are too small
                s = gmres(J,F(x),[],tol);
            else
                J = sparse(J);
                [L,U] = ilu(J);
                restart = 50;
                gm_tol = 1e-8;
                gm_itmax = 20;
                s = gmres(J,F(x),restart,gm_tol,gm_itmax,L,U);
            end
        elseif lsol == -1
            [L,U,P] = lu(-Jac(x));
            y = L\(P*F(x));
            s = U\y;
        end
        alpha = 1.0; %For backtracking
        flag = 0;
        while flag == 0 
            xt = x + alpha*s;
            if norm(F(xt)) < norm(F(x))
                x = xt;
                flag = 1;
            else
                alpha = tau*alpha;
                flag = 0;
            end
        end
        xstar = [xstar,x];
        resvec = [resvec;norm(F(x))];
        stepdiff = [stepdiff;norm(alpha*s)];
        err = norm(F(x)); 
        iter = iter + 1;
    end
end