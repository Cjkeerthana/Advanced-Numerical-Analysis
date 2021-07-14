function [x, resvec, relres, iter] = twogrid(A, b, x0, tol, maxit, smooth_it)
% Two Grid Cycle: Multi Grid Method
% [x, resvec, relres, iter]=twogrid(A,b,x0,tol,maxit)
% OUTPUT parameters: 
% x: solution vector
% resvec: the residual at the end of each iteration
% relres: the relative residual at convergence
% iter: number of iterations
% INPUT parameters: 
% A,b: coefficient matrix and right hand side
% x0: initial guess 
% maxit: maximum number of iterations
% tol: tolerance for the exit test 
% smooth_it: the number of GS iterations to be applied for smoothing

%tic
[P,R] = prolongation(A); % calculation of prolongation &
%ptime = toc;            % restriction matrices

Ac = R * A * P;          % Galerkin coarse grid operator

x = x0;
err = norm(b - A*x);
normb = norm(b);
resvec = [];
iter = 0;
relres = err/normb;
pcg_maxit = 5000;
M = tril(A);        % matrix M = L =>Guass Siedel Smoothing

while err > normb*tol && iter < maxit
    %tic
    x = smooth(A,b,x,M,smooth_it);      % Pre-smoothing
    %stime = toc;
    
    rc = R *(b-A*x);        % Calculation of restricted residual
    
    %tic
    opts.droptol = 10^(-3);
    Lc = ichol(Ac,opts);  %IC factorization with droptol 10(-3)
    [ec,~,~,~,~] = pcg(Ac,rc,tol,pcg_maxit,Lc,Lc'); % solve for restricted error
    %sltime = toc;
    
    x = x + P*ec;           % Prolong the error and 
                            % correct the solution 
    %tic
    x = smooth(A,b,x,M,smooth_it);      % Post-smoothing
    %stime = stime + toc;
    
    err = norm(b-A*x);
    relres = err/normb;
    resvec = [resvec;err];
    iter = iter + 1;
end
    function x = smooth(A,b,x,M,smooth_it)
        for i = 1:smooth_it   %GS iterations applied to x
            r = b - A*x;      %for smooth_it times
            x = x + M\r;
        end
    end
end