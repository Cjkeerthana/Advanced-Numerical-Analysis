function [x, resvec, relres, iter] = Vcycle(A, b, x0, tol, maxit, smooth_it)
% V Cycle: Multi Grid Method
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
[P,R] = prolongation(A); % calculation of prolongation &
                         % restriction matrices
                         
Ac = R * A * P;          % Galerkin coarse grid operator
n = size(Ac,1);

x = x0;
err = norm(b - A*x);
normb = norm(b);
resvec = [];
iter = 0;
relres = err/normb;
pcg_maxit = 5000;
M = tril(A);

while err > normb*tol && iter < maxit
    x = smooth(A,b,x,M,smooth_it);    % Pre-smoothing
    
    rc = R *(b-A*x);      % Calculation of restricted residual
    
    if sqrt(n) <= 100     % If the coarse grid has less than
        opts.droptol = 10^(-3);  % or equal to 200 nodes, solve for the restricted error
        Lc = ichol(Ac,opts); %IC factorization with droptol 10(-3)
        [ec,~,~,~,~] = pcg(Ac,rc,tol,pcg_maxit,Lc,Lc');                           
    else
       ec = Vcycle(Ac,rc,zeros(n,1),tol,maxit,smooth_it); %else recurisvely
                                                %pass to Vcycle
                                                %Vcycle(Ac,rc,0)
    end
    
    x = x + P*ec;          % Prolong the error and
                           % correct the solution
                           
    x = smooth(A,b,x,M,smooth_it);     % Post-smoothing
    
    err = norm(b-A*x);
    relres = err/normb;
    resvec = [resvec;err];
    iter = iter + 1;
end
    function x = smooth(A,b,x,M,smooth_it)
        for i = 1:smooth_it         %GS iterations applied to x
            r = b - A*x;            %for smooth_it times
            x = x + M\r;
        end
    end
end