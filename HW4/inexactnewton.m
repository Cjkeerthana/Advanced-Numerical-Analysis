function [xstar,iter,resvec,stepdiff,liter,ft] = inexactnewton(x0,F,Jac,eta_max,opt_eta,itmax)
% [xstar,iter,resvec,stepdiff] = quasinewton(x0,F,Jac,tol,itmax)
% OUTPUT parameters: 
% xstar: solution vector for all iteration steps
% iter: number of outer iterations
% resvec: vector containing the non-linear residual norm : ||F(x_k)||
% stepdiff: norm of each step: ||s||
% liter: total number of linear iterations
% ft: forcing term in each iteration
% INPUT parameters:  
% x0: initial guess 
% F: The non-linear system of equations F(x) given as a function handle
% returning a vector
% Jac: The Jacobian of the system F'(x) which returns a matrix of n by n
%      The initial B0 matrix = F'(x)
% tol: tolerance for the exit test (||F(x_k)||/||F(x0)||)
% eta_max: Forcing parameter
% opt_eta: Varying the Forcing parameters, options 1-4
% itmax: maximum number of outer iterations
x = x0;
resvec = [norm(F(x0))];
stepdiff = [];
xstar = [x0];
iter = 0;
tol = 1e-13;
err = norm(F(x0));
gm_tol = eta_max;
liter = [0];
ft = [eta_max];
while err > tol*norm(F(x0)) && iter < itmax       
    restart = 50;
    gm_itmax = 20;
    setup.type='ilutp';
    setup.droptol = 0.01;
    [L,U] = ilu(-Jac(x),setup);
    [s,gmflag,gmrelres,gmiter,gmresvec] = gmres(-Jac(x),F(x),restart,gm_tol,gm_itmax,L,U);
    liter = [liter;(gmiter(1)-1)*restart+gmiter(2)];
    xold = x;
    x = x + s;
    xstar = [xstar,x];
    resvec = [resvec;norm(F(x))];
    stepdiff = [stepdiff;norm(s)];
    %options to vary the tolerance for solving the linear system
    if opt_eta == 1
        gm_tol = eta_max;
    elseif opt_eta == 2
        gm_tol = gm_tol/3;
    elseif opt_eta == 3
        t = 0.95*norm(F(x));
        gm_tol = min(eta_max,t);
    else
        t = 0.95*(norm(F(x)))^2/(norm(F(xold)))^2;
        gm_tol = min(eta_max,t);
    end
    err = norm(F(x));
    ft = [ft;gm_tol];
    iter = iter + 1;
end
end