%%
%Adavnced Numerical Analysis - HW 2 - Keerthana C J
%%
% Exercise 6
for nx=200:200:400
    A=delsq(numgrid('S',nx+2));
    n=size(A,1);
    b = A * ones(n,1);
    tol = 10^(-8);
    maxit=500;
    L=ichol(A);
    [W,Lambda]=eigs(A,L*L',10,'sm'); % W matrix
    H=W'*A*W; % H matrix
    % Spectral preconditoner passed as Function handler
    [x,flag,relres, iter,resvec]=pcg(A,b,tol,maxit,@(x)(L*L')\x+W*(H\(W'*x))); 
    [x,flag2,relres, iter2,resvec2]=pcg(A,b,tol,maxit,L,L'); % PCG with IC(0) preconditioner

    %plots
    figure(nx/200)
    semilogy(0:length(resvec)-1,resvec/norm(b),'-o')
    hold on
    semilogy(0:length(resvec2)-1,resvec2/norm(b),'-o')
    yline(tol,'r--');
    legend('spectral','IC(0)')
    xlabel('Iteration number')
    ylabel('Residual')
    title(nx)
end
clear;
