%%
%Exercise 5 - The density of the ILU preconditioner and the dependence of
%the GMRES convergence is studied here. The total time for both the ILU
%factorization and the GMRES convergence gives us more insight.
%%
load ML_laplace.mtx; %Input Matrix
A = spconvert(ML_laplace);
[n,n] = size(A);
b = A*ones(n,1);
x0 = zeros(n,1);
maxit = 550;
tol = 10^(-10);
restart = 50;

tprec = [];
tsol = [];
tcpu = [];
resvec_whole = [];
rho_whole = [];
iters = [];
rs_len = [];
relres_whole = [];

for i = [0.02,0.01,0.003,0.001,10^(-4),10^(-5)]
    tic
    setup.type = 'crout';
    setup.droptol = i;
    [L,U] = ilu(A,setup); %ILU preconditioner
    tprec = [tprec;toc];
    
    %density of the preconditioner
    rho = (nnz(L)+nnz(U)-n)/nnz(A);
    rho_whole = [rho_whole;rho];
    
    tic
    [x,flag,relres,iter,resvec]=gmres(A,b,restart,tol,maxit,L,U); %Left Preconditioning 
    tsol = [tsol;toc];
    temp = (iter(1)-1)*restart+iter(2);
    resvec_whole = [resvec_whole; resvec];
    rs_len = [rs_len; length(resvec)];
    relres_whole = [relres_whole; relres];
    iters = [iters;temp];
    
    j=1;
    tcpu = [tcpu;tsol(j)+tprec(j)];
    j=j+1;
end

%plots
figure(1)
semilogy(0:(rs_len(1)-1),resvec_whole(1:rs_len(1))/norm(U\(L\b)),'-o')
hold on
from = 0;
to = rs_len(1);
for i = [2,3,4,5,6]
    from = from + rs_len(i-1);
    to = from + rs_len(i);
    semilogy(0:(rs_len(i)-1),resvec_whole(from+1:to)/norm(U\(L\b)),'-o')
end
yline(tol,'r--');
legend('0.02','0.01','0.003','10^-^4','10^-^5')
xlabel('Iteration number')
ylabel('Relative Residual')
title('Left Preconditioning')

%table of all times and comparisons
droptol = [0.02;0.01;0.003;0.001;10^(-4);10^(-5)];
T=table(droptol,tprec,tsol,tcpu,iters,relres_whole,rho_whole);
display(T);
clear;