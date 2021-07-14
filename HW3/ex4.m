%%
%Exercise 4 - Using the restart parameter for the GMRES function with
%preconditioner. The effect of the restart parameter is studied.
%%
A = readmatrix('mat13041.rig.txt'); %Input matrix
A = spconvert(A);
[n,n] = size(A);
exact = zeros(n,1);
for i = 1:n
    exact(i)=1/sqrt(i);
end
b = A*exact;
x0 = zeros(n,1);
maxit = 50;
tol = 10^(-12);

%ILU factorization
setup.type = 'crout';
setup.droptol = 0.01;
[L,U] = ilu(A,setup);

iter_L = [];
relres_L =[];
left_time = [];
resvec_L=[];
resvec_L_length=[];

iter_R = [];
relres_R =[];
resvec_R=[];
right_time = [];
resvec_R_length=[];

iter_S = [];
relres_S =[];
resvec_S=[];
split_time = [];
resvec_S_length=[];

for restart = [10,20,30,50,100]
    tic
    [x,flag,relres,iter,resvec]=gmres(A,b,restart,tol,maxit,L,U); %Left preconditioning
    left_time = [left_time;toc];
    temp = (iter(1)-1)*restart+iter(2);
    iter_L = [iter_L;temp];
    relres_L = [relres_L;relres];
    resvec_L=[resvec_L;resvec];
    resvec_L_length=[resvec_L_length;length(resvec)];
    
    tic
    [x,flag,relres,iter,resvec]=gmres(@(x)A*(U\(L\x)),b,restart,tol,maxit,L,U);%Right preconditioning
    right_time = [right_time;toc];
    temp = (iter(1)-1)*restart+iter(2);
    iter_R = [iter_R;temp];
    relres_R = [relres_R;relres];
    resvec_R=[resvec_R;resvec];
    resvec_R_length=[resvec_R_length;length(resvec)];
    
    tic
    [x,flag,relres,iter,resvec]=gmres(@(x)L\(A*(U\x)),L\b,restart,tol,maxit,L,U);%Split preconditioning
    split_time = [split_time;toc];
    temp = (iter(1)-1)*restart+iter(2);
    iter_S = [iter_S;temp];
    relres_S = [relres_S;relres];
    resvec_S=[resvec_S;resvec];
    resvec_S_length=[resvec_S_length;length(resvec)];
end
%plots
%Left Preconditioning
figure(1)
semilogy(0:(resvec_L_length(1)-1),resvec_L(1:resvec_L_length(1))/norm(U\(L\b)),'-o')
hold on
from = 0;
to = resvec_L_length(1);
for i = [2,3,4,5]
    from = from + resvec_L_length(i-1);
    to = from + resvec_L_length(i);
    semilogy(0:(resvec_L_length(i)-1),resvec_L(from+1:to)/norm(U\(L\b)),'-o')
end
yline(tol,'r--');
legend('10','20','30','50','100','tol')
xlabel('Iteration number')
ylabel('Relative Residual')
title('Left Preconditioning')

%Right Preconditioning
figure(2)
semilogy(0:(resvec_R_length(1)-1),resvec_R(1:resvec_R_length(1))/norm(U\(L\b)),'-o')
hold on
from = 0;
for i = [2,3,4,5]
    from = from + resvec_R_length(i-1);
    to = from + resvec_R_length(i);
    semilogy(0:(resvec_R_length(i)-1),resvec_R(from+1:to)/norm(U\(L\b)),'-o')
end
yline(tol,'r--');
legend('10','20','30','50','100','tol')
xlabel('Iteration number')
ylabel('Residual')
title('Right Preconditioning')

%Split preconditioning
figure(3)
semilogy(0:(resvec_S_length(1)-1),resvec_S(1:resvec_S_length(1))/norm(L\(U\b)),'-o')
hold on
from = 0;
for i = [2,3,4,5]
    from = from + resvec_S_length(i-1);
    to = from + resvec_S_length(i);
    semilogy(0:(resvec_S_length(i)-1),resvec_S(from+1:to)/norm(L\(U\b)),'-o')
end
yline(tol,'r--');
legend('10','20','30','50','100','tol')
xlabel('Iteration number')
ylabel('Residual')
title('Split Preconditioning')

%Values
restart = [10;20;30;50;100];
T = table(restart,iter_L,relres_L,left_time,iter_R,relres_R,right_time,iter_S,relres_S,split_time);
display(T);
clear;