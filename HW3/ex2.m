%%
%Exercise 3 - Implementation of preconditioned GMRES and solving mat13041.rig The
%following exercise uses the function myprecgmres.m where the preconditioned GMRES algorithm
%has been implemented for all the three types of preconditoners (Left,
%Right and Split) and solves the equation Ax=b. The results are
%compared with matlab inbuilt gmres()function with respective
%preconditioners.
%An ILU precondtioner is used for all the comparisons.
%%
A = readmatrix('mat13041.rig.txt');
A = spconvert(A);
[n,n] = size(A);
exact = zeros(n,1);
for i = 1:n
    exact(i)=1/sqrt(i);
end
b = A*exact;
x0 = zeros(n,1);
maxit = 800;
tol = 10^(-10);

setup.type = 'crout';
setup.droptol = 0.1;
[L,U] = ilu(A,setup); %% ILU precondtioner

tic
ptype = 'L';
[x,iter,resvec,flag]=myprecgmres(A,b,tol,maxit,x0,ptype,L,U); %myprecgmres with Left preconditioning
myprecgmres_L_time = toc;
iter_L = iter;
relres_L = norm(U\(L\(b-A*x)))/norm(U\(L\b));

tic
ptype = 'R';
[x2,iter,resvec2,flag2]=myprecgmres(A,b,tol,maxit,x0,ptype,L,U); %myprecgmres with Right preconditioning
myprecgmres_R_time = toc;
iter_R = iter;
relres_R = norm(b-A*x2)/norm(b);

tic
ptype = 'S';
[x3,iter,resvec3,flag3]=myprecgmres(A,b,tol,maxit,x0,ptype,L,U); %myprecgmres with Split preconditioning
myprecgmres_S_time = toc;
iter_S = iter;
relres_S = norm(L\(b-A*x3))/norm(L\b);

restart = 0;
tic
[x4,flag4,relres,iter,resvec4]=gmres(A,b,[],tol,maxit,L,U); %matlab left preconditoning
matlab_left_time = toc;
mat_iter_L = (iter(1)-1)*restart+iter(2);
mat_relresL = relres;

tic
[x5,flag5,relres,iter,resvec5]=gmres(@(x)A*(U\(L\x)),b,[],tol,maxit); %matlab right preconditoning
matlab_right_time = toc;
mat_iter_R = (iter(1)-1)*restart+iter(2);
x5 = U\(L\x5);
mat_relresR = relres;

tic
[x6,flag6,relres,iter,resvec6]=gmres(@(x)L\(A*(U\x)),L\b,[],tol,maxit); %matlab split preconditioning
matlab_split_time = toc;
mat_iter_S = (iter(1)-1)*restart+iter(2);
x6 = U\x6;
mat_relresS = relres;

%True residuals
res_L = norm(b-A*x)/norm(b);
res_R = norm(b-A*x2)/norm(b);
res_S = norm(b-A*x3)/norm(b);
mat_resL = norm(b-A*x4)/norm(b);
mat_resR = norm(b-A*x5)/norm(b);
mat_resS = norm(b-A*x6)/norm(b);

figure(1)
semilogy(0:length(resvec)-1,resvec/norm(U\(L\b)),'-o')
hold on
semilogy(0:length(resvec2)-1,resvec2/norm(b),'-o')
semilogy(0:length(resvec3)-1,resvec3/norm((L\b)),'-o')
yline(tol,'r--');
legend('myprec_L','myprec_R','myprec_S','tol')
xlabel('Iteration number')
ylabel('Preconditioned Relative Residual')

figure(2)
semilogy(0:length(resvec4)-1,resvec4/norm(U\(L\b)),'-o')
hold on
semilogy(0:length(resvec5)-1,resvec5/norm(b),'-o')
semilogy(0:length(resvec6)-1,resvec6/norm(L\b),'-o')
yline(tol,'r--')
legend('matlab_L','matlab_R','matlab_S','tol')
xlabel('Iteration number')
ylabel('Preconditioned Relative Residual')

%Table comparing the iterations and timings
preconditioner = ["myprecgmres_left";"matlab_left";"myprecgmres_right";"matlab_right";"myprecgmres_split";"matlab_split"];
iterations = [iter_L;mat_iter_L;iter_R;mat_iter_R;iter_S;mat_iter_S];
real_residuals = [res_L;mat_resL;res_R;mat_resR;res_S;mat_resS];
preconditioned_residuals = [relres_L;mat_relresL;relres_R;mat_relresR;relres_S;mat_relresS];
timings = [myprecgmres_L_time; matlab_left_time;myprecgmres_R_time;matlab_right_time;myprecgmres_S_time;matlab_split_time];
T = table(preconditioner,iterations,timings,real_residuals,preconditioned_residuals);
display(T);
clear;