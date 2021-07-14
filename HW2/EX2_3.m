%%
%Adavnced Numerical Analysis - HW 2 - Keerthana C J
%%
%exercise 2 & 3
expected_it=[];

cg = [];
pcg1 = [];
pcg2 = [];
pcg3 = [];
jpcg = [];

cg_time = [];
pcg1_time = [];
pcg2_time = [];
pcg3_time = [];
jpcg_time = [];

for nx=100:100:400
    A = delsq(numgrid('S',nx+2));
    n=size(A,1);
    exactsol = zeros(n,1);
    for i = 1:n
        exactsol(i)=1/sqrt(i);
    end
    b = A*exactsol;
    
    tol = 10^(-8);
    maxit=2000;
    
    tic
    [x,flag1,relres,iter1,resvec1]=pcg(A,b,tol,maxit); %CG method
    cg_time=[cg_time;toc];
    cg=[cg;iter1];
    
    L=ichol(A);
    tic
    [x,flag2,relres, iter2,resvec2]=pcg(A,b,tol,maxit,L,L'); %PCG with IC(0) preconditioner
    pcg1_time=[pcg1_time;toc];
    pcg1=[pcg1;iter2];
    
    opts.type = 'ict';
    opts.droptol = 10^(-2);
    L=ichol(A,opts);
    tic
    [x,flag3,relres, iter3,resvec3]=pcg(A,b,tol,maxit,L,L'); %PCG with ICT preconditoner 
    pcg2_time=[pcg2_time;toc];                               % tol = 10^-2
    pcg2=[pcg2;iter3];
    
    opts.droptol = 10^(-3);
    L=ichol(A,opts);
    tic
    [x,flag4,relres, iter4,resvec4]=pcg(A,b,tol,maxit,L,L'); %PCG with ICT preconditioner
    pcg3_time=[pcg3_time;toc];                               % tol = 10^-3
    pcg3=[pcg3;iter4];
    
    % Exercise 3
    L=diag(diag(A));
    tic
    [x,flag5,relres, iter5,resvec5]=pcg(A,b,tol,maxit,L,[]); %PCG with Jacobi Preconditioner
    jpcg_time=[jpcg_time;toc];                              
    jpcg=[jpcg;iter5];
    
    expected_it=[expected_it;nx];
end
%Table containing the  iterations and cpu times for each method
T = table(expected_it,cg,cg_time,pcg1,pcg1_time,pcg2,pcg2_time,pcg3,pcg3_time,jpcg,jpcg_time);
display(T);
clear;