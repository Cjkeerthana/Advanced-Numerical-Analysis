%%
%Exercise 6:
%Solving a 2D Bratu problem using Newton's, Bryoden's and Inexact Newton
%methods
%%
h = 5e-3;
nx = (1/h) - 1.0;
B = delsq(numgrid('S',nx+2));
n = size(B,1);
x0 = 0.1*ones(n,1);
tol = 1e-13;
itmax = 50;
lambda = 6.5;
ffun = @(x)(-B*x+h*h*lambda*exp(x)); % F(x)
jfun = @(x)(-B+spdiags((h*h*lambda*exp(x)),0,n,n)); %Jacobian

%Newtons with LU factorization
tic
[nxstar,niter,nresvec,nstepdiff] = newton(x0,@(x)ffun(x),@(x)jfun(x),tol,itmax,-1);
time_n = toc;

%Broyden's method
tic
[qnxstar,qniter,qnresvec,qnstepdiff] = quasinewton(x0,@(x)ffun(x),@(x)jfun(x),tol,itmax);
time_qn = toc;
eta_max = 0.1;
time_in = [];
inxstar = [];
initer = [];
inresvec = [];
instepdiff = [];
inliter = [];
inft = [];
% Inexact Newtom with different forcing parameters
for opt_eta = 1:4
    tic
    [xstar,iter,resvec,stepdiff,liter,ft] = inexactnewton(x0,@(x)ffun(x),@(x)jfun(x),eta_max,opt_eta,itmax);
    time_in = [time_in;toc];
    initer = [initer;iter];
    inresvec = [inresvec;resvec];
    instepdiff = [instepdiff;stepdiff];
    inliter = [inliter;liter];
    inft = [inft;ft];
end

to = 0;
for i = 1:4
    from = to +1;
    to = from + initer(i);
    k = [0:to-from]';
    forcingterm = inft(from:to);
    residual = inresvec(from:to);
    liniter = inliter(from:to);
    T = table(k,forcingterm,residual,liniter);
    disp(T);
end

figure(1)
semilogy(0:length(nresvec)-1,nresvec,'-o')
hold on
semilogy(0:length(qnresvec)-1,qnresvec,'-o')
yline(tol,'r--');
legend('Newton','Broyden')
xlabel('Iteration number')
ylabel('Non-linear Residual')

figure(2)
semilogy(0:initer(1),inresvec(1:initer(1)+1),'-o')
hold on
to = initer(1)+1;
for i = 2:4
    from = to +1;
    to = from + initer(i);
    semilogy(0:initer(i),inresvec(from:to),'-o')
end
yline(tol,'r--');
legend('Newton','opt1','opt2','opt3','opt4')
xlabel('Iteration number')
ylabel('Non-linear Residual')
clear;
