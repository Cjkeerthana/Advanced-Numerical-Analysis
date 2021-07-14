maxit = 100;
pcg_maxit = 5000;
%smooth_it = 5;
tol = 10^(-8);
for smooth_it = 1:2:9
    fprintf("Parameter: smoothing iterations = %d \n\n",smooth_it);
    time_twogrid = [];
    time_vcycle = [];
    time_pcg = [];
    prolongation_time = [];
    nodes = [];
    for nx = 200:200:1000
        fprintf("Parameter: nx = %d \n",nx);
        nodes = [nodes;nx];
        A = delsq(numgrid('S',nx+2));
        n=size(A,1);
        exactsol = zeros(n,1);
        for i = 1:n
          exactsol(i)=1/sqrt(i);
        end
        b = A*exactsol;
        x0 = zeros(n,1);   

        tic
        [P,R] = prolongation(A);
        prolongation_time = [prolongation_time;toc];

        tic
        [x2,resvec2,relres2,iter2] = twogrid(A, b, x0, tol, maxit,smooth_it);
        time_twogrid = [time_twogrid;toc];
        fprintf("Two Grid cycle:\tRELRES= %d;\tITER= %d;\tCPU= %0.2f\n",relres2,iter2,toc);

        tic
        [x3,resvec3,relres3,iter3] = Vcycle(A, b, x0, tol, maxit,smooth_it);
        time_vcycle = [time_vcycle;toc];
        fprintf("V cycle:\t\tRELRES= %d;\tITER= %d;\tCPU= %0.2f\n", relres3,iter3,toc);

        tic
        opts.droptol = 10^(-3);
        L = ichol(A,opts);
        [x4,flag4,relres4,iter4,resvec4] = pcg(A,b,tol,pcg_maxit,L,L');
        time_pcg = [time_pcg;toc;];
        fprintf("PCG:\t\t\tRELRES= %d;\tITER= %d;\tCPU= %0.2f\n",relres4,iter4,toc);
        fprintf("solution for nx = %d completed \n",nx);
        fprintf("************************************\n");
    end
    fprintf("solution for smoothing iterations = %d completed \n",smooth_it);
    fprintf("************************************\n");
    
    figure(smooth_it)
    plot(nodes, time_vcycle, '-o');
    hold on
    plot(nodes, time_twogrid, '-o');
    plot(nodes, time_pcg, '-o');
    plot(nodes, prolongation_time, '-o');
    str = ['Smoothing iterations = ',num2str(smooth_it)];
    title(str);
    xlabel('No. of nodes');
    ylabel('Time in seconds');
    legend('vcycle','twogrid','pcg','prolongation');
    fprintf("************************************\n\n");
end
clear;
