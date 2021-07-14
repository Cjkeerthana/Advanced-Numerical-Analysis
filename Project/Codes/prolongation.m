function [P,R] = prolongation(A)
    nx = sqrt(size(A,1));
    x = [];
    y = [];
    v = [];
    j = 1;
    l = 2;
    if mod(nx,2) == 0
         nxc = nx/2;
    else
         nxc = (nx-1)/2;
    end
    for i = 1: nxc*nxc
        for k = 0:2
             if (k*nx+j > nx*nx)
                break;
             elseif (mod(i,nxc) == 0) && (mod(nx,2) == 0)
                temp = [i;i];
                x = [x;temp];
                temp = [k*nx+j; k*nx+j+1];
                y = [y;temp];
                if k == 1
                    v = [v;2;4];
                else
                    v = [v;1;2];
                end
             else 
                 temp = [i;i;i];
                 x = [x;temp];
                 temp = [k*nx+j; k*nx+j+1; k*nx+j+2];
                 y = [y;temp];
                 if k == 1
                    v = [v;2;4;2];
                 else
                    v = [v;1;2;1];
                 end
             end
        end
        if mod(i,nxc) == 0
            j = l*nx+1;
            l = l + 2;
        else
             j = j+2;
        end
    end
    R = sparse(x,y,v);
    P = R';
end