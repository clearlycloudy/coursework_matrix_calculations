function [u_gmres,flag,relres,iter] = solver_gmres_bv_2ndorder

n=8;
dim=n-1;
h=1/n;

%a u_xxxx + b u_xxyy + c u_yyyy + d u_xx + e u_yy + f u = g

x = sym('x');
y = sym('y');
a = sym('a');
b = sym('b');
c = sym('c');
d = sym('d');
e = sym('e');
f = sym('f');

% a(x,y)=1+exp(x+y);
% b(x,y)=1+1/(x+y);
% c(x,y)=3+x*y;
% d(x,y)=-sin(x)*sin(y);
% e(x,y)=-1-exp(x+y);
% f(x,y)=x+y;

a(x,y)=1;
b(x,y)=2;
c(x,y)=1;
d(x,y)=0;
e(x,y)=0;
f(x,y)=0;

u = sym('u');

% u(x,y)=x^(9/2)*y^(9/2);
u(x,y)=sin(x)*cos(y);

u_xxxx = diff(u,x,4);
u_yyyy = diff(u,y,4);
u_xxyy = diff(diff(u,x,2),y,2);
u_xx = diff(u,x,2);
u_yy = diff(u,y,2);

%boundary values
eta = sym('eta');
zeta = sym('zeta');

eta = diff(u,y,2);
zeta = diff(u,x,2);

gamma = sym('gamma');
gamma = u;

g = sym('g');
g(x,y) = a(x,y) * u_xxxx(x,y) + b(x,y) * u_xxyy(x,y) + c(x,y) * u_yyyy +...
         d(x,y) * u_xx(x,y) + e(x,y) * u_yy(x,y) + f(x,y) * u(x,y);

[A,g] = create_A_g_generic(n,x,y,a,b,c,d,e,f,g,...
    u_xxxx,u_xxyy,u_yyyy,u_xx,u_yy,eta,zeta,gamma);

%solver: banded LU
u_ref = A\g;

%solver: gmres
maxit = 100;
tol = 10^(-9);
restart = 20;

[ret_u,ret_flag,ret_relres,ret_iter] = gmres(A,g,restart,tol,maxit,@precond);

fprintf("abs max error to reference banded LU: %f\n", max(abs(ret_u-u_ref)));

u_gmres = ret_u;
flag = ret_flag;
relres = ret_relres;
iter = ret_iter;

    function u_solve = precond(r)
        h = 1/n;
        dim = n-1;
        [AA, B, T] = create_A(n);
        %compute eta
        f_eta = diff(u,y,2);
        %compute zeta
        f_zeta = diff(u,x,2);
        %compute boundary values using 2nd order derivatives and 0th order values
        g2 = compute_boundary_value_constants(n,f,f_eta,f_zeta);
        %compute g(4th order mixed derivatives) on interior grid points
        f_4th = diff(u,x,4) + 2*diff(diff(u,x,2), y, 2) + diff(u,y,4);
        g1_grid = compute_grid(n,f_4th,x,y);
        assert((size(g1_grid,1)==n-1) && (size(g1_grid,2)==n-1));
        %remap grid points to 1D matrix
        g1 = zeros(n-1*n-1,1);
        for i=1:n-1
            for j=1:n-1
                g1(ij_to_k(dim, i, j),1)=g1_grid(i, j);
            end
        end
        g = g1+g2;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % dst + block LU solver:
        % equivalent to:
        %     u = kron(F^-1_n, I_m) BlockDiag^-1 kron(F_n,I_m)g
        %     BlockDiag := block diagonalization of A
        I = eye(dim);
        j = 1:1:n-1;
    
        %compute eigenvalues of tridiag{1,-2,1}
        eigenvalues = -4*(sin(j*pi/(2*n))).^2;
        D_T = spdiags(eigenvalues', 0:0, dim, dim);
        % alternative: use eigs
        % need to reverse order to match up with ordering of dst, idst
        %D_T = spdiags([flip(eigs(T,dim))], 0:0, dim, dim); %alternative
    
        D_B = D_T^2;
        BlockDiag = (1/h^4)*(kron(I,T*T)+...
                     kron(D_B,I)+...
                     2*(kron(D_T,T)));
                 
        %shape g to m by n matrix, 
        %m being y-axis (fastest increasing axis)
        %of original grid
        g_m_by_n = reshape(r,[dim,dim]);
        %transform by eigenvectors
        g_1_n_by_m = dst(g_m_by_n');
        g_1_stacked = reshape(g_1_n_by_m',[dim*dim,1]);
        %solve BlockDiag g_2 = g_1
        g_2_stacked = BlockDiag\g_1_stacked;
        g_2_m_by_n = reshape(g_2_stacked, [dim, dim]);
        %inverse transform by eigenvectors
        u_n_by_m = idst(g_2_m_by_n');
        u_solve = reshape(u_n_by_m',[dim*dim,1]);
    %     %convert to cartesian coordinate
    %     u = zeros(dim,dim);
    %     for c = 1:dim*dim
    %         [i,j] = k_to_ij(dim,c);
    %         u(i,j) = u_solve(c);
    %     end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

end



