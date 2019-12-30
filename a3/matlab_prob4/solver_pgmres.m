function [u_gmres,flag,relres,iter,resvec,A,g] = solver_pgmres(n,u,x,y,a,b,c,d,e,f)
%assume u(x,y),a,b,c,d,e,f are symbolic functions

dim=n-1;
h=1/n;

%a u_xxxx + b u_xxyy + c u_yyyy + d u_xx + e u_yy + f u = g

u_xxxx = diff(u,x,4);
u_yyyy = diff(u,y,4);
u_xxyy = diff(diff(u,x,2),y,2);
u_xx = diff(u,x,2);
u_yy = diff(u,y,2);

%boundary values
eta = sym('eta');
zeta = sym('zeta');

eta = diff(u,y,1);
zeta = diff(u,x,2);

gamma = sym('gamma');
gamma = u;

g = sym('g');
g(x,y) = a(x,y) * u_xxxx(x,y) + b(x,y) * u_xxyy(x,y) + c(x,y) * u_yyyy +...
         d(x,y) * u_xx(x,y) + e(x,y) * u_yy(x,y) + f(x,y) * u(x,y);

[A,g] = create_A_g_generic(n,x,y,a,b,c,d,e,f,g,...
    u_xxxx,u_xxyy,u_yyyy,u_xx,u_yy,eta,zeta,gamma);

%sanity check
% u_cart = convert_to_cartesian(dim, u_ref);
% u_check = full(compute_grid(n,u,x,y));
% max(max(abs(u_cart-u_check)))

%solver: gmres
maxit = 50;
tol = 10^(-9);
restart = 20;

[ret_u,ret_flag,ret_relres,ret_iter,resvec] = gmres(A,g,restart,tol,maxit,@precond);

u_gmres = ret_u;
flag = ret_flag;
relres = ret_relres;
iter = ret_iter;

    function u_solve = precond(r)

        h = 1/n;
        dim = n-1;
        [AA, B, T, C] = create_A_2(n);
        % dst + block LU solver:
        % equivalent to:
        %     u = kron(F^-1_n, I_m) BlockDiag^-1 kron(F_n,I_m)g
        %     BlockDiag := block diagonalization of A
        I = eye(dim);
        j = 1:1:n-1;
    
        %compute eigenvalues of tridiag{1,-2,1}
        eigenvalues = -4*(sin(j*pi/(2*n))).^2;
        D_T = spdiags(eigenvalues', 0:0, dim, dim);
    
        % kron(I,C) added below for the augmentation of 
        % changed boundary condition from u_yy to u_y
        D_B = D_T^2;
        BlockDiag = (1/h^4)*(kron(I,T*T)+...
                             kron(D_B,I)+...
                             2*(kron(D_T,T))+...
                             kron(I,C));

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
    end
end