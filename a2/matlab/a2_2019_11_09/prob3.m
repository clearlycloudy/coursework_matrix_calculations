%using composition of operators:
%u_xxxx+u_yyyy+2u_xxyy = (u_xx+u_yy)^2

clear

n = 8;
h = 1/n;
dim = n-1;

[A, B, T] = create_A(n);

A = A .* (1/h^4);

%boundary conditions, offset by 1 due to 1-indexing
eta = sparse(n+1,n+1); %along y dimension
zeta = sparse(n+1,n+1); %along x dimension

%grid for 4th derivatives
gamma = sparse(n+1,n+1);

% %test with dummy values for boundary
% %indexing for boundary: (x,y)
% gamma(1:n+1,1)=1;
% gamma(1:n+1,n+1)=1;
% gamma(1,1:n+1)=1;
% gamma(n+1,1:n+1)=1;
% eta(1:n+1,1)=50;
% eta(1:n+1,n+1)=50;
% zeta(1,1:n+1)=50;
% zeta(n+1,1:n+1)=50;

syms f(x,y);
f(x,y)=sin(x)+cos(y);
% f(x,y)=sin(x)*cos(y);

%compute gamma (0th order boundary values)
gamma = compute_boundary_value(n,f,x,y);
assert((size(gamma,1)==n+1) && (size(gamma,2)==n+1));

%compute eta
f_eta = diff(f,y,2);
eta = compute_boundary_value(n,f_eta,x,y);
assert((size(eta,1)==n+1) && (size(eta,2)==n+1));

%compute zeta
f_zeta = diff(f,x,2);
zeta = compute_boundary_value(n,f_zeta,x,y);
assert((size(zeta,1)==n+1) && (size(zeta,2)==n+1));

%compute boundary values using 2nd order derivatives and 0th order values
g2 = compute_boundary_value_constants(n,gamma,eta,zeta);
g2_check = compute_boundary_value_constants_2(n,gamma,eta,zeta);

%debug only
g2_debug = zeros(dim,dim);
for c = 1:dim*dim
    [i,j] = k_to_ij(dim,c);
    g2_debug(i,j) = g2(c);
end

%debug only
g2_check_debug = zeros(dim,dim);
for c = 1:dim*dim
    [i,j] = k_to_ij(dim,c);
    g2_check_debug(i,j) = g2_check(c);
end

%compute g(4th order mixed derivatives) on interior grid points
f_xxxx = diff(f,x,4);
f_yyyy = diff(f,y,4);
f_xxyy = diff(diff(f, x, 2), y, 2);
f_4th = f_xxxx + 2*f_xxyy + f_yyyy;
g1_grid = compute_grid(n,f_4th,x,y);
assert((size(g1_grid,1)==n-1) && (size(g1_grid,2)==n-1));
%remap grid points to 1D matrix
g1 = zeros(n-1*n-1,1);
for i=1:n-1
    for j=1:n-1
        g1(ij_to_k(dim, i, j),1)=g1_grid(i, j);
    end
end
%debug only
g1_debug = zeros(dim,dim);
for c = 1:dim*dim
    [i,j] = k_to_ij(dim,c);
    g1_debug(i,j) = g1(c);
end

g = g1+g2;

%debug only
g_debug = zeros(dim,dim);
for c = 1:dim*dim
    [i,j] = k_to_ij(dim,c);
    g_debug(i,j) = g(c);
end

u=g\A;
%debug only
u_debug = zeros(dim,dim);
for c = 1:dim*dim
    [i,j] = k_to_ij(dim,c);
    u_debug(i,j) = u(c);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% sanity checks

grid_0th = compute_grid(n,f,x,y);
grid_0th_debug = full(grid_0th);

test11=eval(-f_zeta(0,h)*1/h^2-f_eta(h,0)*1/h^2)-...
    1/h^4*(2*gamma(1,1)-6*gamma(2,1)-6*gamma(1,2)+2*gamma(3,1)+2*gamma(1,3))

test12=eval(-f_zeta(0,2*h)*1/h^2)-...
    1/h^4*(2*gamma(1,2)-6*gamma(1,3)+2*gamma(1,4)+1*gamma(2,1))

test13=eval(-f_zeta(0,3*h)*1/h^2)-...
    1/h^4*(2*gamma(1,3)-6*gamma(1,4)+2*gamma(1,5))

test15=eval(-f_zeta(0,5*h)*1/h^2-f_eta(h,6*h)*1/h^2)-...
    1/h^4*(2*gamma(1,5)-6*gamma(1,6)+2*gamma(1,7)-6*gamma(2,7)+2*gamma(3,7))

test55=eval(-f_zeta(6*h,5*h)*1/h^2-f_eta(5*h,6*h)*1/h^2)-...
    1/h^4*(2*gamma(7,5)-6*gamma(7,6)+2*gamma(7,7)-6*gamma(6,7)+2*gamma(5,7))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dst + block LU solver:
% equivalent to:
%     u = kron(F^-1_n, I_m) BlockDiag^-1 kron(F_n,I_m)g
%     BlockDiag := block diagonalization of A
I = eye(dim);
D_T = diag(eig(full(T)));
D_B = D_T^2;
BlockDiag = (1/h^4)*(kron(I,T*T)+...
                     kron(D_B,I)+...
                     2*(kron(D_T,T)));

%shape g to m by n matrix, 
%m being y-axis (fastest increasing axis)
%of original grid
g_m_by_n = reshape(g,[dim,dim]);
%transform by eigenvectors
g_1_n_by_m = dst(g_m_by_n');
g_1_stacked = reshape(g_1_n_by_m',[dim*dim,1]);
%solve BlockDiag g_2 = g_1
g_2_stacked = g_1_stacked \ BlockDiag;
g_2_m_by_n = reshape(g_2_stacked, [dim, dim]);
%inverse transform by eigenvectors
u_n_by_m = idst(g_2_m_by_n');
u_solve = reshape(u_n_by_m',[dim*dim,1]);

u_solve_debug = zeros(dim,dim);
for c = 1:dim*dim
    [i,j] = k_to_ij(dim,c);
    u_solve_debug(i,j) = u_solve(c);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%








