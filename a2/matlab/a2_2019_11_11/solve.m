function u = solve(n, f, x, y)
% n := number of intervals in grid
% assumes f is a symbolic expression in terms of x and y
% returns approximation of u on interior of grid in cartesian coordinates
h = 1/n;
dim = n-1;
[A, B, T] = create_A(n); %A is actually not needed, use 
%compute eta
f_eta = diff(f,y,2);
%compute zeta
f_zeta = diff(f,x,2);
%compute boundary values using 2nd order derivatives and 0th order values
g2 = compute_boundary_value_constants(n,f,f_eta,f_zeta);
%compute g(4th order mixed derivatives) on interior grid points
f_4th = diff(f,x,4) + 2*diff(diff(f, x, 2), y, 2) + diff(f,y,4);
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
g_m_by_n = reshape(g,[dim,dim]);
%transform by eigenvectors
g_1_n_by_m = dst(g_m_by_n');
g_1_stacked = reshape(g_1_n_by_m',[dim*dim,1]);
%solve BlockDiag g_2 = g_1
g_2_stacked = BlockDiag\g_1_stacked;
g_2_m_by_n = reshape(g_2_stacked, [dim, dim]);
%inverse transform by eigenvectors
u_n_by_m = idst(g_2_m_by_n');
u_solve = reshape(u_n_by_m',[dim*dim,1]);
%convert to cartesian coordinate
u = zeros(dim,dim);
for c = 1:dim*dim
    [i,j] = k_to_ij(dim,c);
    u(i,j) = u_solve(c);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%

