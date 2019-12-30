function [AA, u_solve] = my_precond(n,g)
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
    
    BlockDiag = (1/h^4)*(kron(I,T*1/T)+...
                         kron(D_B,I)+...
                         2*(kron(D_T,T))+...
                         kron(I,C));

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
end