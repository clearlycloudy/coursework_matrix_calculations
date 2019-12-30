function [A, B, T, C] = create_A_2(n)
    
    h=1/n;
    dim = n-1;
    I = eye(dim);

    %using composition of operators:
    %u_xxxx+u_yyyy+2u_xxyy = (u_xx+u_yy)^2
    
    %construct tridiagonal matrix used for u_xx, u_yy
    %for each of dimensions x,y
    T = spdiags([ones(dim,1) -2*ones(dim,1) ones(dim,1)],...
                -1:1,...
                dim,dim);

    %let fastest increasing index to be along y dimension

    %kronecker product corresponding to difference operator u_yy
    B_yy = kron(I,T);

    %kronecker product corresponding to difference operator u_xx
    B_xx = kron(T,I);

    %matrix corresponding to (u_xx+u_yy)
    B = B_yy + B_xx;
    
    %alter values near boundary due to change of boundary condition
    %from u_yy = eta to u_y = eta
    C = spdiags([zeros(dim,1)], 0:0, dim, dim);

    %here we change matrix to revert origin boundary condition
    %corresponding to u_yy and apply the new boundary condition 
    %for u_y with o(h^4) accuracy
    %apply 4 * eta B.C. at bottom (j=1)
    D_bottom = 4*[-1/4 -5/6 3/2 -1/2 5/60];
    %apply -4 * eta B.C. at top (j=n-1)
    D_top = -4*flip([1/4 5/6 -3/2 1/2 -5/60]);
    
    C(1,1) = 1; %revert old boundary condition in y direction
    %apply new u_y boundary condition at j=1
    C(1,1:3) = C(1,1:3) + D_bottom(3:5);
    C(dim,dim) = 1; %revert old boundary condition in y direction
    %apply new u_y boundary condition at j=n-1
    C(dim,dim-2:dim) = C(dim,dim-2:dim) + D_top(1:3);

    CC = kron(I,C);
    A = 1/(h^4) * (CC + B^2);

end
    
