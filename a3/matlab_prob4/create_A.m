function [A, B, T] = create_A(n)
    
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
    
    %matrix corresponding to (u_xx+u_yy)^2
    A = 1/(h^4) * B^2;
end
    
