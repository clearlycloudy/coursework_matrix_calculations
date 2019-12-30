%using composition of operators:
%u_xxxx+u_yyyy+2u_xxyy = (u_xx+u_yy)^2

clear
dim = 5;
I = eye(dim);

%construct tridiagonal matrix used for u_xx, u_yy
%for each of dimensions x,y
C = spdiags([ones(dim,1) -2*ones(dim,1) ones(dim,1)],...
            -1:1,...
            dim,dim);

%let fastest increasing index to be along y dimension

%kronecker product corresponding to difference operator u_yy
B_yy = kron(I,C);

%kronecker product corresponding to difference operator u_xx
B_xx = kron(C,I);

%matrix corresponding to (u_xx+u_yy)
B = B_yy + B_xx;

%matrix corresponding to (u_xx+u_yy)^2
A = B^2;

