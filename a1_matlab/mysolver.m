function x = mysolver(n, f_g, f_gamma, f_zeta)

    %assume domain is (0,1)
    dim = n-1;
    h=1/n;
    
    %Ax=b
    
    %construct A
    a = ones(dim,1)*6;
    b = ones(dim,1)*-4;
    c = ones(dim,1)*1;
    A = spdiags([c b a b c], -2:2, dim, dim);
    A(1,1) = 5;
    A(dim,dim) = 5;
    A=A./(h^4);
    
    %get boundary conditions for u(x) at x=0, x=1
    gamma = f_gamma();
    gamma_0 = gamma(1);
    gamma_n = gamma(2);
    
    %get boundary conditions for u(x)'' at x=0, x=1
    zeta = f_zeta();
    zeta_0 = zeta(1);
    zeta_n = zeta(2);
    
    %construct constants arising from boundary conditions
    B0 = zeros(dim,1);
    B0(1,1) = 2*gamma_0/(h^4) - zeta_0/(h^2);
    B0(2,1) = -1*gamma_0/(h^4);
    B0(dim-1,1) = -1*gamma_n/(h^4);
    B0(dim,1) = 2*gamma_n/(h^4) - zeta_n/(h^2);
    
    %construct g(x) = u(x)'''' for non-boundary region
    i = 1:1:dim;
    domain_x = h*i;
    B1 = arrayfun(f_g,domain_x(1:dim))';
    assert(size(B0,1)==dim);
    assert(size(B1,1)==dim);
    
    b=B0+B1;
    
    %solve
    x=A\b;
    
end