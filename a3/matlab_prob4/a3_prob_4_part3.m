ns = [8,16,32,64,128];

diff = zeros(length(ns),1);

u_prev=[];
u_LU_prev=[];

for i=1:length(ns)
    n = ns(i);
    dim = n-1;
    g = ones(dim*dim,1);
    
    fprintf("n: %d\n",n);
    
    [A, u] = my_precond(n,g);
    %0 0 0 0 0 0 0 0 0%
    %  1 2 3 4 5 6 7  %

    mid = n/2;
    k = ij_to_k(dim,mid,mid);

    u_LU = A\g;
    
    u_mid_fst = u(k);
    u_mid_lu = u_LU(k);
    
    fprintf("FST u(midpoint): %.12f\n", u_mid_fst);
    fprintf("LU u(midpoint): %.12f\n", u_mid_lu);
end

%n=8: 0.001880690688

0.001880690688 

