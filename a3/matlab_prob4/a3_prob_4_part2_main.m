syms x y;
syms u(x,y);
u(x,y) = x^(9/2)*y^(9/2);

syms a b c d e f;
a(x,y)=1+exp(x+y);
b(x,y)=1+1/(x+y);
c(x,y)=3+x*y;
d(x,y)=-sin(x)*sin(y);
e(x,y)=-1-exp(x+y);
f(x,y)=x+y;

ns = [8,16,32,64,128];

errors = zeros(length(ns),1);

for i=1:length(ns)
    n = ns(i);
    dim = n-1;
    
    fprintf("n: %d\n",n);
    
    [u_gmres,flag,relres,iter,resvec,A,g] = solver_pgmres(n,u,x,y, ...
                                              a,b,c,d,e,f);

    semilogy(0:length(resvec)-1,resvec./norm(g),'-o');
    xlabel('Iteration number');
    ylabel('Relative residual');

    u_exact = full(compute_grid(n,u,x,y));
    u_exact_stacked = reshape(u_exact', [dim*dim,1]);

    %block LU
    u_block = A\g;

    fprintf("norm of difference of banded LU and PGMRES: %.12f\n", norm(u_gmres-u_block));
    fprintf("banded LU: max abs knot error: %.12f\n", max(abs(u_block-u_exact_stacked)));
    fprintf("PGMRES: max abs knot error: %.12f\n", max(abs(u_gmres-u_exact_stacked)));
    fprintf("PGMRES: converge flag: %d\n", flag);
    fprintf("PGMRES: iter:(%d, %d)\n", iter(1), iter(2));
    errors(i) = max(abs(u_gmres-u_exact_stacked));
end

order_of_convergence = zeros(length(ns)-1,1);

for i=1:length(ns)-1
    order_of_convergence(i)=log2(errors(i)/errors(i+1));
    fprintf("order of convergence: (n: %d, 2n: %d): %f\n", ...
        ns(i), ns(i+1), order_of_convergence(i));
end

%max abs knot error
% n: 8
% norm of difference of banded LU and PGMRES: 1.917901
% banded LU: max abs knot error: 0.399895
% PGMRES: max abs knot error: 0.000277
% PGMRES: converge flag: 0
% PGMRES: iter:(1, 11)
%
% n: 16
% norm of difference of banded LU and PGMRES: 3.737804
% banded LU: max abs knot error: 0.397527
% PGMRES: max abs knot error: 0.000064
% PGMRES: converge flag: 0
% PGMRES: iter:(1, 13)
%
% n: 32
% norm of difference of banded LU and PGMRES: 7.189019
% banded LU: max abs knot error: 0.397540
% PGMRES: max abs knot error: 0.000016
% PGMRES: converge flag: 0
% PGMRES: iter:(1, 15)
%
% n: 64
% norm of difference of banded LU and PGMRES: 14.820621
% banded LU: max abs knot error: 0.451666
% PGMRES: max abs knot error: 0.000004
% PGMRES: converge flag: 0
% PGMRES: iter:(1, 15)
%
% n: 128
% norm of difference of banded LU and PGMRES: 42.247549
% banded LU: max abs knot error: 0.681941
% PGMRES: max abs knot error: 0.000001
% PGMRES: converge flag: 0
% PGMRES: iter:(1, 17)
%
% order of convergence: (n: 8, 2n: 16): 2.103993
% order of convergence: (n: 16, 2n: 32): 2.043303
% order of convergence: (n: 32, 2n: 64): 2.017019
% order of convergence: (n: 64, 2n: 128): 2.009673