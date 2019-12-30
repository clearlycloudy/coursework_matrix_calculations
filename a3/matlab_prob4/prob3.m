clear

ns = [8 16 32 64]';
max_errs_1 = zeros(length(ns),1);
max_errs_2 = zeros(length(ns),1);

for i = 1:length(ns)
    n = ns(i)
    syms f1(x,y);
    f1(x,y)=sin(x)*cos(y);
    u1 = solve(n, f1, x, y);   
    reference_1 = full(compute_grid(n,f1,x,y));
    max_errs_1(i) = max(max(abs(u1-reference_1)));
    
    syms f2(x,y);
    f2(x,y)=x^(7/2)*y^(7/2);
    u2 = solve(n, f2, x, y);
    reference_2 = full(compute_grid(n,f2,x,y));
    max_errs_2(i) = max(max(abs(u2-reference_2)));
end

convergence_1 = zeros(length(ns)-1,1);
convergence_2 = zeros(length(ns)-1,1);

for i = 1:length(ns)-1
    convergence_1(i)=log2(max_errs_1(i)/max_errs_1(i+1));
    convergence_2(i)=log2(max_errs_2(i)/max_errs_2(i+1));
end

fprintf("sin(x)*cos(y):\n");
fprintf("max abs knot error\n");
max_errs_1
fprintf("order of convergence\n");
convergence_1

fprintf("x^(7/2)*y^(7/2):\n");
fprintf("max abs knot error\n");
max_errs_2
fprintf("order of convergence\n");
convergence_2