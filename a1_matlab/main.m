clear
%run solver for each dim for each case and record max knot errors

dims = [ 8 16 32 64 128 ];

%note: reset errors before running
errors_1 = [];
errors_2 = [];
errors_3 = [];

%case 1: u(x)=x^3
for e = 1:length(dims)
    max_error_knots=0;
    n=dims(e);
    t=1:1:n-1;
    xs=1/n.*t(1:n-1);
    %g(x)=u''''(x)=0
    %[gamma_0, gamma_n] = u(x) {eval at x=0,1} = [0,1]
    %[zeta_0, zeta_n] = u''(x)=6x {eval at x=0,1} = [0,6]
    approx_1 = mysolver(n, @(x) 0, @()[0 1], @()[0 6]);
    exact_1 = arrayfun(@(x) x^3, xs)';
    max_error_knots = max(abs(approx_1-exact_1));
    errors_1 = [errors_1 max_error_knots];
    fprintf("case 1, n: %d, max knot error: %f\n", n, max_error_knots);
end

%case 2: u(x)=sin(x)
for e = 1:length(dims)
    max_error_knots=0;
    n=dims(e);
    t=1:1:n-1;
    xs=1/n.*t(1:n-1);
    %g(x)=u''''=sin(x)
    %[gamma_0, gamma_n] = u(x) {eval at x=0,1} = [sin(0),sin(1)]
    %[zeta_0, zeta_n] = u''(x)=-sin(x) {eval at x=0,1} = [-sin(0),-sin(1)]
    approx_2 = mysolver(n, @(x) sin(x), @()[0 sin(1)], @()[0 -sin(1)]);
    exact_2 = arrayfun(@(x) sin(x), xs)';
    max_error_knots = max(abs(approx_2-exact_2));
    errors_2 = [errors_2 max_error_knots];
    fprintf("case: 2, n: %d, max knot error: %f\n", n, max_error_knots);
end

%case 3: u(x)=x^(9/2)
for e = 1:length(dims)
    max_error_knots=0;
    n=dims(e);
    t=1:1:n-1;
    xs=1/n.*t(1:n-1);
    %g(x)=u''''=(9*7*5*3)/(2^4)*x^(1/2)
    %[gamma_0, gamma_n] = u(x) {eval at x=0,1} = [0,1]
    %[zeta_0, zeta_n] = u''(x)=(9*7)/(2^2)*x^(5/2) {eval at x=0,1} = [0,(9*7)/(2^2)]
    approx_3 = mysolver(n, @(x) (9*7*5*3)/(2^4)*x^(1/2), ...
        @()[0 1], @()[0 (9*7)/4]);
    exact_3 = arrayfun(@(x) x^(9/2), xs)';
    max_error_knots = max(abs(approx_3-exact_3));
    errors_3 = [errors_3 max_error_knots];
    fprintf("case: 3, n: %d, max knot error: %f\n", n, max_error_knots);
end

map_dim_to_error_index = containers.Map({8,16,32,64,128},{1,2,3,4,5});

pairs = [8,16; 16,32; 32,64; 64,128];

%calc order of convergence for case 2
for i = 1:size(pairs,1)
    n_a = pairs(i,1);
    n_b = pairs(i,2);
    y = calc_order_convergence(errors_1,errors_2,errors_3,...
        map_dim_to_error_index,...
        2, n_a, n_b);
    fprintf("case: 2, n: %d, n2: %d, order_convergence: %f\n", n_a, n_b, y);
end

%calc order of convergence for case 3
for i = 1:size(pairs,1)
    n_a = pairs(i,1);
    n_b = pairs(i,2);
    y = calc_order_convergence(errors_1,errors_2,errors_3,...
        map_dim_to_error_index,...
        3, n_a, n_b);
    fprintf("case: 3, n: %d, n2: %d, order_convergence: %f\n", n_a, n_b, y);
end