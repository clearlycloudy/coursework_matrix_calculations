clear
%run solver for each dim for each case and record max knot errors

dims = [ 8 16 32 64 128 ];

%note: reset errors before running
errors = [];
iterations = [];

b = 2*pi;
a = 0;
range = b-a;

%case: u(x)=sin(x)
for e = 1:length(dims)

    n=dims(e);
    
    h = range/n;
    t = 0:1:n-1;
    xs = range/n.*t + a;
    
    syms u(x);
    u(x)=sin(x);

    [approx,iter] = mysolver(a, b, n, u, x);
    exact = eval(u(xs)');

    max_error_knots = 0;
    max_error_knots = max(max_error_knots, max(abs(approx-exact)));
    errors = [errors max_error_knots];
    iterations = [iterations iter];
    
    fprintf("case 1, n: %d, max knot error: %f, iterations: %d\n", n, max_error_knots, iter);
end
