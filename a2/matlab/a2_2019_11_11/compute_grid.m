function grid = compute_grid(n,f,x,y)
% assumes f is a symbolic function taking in x and y as inputs
% return evaluations at interior grid points as f(grid(x,y))

h = 1/n;

points = 1:1:n-1;
assert(size(points,2)==n-1);

domain = points .* h;

grid = zeros(n-1,n-1);

for i=1:1:n-1
    grid(:,i) = f(domain, ones(1,n-1).* domain(i));
end

