function grid = compute_grid(n,f,x,y)

% assumes f is a symbolic function taking in x and y as inputs
% return evaluations at interior grid points as f(grid(x,y))

h = 1/n;
points = 0:1:n;
assert(size(points,2)==n+1);
domain = points .* h;

temp = sparse(n+1,n+1); %offset due to 1-indexing

for i=1:n+1
    temp(:,i) = f(domain, ones(1,n+1).* domain(i));
end

%take only interior grid points: index range: [2:n,2:n]
%full(temp)
grid = temp(2:n,2:n);
