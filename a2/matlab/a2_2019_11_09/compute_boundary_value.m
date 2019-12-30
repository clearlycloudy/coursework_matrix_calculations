function boundary = compute_boundary_value(n,f,x,y)

%assumes f is a symbolic function taking in x and y as inputs

h = 1/n;
points = 0:1:n;
assert(size(points,2)==n+1);
domain = points .* h;

boundary = sparse(n+1,n+1); %offset due to 1-indexing

z = zeros(1,n+1);
o = ones(1,n+1);

boundary(1,:) = f(z, domain); %(x=0,y=[0:n]*h)
boundary(n+1,:) = f(o, domain); %(x=1,y=[0:n]*h)

%(x=[0:n]*h,y=0)
for i=2:n
    boundary(i,1) = f(domain(i),0);
end

%(x=[0:n]*h,y=1)
for i=2:n
    boundary(i,n+1) = f(domain(i),1);
end

