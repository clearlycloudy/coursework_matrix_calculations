function out = convert_to_cartesian(dim, input)
%convert to cartesian coordinate
out = zeros(dim,dim);
for c = 1:dim*dim
    [i,j] = k_to_ij(dim,c);
    out(i,j) = input(c);
end