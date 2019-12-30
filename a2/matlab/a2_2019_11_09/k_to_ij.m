function [i,j] = k_to_ij(dim,k)
    assert(k>=1);
    assert(k<=dim*dim);
    i=fix((k-1)/dim)+1;
    j=mod(k-1,dim)+1;
end