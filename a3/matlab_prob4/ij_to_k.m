function k = ij_to_k(dim,i,j)
    assert(i>=1);
    assert(i<=dim);
    assert(j>=1);
    assert(j<=dim);
    k=dim*(i-1)+j;
end