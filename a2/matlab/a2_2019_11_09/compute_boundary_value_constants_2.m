function g2 = compute_boundary_value_constants_2(n,gamma,eta,zeta)
%compute boundary constants that touches or goes out of grid for the right hand side of Au=g
%returns negated value to be used for right side

%assumes gamma, eta, zeta are sparse and of size (n+1,n+1), offset by 1 due to 1-indexing
%gamma: 0th order boundary condition
%eta: 2nd order boundary condition along y dimension
%zeta: 2nd order boundary condition along x dimension

dim = n-1;
h = 1/n;

%dummy setting for testing
%h = 1;

g2 = zeros(dim*dim,1);

off=1; %offset adjustment for boundary value access due to 1-indexing

for i = 1:n-1
    for j = 1:n-1

        if ( (3<=i) && (i<= n-3) && (3<=j) && (j<= n-3))
            continue;
        end
        
        if (i==1)
            g2(ij_to_k(dim,i,j)) = g2(ij_to_k(dim,i,j))...
                                    -1/(h^4)*(+2*gamma(0+off,j-1+off)...
                                              -6*gamma(0+off,j+off)...
                                              +2*gamma(0+off,j+1+off))...
                                    -1/(h^2)*(zeta(0+off,j+off));
        elseif (i==2)
            g2(ij_to_k(dim,i,j)) = g2(ij_to_k(dim,i,j))...
                                    -1/(h^4)*(gamma(0+off,j+off));
        end
        
        if (i==n-1)
            g2(ij_to_k(dim,i,j)) = g2(ij_to_k(dim,i,j))...
                                    -1/(h^4)*(+2*gamma(n+off,j-1+off)...
                                              -6*gamma(n+off,j+off)...
                                              +2*gamma(n+off,j+1+off))...
                                    -1/(h^2)*(zeta(n+off,j+off));
        elseif (i==n-2)
            g2(ij_to_k(dim,i,j))= g2(ij_to_k(dim,i,j))...
                                    -1/(h^4)*(gamma(n+off,j+off));
        end
        
        if (j==1)
            g2(ij_to_k(dim,i,j)) = g2(ij_to_k(dim,i,j))...
                                    -1/(h^4)*(+2*gamma(i-1+off,0+off)...
                                              -6*gamma(i+off,0+off)...
                                              +2*gamma(i+1+off,0+off))...
                                    -1/(h^2)*(eta(i+off,0+off));
        elseif (j==2)
            g2(ij_to_k(dim,i,j)) = g2(ij_to_k(dim,i,j))...
                                    -1/(h^4)*(gamma(i+off,0+off));
        end

        if (j==n-1)
            g2(ij_to_k(dim,i,j)) = g2(ij_to_k(dim,i,j))...
                                    -1/(h^4)*(+2*gamma(i-1+off,n+off)...
                                              -6*gamma(i+off,n+off)...
                                              +2*gamma(i+1+off,n+off))...
                                    -1/(h^2)*(eta(i+off,n+off));
        elseif (j==n-2)
            g2(ij_to_k(dim,i,j)) = g2(ij_to_k(dim,i,j))...
                                    -1/(h^4)*(gamma(i+off,n+off));
        end
    end
end

%take care of overlap in corners in previous loops
g2(ij_to_k(dim,1,1)) = g2(ij_to_k(dim,1,1)) + 1/(h^4)*(2*gamma(0+off,0+off));
g2(ij_to_k(dim,1,n-1)) = g2(ij_to_k(dim,1,n-1)) + 1/(h^4)*(2*gamma(0+off,n+off));
g2(ij_to_k(dim,n-1,1)) = g2(ij_to_k(dim,n-1,1)) + 1/(h^4)*(2*gamma(n+off,0+off));
g2(ij_to_k(dim,n-1,n-1)) = g2(ij_to_k(dim,n-1,n-1)) + 1/(h^4)*(2*gamma(n+off,n+off));
