function g2 = compute_boundary_value_constants_2(n,gamma,eta,zeta)
%compute boundary constants that touches or goes out of grid for the right hand side of Au=g
%returns negated value to be used for right side

%assumes gamma, eta, zeta are symbolic functions
%gamma: 0th order boundary condition
%eta: 1st order boundary condition along y dimension
%zeta: 2nd order boundary condition along x dimension

dim = n-1;
h = 1/n;

g2 = zeros(dim*dim,1);

for i = 1:n-1
    for j = 1:n-1

        if ( (3<=i) && (i<= n-3) && (3<=j) && (j<= n-3))
            continue;
        end
        
        if (i==1)
            g2(ij_to_k(dim,i,j)) = g2(ij_to_k(dim,i,j))...
                                    -1/(h^4)*(+2*gamma(0,(j-1)*h)...
                                              -6*gamma(0,j*h)...
                                              +2*gamma(0,(j+1)*h))...
                                    -1/(h^2)*(zeta(0,j*h));
        elseif (i==2)
            g2(ij_to_k(dim,i,j)) = g2(ij_to_k(dim,i,j))...
                                    -1/(h^4)*(gamma(0,j*h));
        end
        
        if (i==n-1)
            g2(ij_to_k(dim,i,j)) = g2(ij_to_k(dim,i,j))...
                                    -1/(h^4)*(+2*gamma(1,(j-1)*h)...
                                              -6*gamma(1,j*h)...
                                              +2*gamma(1,(j+1)*h))...
                                    -1/(h^2)*(zeta(1,j*h));
        elseif (i==n-2)
            g2(ij_to_k(dim,i,j))= g2(ij_to_k(dim,i,j))...
                                    -1/(h^4)*(gamma(1,j*h));
        end
        
        if (j==1)
            %apply o(h^4) accuracy u_y BC: 
            % 4*eta/h^3 = +4*[-1/4 -5/6 3/2 -1/2 -5/60]*(1/h^4), 
            % where 4*-5/6=-10/3 is the boundary entry
            g2(ij_to_k(dim,i,j)) = g2(ij_to_k(dim,i,j))...
                                    -1/(h^4)*(+2*gamma((i-1)*h,0)...
                                              +(-8-10/3)*gamma(i*h,0)...
                                              +2*gamma((i+1)*h,0))...
                                    +4/(h^3)*(eta(i*h,0));
        elseif (j==2)
            g2(ij_to_k(dim,i,j)) = g2(ij_to_k(dim,i,j))...
                                    -1/(h^4)*(gamma(i*h,0));
        end

        if (j==n-1)
            %apply o(h^4) accuracy u_y BC:
            % -4*eta/h^3 = -4*flip([1/4 5/6 -3/2 1/2 -5/60])*(1/h^4), 
            % where -4*5/6=-10/3 is the boundary entry
            g2(ij_to_k(dim,i,j)) = g2(ij_to_k(dim,i,j))...
                                    -1/(h^4)*(+2*gamma((i-1)*h,1)...
                                              +(-8-10/3)*gamma(i*h,1)...
                                              +2*gamma((i+1)*h,1))...
                                    -4/(h^3)*(eta(i*h,1));
        elseif (j==n-2)
            g2(ij_to_k(dim,i,j)) = g2(ij_to_k(dim,i,j))...
                                    -1/(h^4)*(gamma(i*h,1));
        end
    end
end

%take care of overlap in corners in previous loops
g2(ij_to_k(dim,1,1)) = g2(ij_to_k(dim,1,1)) + 1/(h^4)*(2*gamma(0,0));
g2(ij_to_k(dim,1,n-1)) = g2(ij_to_k(dim,1,n-1)) + 1/(h^4)*(2*gamma(0,1));
g2(ij_to_k(dim,n-1,1)) = g2(ij_to_k(dim,n-1,1)) + 1/(h^4)*(2*gamma(1,0));
g2(ij_to_k(dim,n-1,n-1)) = g2(ij_to_k(dim,n-1,n-1)) + 1/(h^4)*(2*gamma(1,1));
