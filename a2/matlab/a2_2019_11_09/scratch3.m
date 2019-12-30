clear;

n = 8;
dim = n-1;

h = zeros(dim,1);

%offset by 1 due to 1-indexing
%boundary conditions
gamma = sparse(n+1,n+1);
eta = sparse(n+1,n+1); %along y dimension
zeta = sparse(n+1,n+1); %along x dimension

%test with dummy values for boundary
%indexing for boundary: (x,y)
gamma(1:n+1,1)=1;
gamma(1:n+1,n+1)=1;
gamma(1,1:n+1)=1;
gamma(n+1,1:n+1)=1;
eta(1:n+1,1)=50;
eta(1:n+1,n+1)=50;
zeta(1,1:n+1)=50;
zeta(n+1,1:n+1)=50;

g = zeros(dim*dim,1);

domain_range = 1;
h = domain_range/n;

%dummy setting for testing
h = 1;

%set boundary values by going through each column
g2 = zeros(dim*dim,1);

off=1; %offset adjustment for boundary value access due to 1-indexing

for i = 1:2
    for j = 1:n-1
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
for i = n-2:n-1
    for j = 1:n-1
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

for i = 3:n-3
    for j = 1:2
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
    end
    for j = n-2:n-1
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

g2_debug_new = zeros(dim,dim);
for c = 1:dim*dim
    [i,j] = k_to_ij(dim,c);
    g2_debug_new(i,j) = g2(c);
end

g2_debug_new
