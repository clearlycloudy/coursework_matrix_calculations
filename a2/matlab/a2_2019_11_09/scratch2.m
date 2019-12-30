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

%(1,1):
    g2(ij_to_k(dim,1,1))=-1/(h^4)*( 2*gamma(0+off,0+off)-6*gamma(0+off,1+off)+2*gamma(0+off,2+off)...
                                   -6*gamma(1+off,0+off)+2*gamma(2+off,0+off))...
                         -1/(h^2)*(eta(1+off,0+off) + zeta(0+off,1+off));
%(1,2):
    g2(ij_to_k(dim,1,2))=-1/(h^4)*( 2*gamma(0+off,1+off)-6*gamma(0+off,2+off)+2*gamma(0+off,3+off)...
                                     +gamma(1+off,0+off))...
                         -1/(h^2)*(zeta(0+off,2+off));
%(1,3:n-3):
for j = 3:1:n-3
    g2(ij_to_k(dim,1,j))=-1/(h^4)*(2*gamma(0+off,j-1+off)-6*gamma(0+off,j+off)+2*gamma(0+off,j+1+off))...
                         -1/(h^2)*(zeta(0+off,j+off));
end
%(1,n-2):
    g2(ij_to_k(dim,1,n-2))=-1/(h^4)*(2*gamma(0+off,n-3+off)-6*gamma(0+off,n-2+off)+2*gamma(0+off,n-1+off)...
                                      +gamma(1+off,n+off))...
                           -1/(h^2)*(zeta(0+off,n-2+off));
%(1,n-1):
    g2(ij_to_k(dim,1,n-1))=-1/(h^4)*( 2*gamma(0+off,n-2+off)-6*gamma(0+off,n-1+off)+2*gamma(0+off,n+off)...
                                     -6*gamma(1+off,n+off)+2*gamma(2+off,n+off))...
                           -1/(h^2)*(eta(1+off,n+off) + zeta(0+off,n-1+off));

%(2,1):
g2(ij_to_k(dim,2,1))=-1/(h^4)*(2*gamma(1+off,0+off)-6*gamma(2+off,0+off)+2*gamma(3+off,0+off)...
                                +gamma(0+off,1+off))...
                     -1/(h^2)*(eta(2+off,0+off));
%(2,2):
g2(ij_to_k(dim,2,2))=-1/(h^4)*(gamma(2+off,0+off)+gamma(0+off,2+off));
%(2,3:n-3):
for j = 3:1:n-3
    g2(ij_to_k(dim,2,j))=-1/(h^4)*(gamma(0+off,j+off));
end
%(2,n-2):
g2(ij_to_k(dim,2,n-2))=-1/(h^4)*(gamma(0+off,n-2+off)+gamma(2+off,n+off));
%(2,n-1):
g2(ij_to_k(dim,2,n-1))=-1/(h^4)*(  gamma(0+off,n-1+off)...
                                +2*gamma(1+off,n+off)-6*gamma(2+off,n+off)+2*gamma(3+off,n+off))...
                       -1/(h^2)*(eta(2+off,n+off));

%(3:n-3, :):
for i = 3:1:n-3
    g2(ij_to_k(dim,i,1))=-1/(h^4)*(2*gamma(i-1+off,0+off)-6*gamma(i+off,0+off)+2*gamma(i+1+off,0+off))...
                         -1/(h^2)*(eta(i+off,0+off));
    g2(ij_to_k(dim,i,2))=-1/(h^4)*2*gamma(i+off,0+off);
    g2(ij_to_k(dim,i,n-2))=-1/(h^4)*2*gamma(i+off,n+off);
    g2(ij_to_k(dim,i,n-1))=-1/(h^4)*(2*gamma(i-1+off,n+off)-6*gamma(i+off,n+off)+2*gamma(i+1+off,n+off))...
                           -1/(h^2)*(eta(i+off,n+off));
end

%(n-2,1):
g2(ij_to_k(dim,n-2,1))=-1/(h^4)*(2*gamma(n-3+off,0+off)-6*gamma(n-2+off,0+off)+2*gamma(n-1+off,0+off)...
                                  +gamma(n+off,1+off))...
                       -1/(h^2)*(eta(n-2+off,0+off));
%(n-2,2):
g2(ij_to_k(dim,n-2,2))=-1/(h^4)*(gamma(n-2+off,0+off)+gamma(n+off,2+off));
%(n-2,3:n-3):
for j = 3:1:n-3
    g2(ij_to_k(dim,n-2,j))=-1/(h^4)*(gamma(n+off,j+off));
end
%(n-2,n-2):
g2(ij_to_k(dim,n-2,n-2))=-1/(h^4)*(gamma(n+off,n-2+off)+gamma(n-2+off,n+off));
%(n-2,n-1):
g2(ij_to_k(dim,n-2,n-1))=-1/(h^4)*(   gamma(n+off,n-1+off)...
                                   +2*gamma(n-3+off,n+off)-6*gamma(n-2+off,n+off)+2*gamma(n-1+off,n+off))...
                         -1/(h^2)*(eta(n-2+off,n+off));

%i=n-1:
%(n-1,1):
    g2(ij_to_k(dim,n-1,1))=-1/(h^4)*( 2*gamma(n+off,0+off)-6*gamma(n+off,1+off)+2*gamma(n+off,2+off)...
                                     +2*gamma(n-2+off,0+off)-6*gamma(n-1+off,0+off))...
                           -1/(h^2)*(eta(n-1+off,0+off) + zeta(n+off,1+off));
%(n-1,2):
    g2(ij_to_k(dim,n-1,2))=-1/(h^4)*(2*gamma(n+off,1+off)-6*gamma(n+off,2+off)+2*gamma(n+off,3+off)...
                                      +gamma(n-1+off,0+off))...
                           -1/(h^2)*(zeta(n+off,2+off));
%(n-1,3:n-3):
for j = 3:1:n-3
    g2(ij_to_k(dim,n-1,j))=-1/(h^4)*(2*gamma(n+off,j-1+off)-6*gamma(n+off,j+off)+2*gamma(n+off,j+1+off))...
                           -1/(h^2)*(zeta(n+off,j+off));
end
%(n-1,n-2):
    g2(ij_to_k(dim,n-1,n-2))=-1/(h^4)*(2*gamma(n+off,n-3+off)-6*gamma(n+off,n-2+off)+2*gamma(n+off,n-1+off)...
                                        +gamma(n-1+off,n+off))...
                             -1/(h^2)*(zeta(n+off,n-2+off));
%(n-1,n-1):
    g2(ij_to_k(dim,n-1,n-1))=-1/(h^4)*( 2*gamma(n+off,n-2+off)-6*gamma(n+off,n-1+off)+2*gamma(n+off,n+off)...
                                       +2*gamma(n-2+off,n+off)-6*gamma(n-1+off,n+off))...
                             -1/(h^2)*(eta(n-1+off,n+off) + zeta(n+off,n-1+off));

g2_debug = zeros(dim,dim);
for c = 1:dim*dim
    [i,j] = k_to_ij(dim,c);
    g2_debug(i,j) = g2(c);
end

g2_debug
