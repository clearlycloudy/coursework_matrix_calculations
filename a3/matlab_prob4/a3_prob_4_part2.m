n=8;
dim=n-1;
h=1/n;

syms f(x,y);
f(x,y)=x^(9/2)*y^(9/2);

%a u_xxxx + b u_xxyy + c u_yyyy + d u_xx + e u_yy + f u = g

syms a(x,y);
syms b(x,y);
syms c(x,y);
syms d(x,y);
syms e(x,y);
syms f(x,y);

a(x,y)=1+exp(x+y);
b(x,y)=1+1/(x+y);
c(x,y)=3+x*y;
d(x,y)=-sin(x)*sin(y);
e(x,y)=-1-exp(x+y);
f(x,y)=x+y;

I = eye(dim);

aa = zeros(dim*dim,1);
bb = zeros(dim*dim,1);
cc = zeros(dim*dim,1);
dd = zeros(dim*dim,1);
ee = zeros(dim*dim,1);
ff = zeros(dim*dim,1);
for i=1:n-1
    for j=1:n-1
        aa(ij_to_k(dim,i,j),1) = double(a(i/n*h,j/n*h));
        bb(ij_to_k(dim,i,j),1) = double(b(i/n*h,j/n*h));
        cc(ij_to_k(dim,i,j),1) = double(c(i/n*h,j/n*h));
        dd(ij_to_k(dim,i,j),1) = double(d(i/n*h,j/n*h));
        ee(ij_to_k(dim,i,j),1) = double(e(i/n*h,j/n*h));
        ff(ij_to_k(dim,i,j),1) = double(f(i/n*h,j/n*h));
    end
end

%u_xxxx block pentadiag [. . . . .]
T = spdiags([ones(dim,1) -4*ones(dim,1) 6*ones(dim,1) -4*ones(dim,1), ones(dim,1)],...
              -2:2,...
              dim,dim);
Txxxx = kron(T, I);

%u_yyyy pentadiag [.....]
T = spdiags([ones(dim,1) -4*ones(dim,1) 6*ones(dim,1) -4*ones(dim,1), ones(dim,1)],...
              -2:2,...
              dim,dim);
Tyyyy = kron(I, T);

%u_xxyy block tridiag [... ... ...]

%[      1     ]
%[   2 -8  2  ]
%[1 -8 20 -8 1]
%[   2 -8  2  ]
%[      1     ]

%[            ]
%[            ]
%[1 -4  6 -4 1]
%[            ]
%[            ]

%[      1     ]
%[     -4     ]
%[      6     ]
%[     -4     ]
%[      1     ]

%[      1     ]
%[     -4     ]
%[1 -4 12 -4 1]
%[     -4     ]
%[      1     ]

% 2 u_xxyy =
%[            ]
%[   2 -4  2  ]
%[  -4  8 -4  ]
%[   2 -4  2  ]
%[            ]

% u_xxyy =
%[            ]
%[   1 -2  1  ]
%[  -2  4 -2  ]
%[   1 -2  1  ]
%[            ]

T = spdiags([-2*ones(dim,1), 4*ones(dim,1), -2*ones(dim,1)],...
              -1:1,...
              dim,dim);
Txxyy_a = kron(I, T);

T = spdiags([ones(dim,1) -2*ones(dim,1) ones(dim,1)],...
              -1:1,...
              dim,dim);
          
II = spdiags([ones(dim,1)],...
              1,...
              dim,dim);
Txxyy_b = kron(II, T);

II = spdiags([ones(dim,1)],...
              -1,...
              dim,dim);
Txxyy_c = kron(II, T);

Txxyy = Txxyy_a + Txxyy_b + Txxyy_c;

%u_xx tridiag [. . .]
T = spdiags([ones(dim,1) -2*ones(dim,1) ones(dim,1)],...
              -1:1,...
              dim,dim);
Txx = kron(T, I);

%u_yy tridiag [...]
T = spdiags([ones(dim,1) -2*ones(dim,1) ones(dim,1)],...
              -1:1,...
              dim,dim);
Tyy = kron(I, T);

Tu = spdiags([ones(dim*dim,1)],...
              0,...
              dim*dim,dim*dim);
          
%----------

%multiply each row in Txxxx by scalar in aa
Txxxx_scaled = bsxfun(@times,Txxxx,aa);

%multiply each row in Txxyy by scalar in bb
Txxyy_scaled = bsxfun(@times,Txxyy,bb);

%multiply each row in Tyyyy by scalar in cc
Tyyyy_scaled = bsxfun(@times,Tyyyy,cc);

%[      1     ]
%[     -4     ]
%[      6     ]
%--------------
%[     -4     ] boundary
%--------------
%[      1     ] outside of boundary

full(Txxxx + 2*Txxyy + Tyyyy)

%multiply each row in Txx by scalar in dd
Txx_scaled = bsxfun(@times,Txx,dd);

%multiply each row in Tyy by scalar in ee
Tyy_scaled = bsxfun(@times,Tyy,ee);

%expand inner grid matrix by 2 at boundary of each direction

%o o o o o o o o o%
%o o o o o o o o o%
%o o . . . . . o o%
%o o . . . . . o o%
%o o . . . . . o o%
%o o . . . . . o o%
%o o . . . . . o o%
%o o o o o o o o o%
%o o o o o o o o o%

%resulting in each row of A matrix:
%o o o o o o o o o | o o o o o o o o o | o o . . . . . o o | ... | o o o o o o o o o | o o o o o o o o o

%[            ]
%[      1ee   ]
%[     -2ee   ]
%--------------
%[      1ee   ] boundary
%--------------
%[            ] outside of boundary
%dot product weighted coefficents at boundary and 
%with gamma(x,y) and move to right hand side via negation
%finally scale by 1/h^m if necessary, m=0/2/4
%do the same for Txx, Tyyyy, Txxxx, Txxyy

%take care of points on outside of boundary:
%for Tyyyy:
%obtain value of point at outside of boundary, k
%apply 2nd order boundary condition stencil for Txxxx by scaling with k and
%cancelling stencil values, move zeta values scaled by 1/h^2 to right hand
%side via negation
%dot product weighted coefficents at boundary and 
%with gamma(x,y) and move to right hand side via negation scaled by 1/h^m
%scale extract inner grid matrix by and scale by 1/h^m
%
%apply 1st order boundary condition stencil for Tyyyy by scaling with k and
%cancelling stencil values, move eta values scaled by  1/h^3 to right hand
%side vai negation
%dot product weighted coefficents at boundary and 
%with gamma(x,y) and move to right hand side via negation scaled by 1/h^m
%scale extract inner grid matrix by and scale by 1/h^m

%multiply each row in Tu by scalar in ff
Tu_scaled = bsxfun(@times,Tu,ff);

A =   1/(h^4).*(Txxxx_scaled + Txxyy_scaled + Tyyyy_scaled)...
    + 1/(h^2).*(Txx_scaled + Tyy_scaled)...
    + Tu;

syms u(x,y)
u(x,y)=x^(9/2)*y^(9/2);

u_xxxx = diff(u,x,4);
u_yyyy = diff(u,y,4);
u_xxyy = diff(diff(u,x,2),y,2);
u_xx = diff(u,x,2);
u_yy = diff(u,y,2);

[AA,gg] = create_A_g_generic(n,x,y,a,b,c,d,e,f,u_xxxx,u_xxyy,u_yyyy,u_xx,u_yy);
