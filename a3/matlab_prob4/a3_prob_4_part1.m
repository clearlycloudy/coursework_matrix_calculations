n = 8;
dim = n-1;

syms f1(x,y);
%f1(x,y)=sin(x)*cos(y);
f1(x,y)=x^(9/2)*y^(9/2);

[u1, A, g] = solve_2(n, f1, x, y);
u = A\g;

u_cart = convert_to_cartesian(dim, u);

max(max(abs(u_cart-u1)));

reference_1 = full(compute_grid(n,f1,x,y));
max(max(abs(u_cart-reference_1)))