%assume domain is (0,1)
n=8;
dim = n-1;
h=1/n;

%Ax=b

%construct A
a = ones(dim,1)*-2;
b = ones(dim,1)*1;
B = spdiags([b a b], -1:1, dim, dim);
B=B./(h^2);
A=B^2;

syms w(x);
w(x) = x*(x-1)/2;
eval_w = eval(w([1:n-1]*h))';
eval_w_norm_inf = max(abs(eval_w));

y=B*eval_w;

norm(inv(B),inf)

%inf norm (eval_w) = 0.125
%inf norm (y) = 1
%inf norm (inv(B)) = 0.125

syms v(x);
v(x) = (x-1)*x*(2*(x-1)+1)/6;
v(x) = x*(x-1)*(x-2)*(x-3)/4;
v(x) = (x-1)*(x)*(x+1)*(x+2)/4;
eval_v = eval(v([1:n-1]*h))';
eval_v_norm_inf = max(abs(eval_v));
y2=A*eval_v;

norm(inv(A),inf)