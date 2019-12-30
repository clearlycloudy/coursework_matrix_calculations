[u_gmres,flag,relres,iter] = solver_gmres_bv_2ndorder();

n = 8;
dim = n-1;
syms f1(x,y);
f1(x,y)=sin(x)*cos(y);
% f1(x,y)=x^(9/2)*y^(9/2);

reftemp = full(compute_grid(n,f1,x,y));
ref = reshape(reftemp', [dim*dim,1]);

max(max(abs(u_gmres-ref)))