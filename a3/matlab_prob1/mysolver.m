function [x,it] = mysolver(a, b, n, u, x)
    % assumes u(x) is a symbolic function
    
    range = b-a;
    h = range/n;
    t = 0:1:n-1;
    xs = h.*t + a;
    
    g_symbol = diff(u,x,4) + u;
    g = eval(g_symbol(xs)');
    
    %construct A
    r = [6 -4 1 zeros(1,n-5) 1 -4];
    A = sparse(toeplitz(r));
    A=A./(h^4); % for u''''
    I=spdiags([ones(n)],0:0,n,n); %for u
    A=A+I;
    A(:,1) = A(:,1)*2; %for periodic boundary condition
    
    %solve with CG
    tol = 10^(-8);
    maxit = 200;
    %x0 zero vector, preconditioner = I
    [x,flag,relres,iter] = pcg(A,g,tol,maxit,[],[],[]);
    it=iter;
    
end