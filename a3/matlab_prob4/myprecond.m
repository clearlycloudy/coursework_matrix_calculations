function u_solve = myprecond()
    % n := number of intervals in grid
    % assumes f is a symbolic expression in terms of x and y
    % returns approximation of u on interior of grid in cartesian coordinates
    n=8;
    h = 1/n;
    dim = n-1;
%     [AA, B, T, C] = create_A_2(n); %A is actually not needed
    [AA, B, T] = create_A(n);
    u_solve=AA;
end