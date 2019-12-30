function y = calc_order_convergence(errors_1, errors_2, errors_3, map_dim_to_error_index, whichcase, n1, n2)
    if whichcase==1
        e1=errors_1(map_dim_to_error_index(n1));
        e2=errors_1(map_dim_to_error_index(n2));
        order_convergence=log2(e1/e2);
    elseif whichcase==2
        e1=errors_2(map_dim_to_error_index(n1));
        e2=errors_2(map_dim_to_error_index(n2));
        order_convergence=log2(e1/e2);
    else
        e1=errors_3(map_dim_to_error_index(n1));
        e2=errors_3(map_dim_to_error_index(n2));
        order_convergence=log2(e1/e2);
    end
    y = order_convergence;
end