function [A, rhs] = create_A_g_generic(n,x,y,a,b,c,d,e,f,g,u_xxxx,u_xxyy,u_yyyy,u_xx,u_yy,eta,zeta,gamma)
    %assumes a,b,c,d,e,f,u_xxxx,u_xxyy,u_yyyy,u_xx,u_yy are symbolic functions in terms of x and y
    %a u_xxxx + b u_xxyy + c u_yyyy + d u_xx + e u_yy + f u = g
    
    h=1/n;
    dim = n-1;
    
    %padd 4 in each direction so that values at boundary and 1 grid outside
    %of boundary are captured
    pad = 4;
    pad_offset = 2;
    dim_mod = dim+4;
    
    rhs = sparse(zeros(dim_mod*dim_mod,1));

    gamma_aug = zeros(dim_mod,dim_mod);
    
    aa = zeros(dim,dim);
    bb = zeros(dim,dim);
    cc = zeros(dim,dim);
    dd = zeros(dim,dim);
    ee = zeros(dim,dim);
    ff = zeros(dim,dim);
    gg = zeros(dim,dim);
    
    for i=1:n-1
        for j=1:n-1
            aa(i,j) = double(a(i*h,j*h));
            bb(i,j) = double(b(i*h,j*h));
            cc(i,j) = double(c(i*h,j*h));
            dd(i,j) = double(d(i*h,j*h));
            ee(i,j) = double(e(i*h,j*h));
            ff(i,j) = double(f(i*h,j*h));
            gg(i,j) = double(g(i*h,j*h));
        end
    end
    
    aa_aug = zeros(dim_mod,dim_mod);
    bb_aug = zeros(dim_mod,dim_mod);
    cc_aug = zeros(dim_mod,dim_mod);
    dd_aug = zeros(dim_mod,dim_mod);
    ee_aug = zeros(dim_mod,dim_mod);
    ff_aug = zeros(dim_mod,dim_mod);
    gg_aug = zeros(dim_mod,dim_mod);
    
    for i=1:n-1
       aa_aug(:,2+i) = [ 0; 0; aa(:,i); 0; 0];
       bb_aug(:,2+i) = [ 0; 0; bb(:,i); 0; 0];
       cc_aug(:,2+i) = [ 0; 0; cc(:,i); 0; 0];
       dd_aug(:,2+i) = [ 0; 0; dd(:,i); 0; 0];
       ee_aug(:,2+i) = [ 0; 0; ee(:,i); 0; 0];
       ff_aug(:,2+i) = [ 0; 0; ff(:,i); 0; 0];
       gg_aug(:,2+i) = [ 0; 0; gg(:,i); 0; 0];
    end
    
    eta_eval = zeros(dim_mod,dim_mod);
    zeta_eval = zeros(dim_mod,dim_mod);
    
    for i=1:n-1+2
        point = (i-1)*h;
        eta_eval(i+1,2) = eta(point,0);
        eta_eval(i+1,dim_mod-1) = eta(point,1);
        
        zeta_eval(2,i+1) = zeta(0,point);
        zeta_eval(dim_mod-1,i+1) = zeta(1,point);
        
        gamma_aug(i+1,2) = gamma(point,0);
        gamma_aug(i+1,dim_mod-1) = gamma(point,1);
        gamma_aug(2,i+1) = gamma(0,point);
        gamma_aug(dim_mod-1,i+1) = gamma(1,point);
    end
    
    %stack by column (y-direction increasing fastest)
    aa_stack = sparse(reshape(aa_aug',[dim_mod*dim_mod,1]));
    bb_stack = sparse(reshape(bb_aug',[dim_mod*dim_mod,1]));
    cc_stack = sparse(reshape(cc_aug',[dim_mod*dim_mod,1]));
    dd_stack = sparse(reshape(dd_aug',[dim_mod*dim_mod,1]));
    ee_stack = sparse(reshape(ee_aug',[dim_mod*dim_mod,1]));
    ff_stack = sparse(reshape(ff_aug',[dim_mod*dim_mod,1]));
    gg_stack = sparse(reshape(gg_aug',[dim_mod*dim_mod,1]));
    
    eta_stack = sparse(reshape(eta_eval',[dim_mod*dim_mod,1]));
    zeta_stack = sparse(reshape(zeta_eval',[dim_mod*dim_mod,1]));
    
    gamma_stack = sparse(reshape(gamma_aug',[dim_mod*dim_mod,1]));
    
    op = ones(dim_mod,1);
    op_mod = ones(dim_mod,1);
    op_mod(1:2,1)=0;
    op_mod(dim_mod-1:dim_mod,1)=0;
    
    %u_xxxx block pentadiag [. . . . .]
    T = spdiags([op -4*op 6*op -4*op, op],...
                  -2:2,...
                  dim_mod,dim_mod);
    T(1:2,:)=0;
    T(dim_mod-1:dim_mod,:)=0;
    
    Txxxx = kron(T, diag(op_mod));

    %u_yyyy pentadiag [.....]
    T = spdiags([op -4*op 6*op -4*op, op],...
                  -2:2,...
                  dim_mod,dim_mod);
    T(1:2,:)=0;
    T(dim_mod-1:dim_mod,:)=0;
    
    Tyyyy = kron(diag(op_mod), T);

    %u_xxyy:
    T = spdiags([-2*op, 4*op, -2*op],...
                  -1:1,...
                  dim_mod,dim_mod);
    T(1:2,:)=0;
    T(dim_mod-1:dim_mod,:)=0;
    Txxyy_a = kron(diag(op_mod), T); %same column

    T = spdiags([op -2*op op],...
                  -1:1,...
                  dim_mod,dim_mod);
    T(1:2,:)=0;
    T(dim_mod-1:dim_mod,:)=0;

    II = spdiags([op],...
                  1,...
                  dim_mod,dim_mod);
    II(1:2,:)=0;
    II(end-1:end,:)=0;
    Txxyy_b = kron(II, T); %column to the right

    II = spdiags([op],...
                  -1,...
                  dim_mod,dim_mod);
    II(1:2,:)=0;
    II(end-1:end,:)=0;
    Txxyy_c = kron(II, T); %column to the left

    Txxyy = Txxyy_a + Txxyy_b + Txxyy_c;
    
    %u_xx tridiag [. . .]
    T = spdiags([op -2*op op],...
                  -1:1,...
                  dim_mod,dim_mod);
    T(1:2,:)=0;
    T(dim_mod-1:dim_mod,:)=0;
    Txx = kron(T, diag(op_mod));

    %u_yy tridiag [...]
    T = spdiags([op -2*op op],...
                  -1:1,...
                  dim_mod,dim_mod);
    T(1:2,:)=0;
    T(dim_mod-1:dim_mod,:)=0;
    Tyy = kron(diag(op_mod), T);

    %u
    T = spdiags([op],...
                 0,...
                 dim_mod,dim_mod);
    T(1:2,:)=0;
    T(dim_mod-1:dim_mod,:)=0;
    Tu = kron(diag(op_mod),T);

    %multiply by provided coefficients a,b,c,d,e,f
    %a u_xxxx + b u_xxyy + c u_yyyy + d u_xx + e u_yy + f u = g
    Txxxx_scaled = bsxfun(@times,Txxxx,aa_stack);
    Txxyy_scaled = bsxfun(@times,Txxyy,bb_stack);
    Tyyyy_scaled = bsxfun(@times,Tyyyy,cc_stack);
    Txx_scaled = bsxfun(@times,Txx,dd_stack);
    Tyy_scaled = bsxfun(@times,Tyy,ee_stack);
    Tu_scaled = bsxfun(@times,Tu,ff_stack);
    
    %gamma_stack contains the values at boundary, dot product with each row
    %to get boundary constants, move to right hand side of equation
    rhs = rhs - 1/(h^4)*(Txxyy_scaled * gamma_stack);
    rhs = rhs - 1/(h^2)*(Txx_scaled * gamma_stack);
    rhs = rhs - 1/(h^2)*(Tyy_scaled * gamma_stack);
    rhs = rhs - (Tu_scaled * gamma_stack);
    
    %cancel value outside of grid for boundary condition in x-direction:
    %Txxxx_scaled(:,1:dim_mod)
    %Txxxx_scaled(:,dim_mod*dim_mod-dim_mod+1:dim_mod*dim_mod)
    Tboundary_x = [1 -2 1];
    for row=1:dim_mod*dim_mod
        for i=1:dim_mod
            val = Txxxx_scaled(row,i);
            if val ~= 0
                val_negate = -val;
                scale = val_negate/Tboundary_x(1);
                rhs(row) = rhs(row) + scale * 1/(h^2) * zeta_stack(i+dim_mod);
                Tboundary_x_scaled = scale.* Tboundary_x;
                for j=1:length(Tboundary_x_scaled)
                    Txxxx_scaled(row,i+(j-1)*dim_mod) = Txxxx_scaled(row,i+(j-1)*dim_mod) + Tboundary_x_scaled(j);
                end
                break;
            end
        end
        for i=dim_mod*dim_mod:-1:dim_mod*dim_mod-dim_mod+1
            val = Txxxx_scaled(row,i);
            if val ~= 0
                val_negate = -val;
                scale = val_negate/Tboundary_x(3);
                rhs(row) = rhs(row) + scale * 1/(h^2) * zeta_stack(i-dim_mod);
                Tboundary_x_scaled = scale.* Tboundary_x;
                for j=1:length(Tboundary_x_scaled)
                    Txxxx_scaled(row,i-(j-1)*dim_mod) = Txxxx_scaled(row,i-(j-1)*dim_mod) + Tboundary_x_scaled(3-j+1);
                end
                break;
            end
        end
    end
    %sanity check that all exterior grid points are zeroed out
    assert(sum(sum(Txxxx_scaled(:,1:dim_mod)~=0))==0);
    assert(sum(sum(Txxxx_scaled(:,dim_mod*dim_mod-dim_mod+1:dim_mod*dim_mod)~=0))==0);
    
    %gamma_stack contains the values at boundary, dot product with each row
    %to get boundary constants, move to right hand side of equation
    rhs = rhs - 1/(h^4)*(Txxxx_scaled * gamma_stack);
    
    %cancel value outside of grid for boundary condition in y-direction:
    %col_indices = 1:dim_mod:dim_mod*dim_mod;
    %Tyyyy_scaled(:,col_indices)
    %Tyyyy_scaled(:,dim_mod:dim_mod:dim_mod*dim_mod)
    %bottom: eta = 1/h*[-1/4 -5/6 3/2 -1/2 5/60]
    %top: eta = 1/h*flip([1/4 5/6 -3/2 1/2 -5/60])
    Tboundary_y_bottom = [-1/4 -5/6 3/2 -1/2 5/60];
    Tboundary_y_top = flip([1/4 5/6 -3/2 1/2 -5/60]);
    for row=1:dim_mod*dim_mod
        for i=1:dim_mod:dim_mod*dim_mod
            val = Tyyyy_scaled(row,i);
            if val ~= 0
                val_negate = -val;
                scale = val_negate/Tboundary_y_bottom(1);
                rhs(row) = rhs(row) + scale * 1/(h^3) * eta_stack(i+1); %apply scale * 1/h^3 * eta to rhs of equation
                Tboundary_y_scaled = scale.* Tboundary_y_bottom;
                %apply BC stencil values
                for j=1:length(Tboundary_y_scaled)
                    Tyyyy_scaled(row,i+j-1) = Tyyyy_scaled(row,i+j-1) + Tboundary_y_scaled(j);
                end
                break; %apply only once when first non-zero exterior grid point is detected
            end
        end
        for i=flip(dim_mod:dim_mod:dim_mod*dim_mod)
            val = Tyyyy_scaled(row,i);
            if val ~= 0
                val_negate = -val;
                scale = val_negate/Tboundary_y_top(end);
                rhs(row) = rhs(row) + scale * 1/(h^3) * eta_stack(i-1); %apply scale * 1/h^3 * eta to rhs of equation
                Tboundary_y_scaled = scale.* Tboundary_y_top;
                %apply BC stencil values
                for j=1:length(Tboundary_y_scaled)
                    Tyyyy_scaled(row,i-(j-1)) = Tyyyy_scaled(row,i-(j-1)) + Tboundary_y_scaled(length(Tboundary_y_scaled)-j+1);
                end
                break; %apply only once when first non-zero exterior grid point is detected
            end
        end
    end
    %sanity check that all exterior grid points are zeroed out
    assert(sum(sum(Tyyyy_scaled(:,1:dim_mod:dim_mod*dim_mod)~=0))==0);
    assert(sum(sum(Tyyyy_scaled(:,dim_mod:dim_mod:dim_mod*dim_mod)~=0))==0);
    
    %gamma_stack contains the values at boundary, dot product with each row
    %to get boundary constants, move to right hand side of equation
    rhs = rhs - 1/(h^4)*(Tyyyy_scaled * gamma_stack);
    
    rhs = rhs + gg_stack;
    
    Atemp = 1/(h^4)*(Txxxx_scaled + Txxyy_scaled + Tyyyy_scaled) +...
        1/(h^2)*(Txx_scaled + Tyy_scaled) + Tu_scaled;
    
    %get only interior points
    Atemp = Atemp(2*dim_mod+1:dim_mod*dim_mod-dim_mod*2,...
                  2*dim_mod+1:dim_mod*dim_mod-dim_mod*2);
    A = sparse(zeros(dim*dim,dim*dim));
    
    for i =1:dim
         offset_i = (i-1) * dim_mod + 1;
         for j =1:dim
             offset_j = (j-1) * dim_mod + 1;
             A((i-1)*dim+1:(i-1)*dim+1+dim-1,...
               (j-1)*dim+1:(j-1)*dim+1+dim-1) = Atemp(offset_i+2:offset_i+2+dim-1,...
                                                      offset_j+2:offset_j+2+dim-1);
            
         end
    end
    
    rhs_m_by_n = reshape(rhs,[dim_mod,dim_mod]);
    rhs_temp = rhs_m_by_n(3:3+dim-1, 3:3+dim-1);
    rhs = reshape(rhs_temp,[dim*dim,1]);

end
    
