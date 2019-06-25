function f = get_f( param, x, t)

    f = zeros(param.n,1);

    if strcmp( param.model, 'reynolds')

        A = get_A(param, x);
        for i = 1:param.Na
            p_i = get_p( param, x, i);
            v_i = get_v( param, x, i);
            
            a_i = zeros(param.Nd,1);
            for j = 1:param.N
                if ~ (i == j)
                    p_j = get_p( param, x, j);
                    v_j = get_v( param, x, j);

                    r_ij = p_j - p_i;
                    a_i = a_i + A(i,j)*(...
                        param.kv*(v_j - v_i) + ...
                        param.kx*r_ij*(1 - param.R_des/norm(r_ij)^3));
                end
            end
            
            idx = (i-1)*2*param.Nd + 1;
            f( idx: idx + 2*param.Nd - 1) = [ v_i; a_i];
        end
        
        for i = param.Na+1: param.N
            v_i = get_v( param, x, i);
            idx = (i-1)*2*param.Nd + 1;
            f( idx: idx + 2*param.Nd - 1) = [ v_i; zeros(param.Nd,1)];
        end
        
    end
end
