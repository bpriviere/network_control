function g = get_g( param, x, t)

    g = zeros( param.n, param.m);
    
    for i = 1:param.Nb
        row_idx = 2*param.Na*param.Nd + param.Nd + ...
            (i-1)*2*param.Nd + 1;
        col_idx = (i-1)*param.Nd + 1;
        g(row_idx:row_idx + param.Nd-1, ...
            col_idx:col_idx + param.Nd-1) = eye(param.Nd);
    end
    
    
end