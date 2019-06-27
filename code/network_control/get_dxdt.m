function dxdt = get_dxdt( param, x, t)

    f = get_f( param, x, t);
    g = get_g( param, x, t);
    u = get_u( param, x, t);
    
    dxdt = f + g*u;
    
    fprintf('%d/%d\n', t, param.t(end));
end

