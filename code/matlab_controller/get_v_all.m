function [vel_a, vel_b] = get_v_all( param, x)

    vel_a = nan(param.Nd*param.Na,1);
    vel_b = nan(param.Nd*param.Nb,1);

    for i = 1:param.Na
        idx = (i-1)*param.Nd + 1;
        vel_a(idx:idx+param.Nd-1) = get_v( param, x, i);        
    end

    for i = 1:param.Nb
        idx = (i-1)*param.Nd + 1;
        vel_b(idx:idx+param.Nd-1) = get_v( param, x, i+param.Na);  
    end
end