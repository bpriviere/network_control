function vdot_a = get_vdot_all(param,x,t)

    vdot_a_idx = nan(param.Nd*param.Na,1);
    for i = 1:param.Na
        idx = (i-1)*2*param.Nd + 1 + param.Nd;
        vdot_a_idx(idx:idx+param.Nd-1) = ones(param.Nd,1);
    end
    f       = get_f(param, x, t);
    vdot_a  = f(vdot_a_idx==1);
end