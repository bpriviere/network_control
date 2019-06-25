function eta = get_eta(param, x, t)

    k = floor(t/param.dt);
    if k > param.nt, k = param.nt; end
    
    [pose_a, ~] = get_p_all(param,x);
    [vel_a, ~] = get_v_all(param,x);
    vdot_a = get_vdot_all(param,x,t);
    
    y   = param.my_1'*pose_a - param.xd(:,k);
    yp  = param.my_1'*vel_a - param.vd(:,k);
    ypp = param.my_1'*vdot_a - param.ad(:,k);
    
    eta = [y; yp; ypp];
    eta = my_permute(param, eta);
    
end