function u = get_u( param, x, t)

    u = zeros( param.m,1);
%     u = ones(param.m,1);
    if and( strcmp( param.controller, 'fdbk'), param.Nb ~= 0)
        
        k = floor(t/param.dt);
        if k > param.nt, k = param.nt; end
        
        eta = get_eta(param,x,t);
        
        [v_a, v_b] = get_v_all(param,x);
        vdot_a = get_vdot_all(param,x,t);
        xtilde_dot = [ v_a; vdot_a; v_b];
        
        dvdotdpa = get_dvdota_dpa(param, x, t);
        dvdotdva = get_dvdota_dva(param, x, t);
        dvdotdpb = get_dvdota_dpb(param, x, t);
        dvdotdvb = get_dvdota_dvb(param, x, t);
        dvdotdxtilde = [dvdotdpa, dvdotdva, dvdotdpb];
        
        A = param.my_1' * dvdotdxtilde * xtilde_dot - param.ad(:,k);
        B = param.my_1' * dvdotdvb; 
        K = param.k_fdbk*kron( eye(param.Nd), ones(1,param.gamma));
        
        u = pinv(B)*(-A - K*eta);
        
%     elseif strcmp( param.controller, 'QP')
%         
%         f = get_f(param,x,t);
%         g = get_g(param,x,t);
%         
%         % stabilization
%         V = get_V(param,x,t);
%         dVdeta = get_dVdeta(param,x,t);
%         detadx = get_detadx(param,x,t);
%         detadt = get_detadt(param,x,t);
%         dVdx = dVdeta*detadx;
%         dVdt = dVdeta*detadt; % (partial time derivative)
%         LgV = dVdx*g;
%         LfV = dVdx*f + dVdt;
%         
%         % safety
%         h = get_h(param,x,t);
%         dhdL = get_dhdL(param,x,t);
%         dLdx = get_dLdx(param,x,t);
%         dhdx = dhdL*dLdx;
%         Lfh  = dhdx*f;
%         Lgh  = dhdx*g;
%         
%         % optimization
%         m = param.m;
%         p_v = param.p_v;
%         p_u = param.p_u;
%         lambda_v = param.lambda_v;
%         lambda_h = param.lambda_h;
%         cvx_begin 
%             variables u(m,1)
%             variables delta_v
%             minimize norm(u) + p_v*delta_v + p_u*(u - u_prev)
%             subject to
%                 Lfh + Lgh*u + lambda_h*h >= 0
%                 LfV + LgV*u + lambda_v*V <= delta_v
%                 delta_v >= 0
%         cvx_end
    end

end