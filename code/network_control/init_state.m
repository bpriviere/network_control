function x0 = init_state(param)

    x0 = nan(param.n,1);
    
    min_dist = 0;
    while min_dist < param.min_dist_constraint

        for i_na = 1:param.Na    
            p_i = param.p_lim*rand(param.Nd,1) - param.p_lim/2;
            v_i = param.v_lim*rand(param.Nd,1) - param.v_lim/2;

            state_idx = (i_na-1)*2*param.Nd + 1;
            x0( state_idx:state_idx + 2*param.Nd -1 ) = [p_i; v_i];
        end

        for i_nb = 1:param.Nb
            p_i = param.p_lim*rand(param.Nd,1) - param.p_lim/2;
            v_i = param.v_lim*rand(param.Nd,1) - param.v_lim/2;

            state_idx = (i_nb-1)*2*param.Nd + param.Na*2*param.Nd + 1;
            x0( state_idx:state_idx + 2*param.Nd -1 ) = [p_i; v_i];
        end
        
        min_dist = get_min_dist(param,x0);
    end
    
    
%     p_i = [-0.5;-0.5];
%     x0( state_idx:state_idx + 2*param.Nd -1 ) = [p_i; v_i];
%     x0 = [ 0,0,1,1,1,1,1,1,4,4,1,1]';
%     x0 = [ 0,0,0,0,1,1,0,0,4,4,0,0]';
    
end