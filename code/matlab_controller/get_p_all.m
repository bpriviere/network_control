function [pose_a, pose_b] = get_p_all(param,x)

    pose_a = nan(param.Nd*param.Na,1);
    pose_b = nan(param.Nd*param.Nb,1);

    for i = 1:param.Na
        idx = (i-1)*param.Nd + 1;
        pose_a(idx:idx+param.Nd-1) = get_p( param, x, i);        
    end

    for i = 1:param.Nb
        idx = (i-1)*param.Nd + 1;
        pose_b(idx:idx+param.Nd-1) = get_p( param, x, i+param.Na);  
    end
end