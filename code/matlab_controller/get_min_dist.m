function min_dist = get_min_dist(param,x)

    min_dist = Inf;
    for i = 1:param.N
        pose_i = get_p(param,x,i);
        for j = 1:param.N
            if i~=j
                pose_j = get_p(param,x,j);
                dist = norm(pose_i-pose_j);
                if min_dist > dist
                    min_dist = dist;
                end
            end
        end
    end
end