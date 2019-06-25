function plot_state_instant(param, x, k)

cla(gca);

if param.Nd == 2
    
    A = get_A(param,x);
    for i = 1:param.N
        pose_i = get_p(param,x,i);
        for j = i+1:param.N
            pose_j = get_p(param,x,j);
            if A(i,j) ~= 0
                plot( [pose_i(1) pose_j(1)], [pose_i(2) pose_j(2)],...
                    'linewidth',param.linewidth/2,...
                    'color', param.color_a)
            end
        end
    end
    
    for i = 1:param.Na
        pose = get_p(param,x,i);
        
        plot( pose(1), pose(2), ...
            'color', param.color_a, ...
            'Marker',param.marker_free, ...
            'MarkerSize', 2*param.markersize,...
            'MarkerFaceColor',[1 1 1])
    end
    
    for i = param.Na+1:param.N
        pose = get_p(param,x,i);
        
        plot( pose(1), pose(2), ...
            'color', param.color_b, ...
            'Marker',param.marker_driver, ...
            'MarkerSize', 2*param.markersize, ...
            'MarkerFaceColor',[1 1 1])
    end
    
    [pose_a, ~] = get_p_all( param,x);
    xbar = param.my_1'*pose_a;
    
    plot( xbar(1), xbar(2), 'color', param.color_c,...
        'Marker','.','Markersize', 2*param.markersize)
    
    plot( param.xd(1,k), param.xd(2,k), 'color', param.color_d,...
        'Marker','.','Markersize', 2*param.markersize)
end

end