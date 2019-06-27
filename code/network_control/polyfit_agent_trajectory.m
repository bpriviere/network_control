
function data = polyfit_agent_trajectory( param, pose_i, t)

    % FOR 2D!!!
    
    polyfit_t = param.polyfit_dt:param.polyfit_dt:t(end);
%     
%     figure()
%     subplot(1,2,1)
%     plot(pose_i(1,:), pose_i(2,:),'linewidth',param.linewidth);
%     title('Pose')
%     subplot(1,2,2)
%     plot(pose_i(1,:), pose_i(2,:),'linewidth',param.linewidth);
%     title('Polyfit')
    
    data = zeros(length(polyfit_t),33);
    for i_t = 1:length(polyfit_t)
        
        % find corresponding time index
        try
            idx_t_start = find( t > polyfit_t(i_t-1),1);
        catch
            idx_t_start = 1;
        end
        
        idx_t_end = find( t > polyfit_t(i_t), 1);
        if isempty(idx_t_end)
            idx_t_end = length(t);
        end
        idxs = idx_t_start:idx_t_end;
        
        t_clip = t(idxs) - t(idx_t_start);
        pose_x_clip = pose_i(1,idxs)';
        pose_y_clip = pose_i(2,idxs)';
        
        poly_x =  polyfit(t_clip, pose_x_clip,7);
        poly_y =  polyfit(t_clip, pose_y_clip,7);
        duration = t(idx_t_end) - t(idx_t_start);
        
        data(i_t, 1:17) = [duration, fliplr(poly_x), fliplr(poly_y)];
        
        t_test = t_clip(1):0.01:t_clip(end);
        
%         subplot(1,2,1);
%         plot( pose_x_clip, pose_y_clip, 'x--','linewidth',param.linewidth);
%         
%         subplot(1,2,2);
%         plot( polyval(poly_x,t_test), polyval(poly_y, t_test), ...
%             'x--', 'linewidth',param.linewidth);
    end
end