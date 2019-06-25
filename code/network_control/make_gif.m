function make_gif(param,x,t,P)

    for k = 1:param.nt - 1
        if k == 1
            xl = [min( [min(min(P(1,:,:))), min(param.xd(1,:))])-param.gif_buffer, ...
                max([max(max(P(1,:,:))), max(param.xd(1,:))])+param.gif_buffer];
            yl = [min( [min(min(P(2,:,:))), min(param.xd(2,:))])-param.gif_buffer, ...
                max([max(max(P(2,:,:))), max(param.xd(2,:))])+param.gif_buffer];
            
            param.gif_fig = figure();
            axis off;
            set(gca,'ylim',yl);
            set(gca,'xlim',xl);
            gif( param.gif_fn);
        end
        plot_state_instant(param,x(:,k), k)
        gif;
    end
end