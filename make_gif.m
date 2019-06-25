function make_gif(param,x,t,P)

    for k = 1:param.nt - 1
        if k == 1
            xl = [min(min(P(1,:,:)))-param.gif_buffer, ...
                max(max(P(1,:,:)))+param.gif_buffer];
            yl = [min(min(P(2,:,:)))-param.gif_buffer, ...
                max(max(P(2,:,:)))+param.gif_buffer];
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