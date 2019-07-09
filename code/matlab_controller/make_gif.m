function make_gif(param,x,t,P)

    for i_gif_nt = 1:param.gif_nt 
        if i_gif_nt == 1
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
        
        i_k = round( i_gif_nt/param.gif_nt*param.nt);
        [~,i_t] = min( abs( i_k*param.dt - t) );
        plot_state_instant(param,x(:,i_t), t(i_t))
        gif;
    end
end