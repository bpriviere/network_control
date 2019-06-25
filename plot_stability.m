function plot_stability(param,V)

    Vdot = gradient(V,param.t);

    figure()
    h1 = plot( param.t, V, 'linewidth',param.linewidth);
    h2 = plot( param.t, Vdot, 'linewidth',param.linewidth);
    hleg = legend([h1,h2],'$$V$$','$$\dot{V}$$','location','best');
    set(hleg,'Interpreter','latex');
    set(gca,'fontsize',param.fontsize);
    grid on

end