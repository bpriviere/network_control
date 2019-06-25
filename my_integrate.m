function [x,V] = my_integrate( param, x0)

    x = nan( param.n, param.nt);
    V = nan( param.nt, 1);
    x(:,1) = x0;
    V(1) = get_V(param,x(:,1),1);
    for k = 1:param.nt-1
        dxdt = get_dxdt( param, x(:,k), param.t(k));
        x(:,k+1) = x(:,k) + dxdt * param.dt;    
        V(k+1) = get_V(param,x(:,k+1),k+1);
    end

end