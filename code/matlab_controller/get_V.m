function V = get_V(param,x,t)
    eta = get_eta(param,x,t);
    P = get_P(param,x,t);
    V = eta'*P*eta;
end