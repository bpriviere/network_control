function dVdeta = get_dVdeta(param, x, t)
    eta = get_eta(param, x, t);
    P = get_P(param, x, t);
    dVdeta = 2*eta'*P;
end
