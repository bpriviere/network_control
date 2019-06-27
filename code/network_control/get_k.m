function k = get_k( param, t)
    [~,k] = min(abs( param.t - t));
end