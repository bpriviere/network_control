function detadt = get_detadt(param,x,t)
    k = floor( t/param.dt);
    if k > param.nt, k = param.nt; end
    detadt = -[...
        param.vd(:,k); 
        param.ad(:,k); 
        param.jd(:,k)];
    detadt = my_permute(param, detadt);
end
