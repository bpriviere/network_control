function detadx = get_detadx(param, x, k)

    dvdotdpa = get_dvdota_dpa(param, x, k);
    dvdotdva = get_dvdota_dva(param, x, k);
    dvdotdpb = get_dvdota_dpb(param, x, k);
    dvdotdvb = get_dvdota_dvb(param, x, k);
    
    my_1 = param.my_1;
    O_a = zeros( param.Nd, param.Nd*param.Na);
    O_b = zeros( param.Nd, param.Nd*param.Nb);
    detadx = [...
        my_1', O_a, O_b, O_b;
        O_a, my_1', O_b, O_b; 
        my_1'*dvdotdpa, my_1'*dvdotdva, ...
        my_1'*dvdotdpb, my_1'*dvdotdvb];
    
    
    detadx = my_permute(param, detadx);
    
        
end
