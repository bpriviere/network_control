function dvdotdvb = get_dvdota_dvb(param, x, k)


    dvdotdvb = nan(param.Na*param.Nd, param.Nb*param.Nd);
    A = get_A(param, x);
    for i = 1:param.Na
        
        for j = 1:param.Nb
            
            temp = zeros(param.Nd);
            if 1 % ~(i==j)
                temp = A(i,j)*param.kv*eye(param.Nd);
            end
            
            i_idx = (i-1)*param.Nd + 1;
            j_idx = (j-1)*param.Nd + 1;
            dvdotdvb(i_idx:i_idx + param.Nd - 1,...
                j_idx:j_idx + param.Nd - 1) = temp;
        end
    end

end