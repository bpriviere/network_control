function dvdotdva = get_dvdota_dva(param, x, k)

    dvdotdva = nan(param.Na*param.Nd, param.Na*param.Nd);
    A = get_A(param, x);
    for i = 1:param.Na
        
        for j = 1:param.Na
            
            temp = zeros(param.Nd);
            if ~(i==j)
                temp = A(i,j)*param.kv*eye(param.Nd);
            end
            
            i_idx = (i-1)*param.Nd + 1;
            j_idx = (j-1)*param.Nd + 1;
            dvdotdva(i_idx:i_idx + param.Nd - 1,...
                j_idx:j_idx + param.Nd - 1) = temp;
        end
    end
end