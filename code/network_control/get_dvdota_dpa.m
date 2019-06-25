function dvdotdpa = get_dvdota_dpa(param, x, k)

    dvdotdpa = nan(param.Na*param.Nd, param.Na*param.Nd);
    A = get_A(param, x);
    for i = 1:param.Na
        p_i = get_p( param, x, i);
        
        for j = 1:param.Na
            p_j = get_p(param, x, j);
            r_ij = p_j - p_i;
            
            temp = zeros(param.Nd);
            if ~(i==j)
                
                temp = A(i,j)*param.kx*...
                    ( 1 + param.R_des/norm(r_ij)*(...
                    eye(param.Nd) - (r_ij*r_ij')/norm(r_ij)));

%                 temp = temp + A(i,j)*param.kx*(...
%                    eye(param.Nd) - param.R_des*( eye(param.Nd)/norm(r_ij)...
%                    - 3*(r_ij*r_ij')/(norm(r_ij)^7)));
            end
        
            i_idx = (i-1)*param.Nd + 1;
            j_idx = (j-1)*param.Nd + 1;
            dvdotdpa(i_idx:i_idx + param.Nd - 1,...
                j_idx:j_idx + param.Nd - 1) = temp;
        end
    end
end