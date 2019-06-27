function dvdotdpb = get_dvdota_dpb(param, x, k)

    dvdotdpb = nan(param.Na*param.Nd, param.Nb*param.Nd);
    A = get_A(param, x);
    for i = 1:param.Na
        p_i = get_p( param, x, i);
        
        for j = 1:param.Nb
            p_j = get_p(param, x, j + param.Na);
            r_ij = p_j - p_i;
                
            temp = A(i,j)*param.kx*...
                ( 1 - param.R_des/norm(r_ij)*(...
                eye(param.Nd) - (r_ij*r_ij')/norm(r_ij)^2));

%             temp = temp + A(i,j)*param.kx*(...
%                eye(param.Nd) + param.R_des*( eye(param.Nd)/norm(r_ij)...
%                - 3*(r_ij*r_ij')/(norm(r_ij)^7)));
            
        
            i_idx = (i-1)*param.Nd + 1;
            j_idx = (j-1)*param.Nd + 1;
            dvdotdpb(i_idx:i_idx + param.Nd - 1,...
                j_idx:j_idx + param.Nd - 1) = temp;
        end
    end

end