function eta_new = my_permute( param, eta_old)

    idx = zeros( size(eta_old,1),1);
    count = 0;
    for i_d = 1:param.Nd
        for i_g = 1:param.gamma
            count = count + 1;
            idx(count) = (i_g-1)*param.Nd + i_d;
        end
    end
    
    eta_new = eta_old( idx,:);
%     eta_new = eta_old;
end