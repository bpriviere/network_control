function P = get_P(param, x, t)

    a1 = eye(param.gamma-1);
    a1 = [zeros(param.gamma-1,1) a1];
    a1 = [a1; -1*param.k_fdbk*ones(1,param.gamma)];

    A_cl = a1;
    for i = 1:param.Nd-1
        A_cl = blkdiag(A_cl,a1);
    end

    P = sylvester( A_cl', A_cl, -eye(size(A_cl)));
    % c = 1 / eigs(param.P,1);
end