function A = get_A(param, x)

A = nan(param.N);

for i = 1:param.N
    p_i = get_p( param, x, i);
    
    for j = 1:param.N 
        p_j = get_p( param, x, j);
        
        if norm( p_j - p_i) < param.R_comm
            A(i,j) = exp( - param.lambda_a * norm( p_j - p_i));
        else
            A(i,j) = 0;
        end
    end
end

% A = ones(param.N);

end