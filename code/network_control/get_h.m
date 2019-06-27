
function h = get_h( param,x,t)

    A = get_A(param,x);
    D = zeros(size(A));
    for i = 1:param.N
        D(i,i) = sum(A(i,:));
    end
    L = D-A;
    [V,D] = eigs(L);
    
end