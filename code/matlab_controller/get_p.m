function p_i = get_p( param, x, agent_idx)
    idx = (agent_idx-1)*2*param.Nd + 1;
    p_i = x(idx: idx + param.Nd - 1);
end