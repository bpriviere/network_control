function v_i = get_v( param, x, agent_idx)
    idx = (agent_idx-1)*2*param.Nd + 1;
    v_i = x(idx + param.Nd: idx + 2*param.Nd - 1);
end