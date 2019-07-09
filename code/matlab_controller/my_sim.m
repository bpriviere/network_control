
clear all; close all; clc;

rng(3);
param = my_set_param();

% state = [p^1_a; v^1_a;, ..., p^1_b, v^1_b, ...];
x0 = init_state( param);
t = param.t;

% integrate
fprintf('Simulation...')
if param.ode45_on
    [t,x] = ode45( @(t,x) get_dxdt( param, x, t), [param.t(1), param.t(end)], x0);
    x = x';
else
    [x,V] = my_integrate( param, x0);
    plot_stability(param,V);
end
fprintf('Complete!\n')

% extract position
P = zeros( param.Nd, param.N, length(t));
for k = 1:length(t)
    [p_a, p_b] = get_p_all( param, x(:,k));
    P( :, :, k) = reshape( [p_a; p_b], [param.Nd param.N]);
end

min_dist = Inf;
for k = 1:length(t)
    dist = get_min_dist(param,x(:,k));
    if min_dist > dist
        min_dist = dist;
    end
end
fprintf('Min Dist: %d\n', min_dist);

% plots
plot_state_space(param, P, t);

if param.gif_on 
    try 
        delete(param.gif_fn)
    catch
        '';
    end
    
    fprintf('GIF...')
    make_gif(param,x,t,P);
    fprintf('Complete!\n')
end

if param.write_traj_file_on
    fprintf('Writing Trajectory Files...');
    for agent_i = 1:param.N
        fn = strcat(param.polyfit_base_fn,sprintf('%d',agent_i),'.csv');
        data = polyfit_agent_trajectory(param,squeeze(P(:,agent_i,:)),t);
        
        % 
        try
            delete(fn);
        catch
            '';
        end
        
        % write
        fid = fopen(fn,'w');
        for header_i = 1:length(param.polyfit_header)-1
            fprintf(fid,strcat( string(param.polyfit_header(header_i)),','));
        end
        fprintf(fid,'%s\n',string(param.polyfit_header(end)));
        fclose(fid);
        dlmwrite(fn,data,'-append','precision','%.16f');
        
%         test_traj_file( param,agent_i,fn);
    end
    fprintf('Complete!\n');
end
