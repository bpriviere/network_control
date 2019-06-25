
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

% plots
plot_state_space(param, P, t);

if param.gif_on 
    fprintf('GIF...')
    make_gif(param,x,t,P);
    fprintf('Complete!\n')
end
