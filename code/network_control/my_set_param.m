

function param = my_set_param()

    % see documentation for parameters

    % switch
    param.ode45_on = 0;
    param.gif_on = 1;
    
    % discrete
    param.dt = 0.01;
    param.tf = 10;
    param.t = (param.dt:param.dt:param.tf)';
    param.nt = length(param.t);
    
    % analysis
    param.Nd = 2;
    param.Na = 4;
    param.Nb = 2;
    param.R_comm = 3;
    param.R_des = 1;    
    param.lambda_h = 1;
    param.delta_h = 1;
    param.lambda_v = 1;
    param.delta_v = 1;
    param.p_v = 100;
    param.p_u = 1;
    param.lambda_a = 0.01;
    param.kx = 1;
    param.kv = 1;
    param.k_fdbk = 10; 
    param.gamma = 3;
    param.model = 'reynolds';
    param.controller = 'fdbk';

    % environment
    param.p_lim = 2;
    param.v_lim = 0.01;
    
    % plotting
    param.fontsize = 20;
    param.linewidth = 3;
    param.color_a = [0, 0.4470, 0.7410];
    param.color_b = [0.8500, 0.3250, 0.0980];
    param.color_d = [0.4660, 0.6740, 0.1880];
    param.color_c = [0.4940, 0.1840, 0.5560];
    param.marker_start = '^';
    param.marker_stop = 'o';
    param.marker_free = 'o';
    param.marker_driver = 's';
    param.markersize = 8;
    
    % calculated parameters
    param.N = param.Na + param.Nb;
    param.n = 2*param.Nd*param.N;
    param.m = param.Nd*param.Nb;
    param.H = eye(param.m);
    param.my_1 = kron(ones(param.Na,1),eye(param.Nd))/param.Na;
    
    % desired traj
    param.case_xd = 1; 
    param.R_xd = 3;
    param = init_xd(param);
    
    % filenames
    param.gif_fn = strcat(pwd,'/','animated.gif');
    param.gif_buffer = 0.5;
    
end
