

function param = my_set_param()

    % see documentation for parameters

    % switch
    param.ode45_on = 0;
    param.gif_on = 1;
    param.write_traj_file_on = 0;
    
    % discrete
    param.dt = 0.05;
    param.tf = 100;
    param.t = (param.dt:param.dt:param.tf)';
    param.nt = length(param.t);
    param.polyfit_dt = 10; 
    
    % analysis
    param.Nd = 2;
    param.Na = 6;
    param.Nb = 1;
    param.R_comm = 1.;
    param.R_des  = 0.8; 
    param.lambda_h = 1;
    param.delta_h = 1;
    param.lambda_v = 1;
    param.delta_v = 1;
    param.p_v = 100;
    param.p_u = 1;
    param.lambda_a = 1;
    param.kx = 2;
    param.kv = 2;
    param.k_fdbk = 20; 
    param.gamma = 3;
    param.model = 'reynolds';
    param.controller = 'fdbk';
    param.min_dist_constraint = 0.3;

    % environment
    param.p_lim = 1;
    param.v_lim = 0.2;
    
    % plotting
    param.fontsize = 20;
    param.linewidth = 3;
    param.color_a = [0, 0.4470, 0.7410];
    param.color_b = [0.8500, 0.3250, 0.0980];
    param.color_d = [0.4660, 0.6740, 0.1880];
    param.color_c = [0.4940, 0.1840, 0.5560];
    param.marker_start = '^';
    param.marker_stop = 's';
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
    param.R_xd = .75;
    param = init_xd(param);
    
    % gif
    param.gif_fn = strcat(pwd,'/','animated.gif');
    param.gif_buffer = 2;
    param.gif_nt = 300;
    
    % trajectory file output
    param.polyfit_base_fn = 'polyfit_traj_';
    
    % make header for Wolfgang :)
    param.polyfit_header = {'duration', 'x^0', 'x^1', 'x^2', 'x^3', 'x^4', 'x^5', 'x^6', ...
        'x^7', 'y^0', 'y^1', 'y^2', 'y^3', 'y^4', 'y^5', 'y^6', 'y^7', 'z^0', 'z^1', 'z^2', 'z^3', ...
        'z^4', 'z^5', 'z^6', 'z^7', 'yaw^0', 'yaw^1', 'yaw^2', 'yaw^3', 'yaw^4', 'yaw^5', 'yaw^6', ...
        'yaw^7'}; 
    
    % fig output name
    param.ss_fig_fn = 'state_space';
    
    
end
