

function test_polyfit(param,P,t)
    
    poly_traj = [];
    for agent_i = 1:param.N
        
        fn = strcat(param.polyfit_base_fn,sprintf('%d',agent_i),'.csv');
        data_agent = csvread( fn,1,0);
        
        x = [];
        y = [];
        for row_i = 1:size(data_agent,1)
            poly_t = param.polyfit_dt:param.polyfit_dt:data_agent(row_i,1);
            
            poly_x = fliplr( data_agent(row_i,2:9));
            poly_y = fliplr( data_agent(row_i,10:17));
            
            x_i = polyval( poly_x, poly_t);
            y_i = polyval( poly_y, poly_t);
            
            x = horzcat( x, x_i);
            y = horzcat( y, y_i); 
        end
        
        poly_traj(agent_i, :,1) = x; 
        poly_traj(agent_i, :,2) = y; 
    end

    t_poly = (1:size(poly_traj,2))*param.polyfit_dt;
    
    for agent_i = 1:param.N
        figure()
        title(sprintf('Agent %d\n', agent_i))
        subplot(1,2,1)
        plot( t_poly, poly_traj( agent_i, :, 1),'b');
        plot( t, P(1, 
    

end