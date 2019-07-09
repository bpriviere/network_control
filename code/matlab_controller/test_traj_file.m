function test_traj_file(param,i,fn)

    data = csvread(fn,1,0);

    dt = 0.01;
    x = [];
    y = [];
    t = [];
    for i_data = 1:size(data,1)
        duration_i = data(i_data,1);
        t_i = 0: dt: duration_i;

        poly_x = fliplr(data(i_data,2:9));
        poly_y = fliplr(data(i_data,10:17));

        x_i = polyval(poly_x, t_i);
        y_i = polyval(poly_y, t_i);

        x = horzcat(x, x_i);
        y = horzcat(y, y_i);
%         t = horzcat(t, t_i);

    end

    figure()
    plot( x, y, 'linewidth',5);
    axis equal
    
%     figure()
%     subplot(1,2,1)
%     plot( t, x, 'linewidth',5);
%     subplot(1,2,2)
%     plot( t, y, 'linewidth',5);
%     axis equal
    

end