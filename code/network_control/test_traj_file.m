function test_traj_file(fn)

    data = csvread(fn,1,0);

    dt = 0.01;
    nt = ceil(sum(data(:,1)/dt));
    x = [];
    y = [];
    for i_data = 1:size(data,1)
        duration_i = data(i_data,1);
        t = 0: dt: duration_i;

        poly_x = fliplr(data(i_data,2:9));
        poly_y = fliplr(data(i_data,10:17));

        x_i = polyval(poly_x, t);
        y_i = polyval(poly_y, t);

        x = horzcat(x, x_i);
        y = horzcat(y, y_i);

    end

    figure()
    plot( x, y, 'linewidth',5);
    axis equal

end