function param = init_xd(param)

    t = param.t;
    switch param.case_xd
        
        case 0 % straight line in 2d
            param.xd = ones( length(t),2)';
            
        case 1 % circle in 2d
            param.xd = param.R_xd*[ ...
                cos(t./param.tf*2*pi), sin(t./param.tf*2*pi)]';
        
        case 2 % circle in 3d
            param.xd = param.R_xd*[ ...
                cos(t./param.tf*2*pi), sin(t./param.tf*2*pi), t]';
        
        case 3 % sin wave
            param.xd = [ ...
                param.R_xd*10*sin(t./param.tf*2*pi), t]';
            
    end
    
    if param.Nd == 2
        param.vd = gradient(param.xd, t, t);
        param.ad = gradient(param.vd, t, t);
        param.jd = gradient(param.ad, t, t);
    elseif param.Nd == 3
        for i = 1:param.Nd
            param.vd(i,:) = gradient(param.xd(i,:), t);
            param.ad(i,:) = gradient(param.vd(i,:), t);
            param.jd(i,:) = gradient(param.ad(i,:), t);
        end
    end
end