function [bool,int_point,path_length] = intersect_plane(point,dir)
    

    dir = dir./(norm(dir)); 
    
    % Ray is either pointing away from the substrate or parallel to it, no
    % intersection of concern
    if dir(3) >= 0

        bool = 0;
        
        int_point = [NaN NaN NaN];
        
        path_length = inf;
        
        %error('Ray does not intersect plane at all, or points away from the plane (backwards intersection)')
    
    % Ray intersects the substrate
    else
        
        bool = 1;
        
        D = -point(3) / (dir(3));

        x2 = D*(dir(1)) + point(1);

        y2 = D*(dir(2)) + point(2);

        z2 = 0;

        int_point = [x2 y2 z2];

        path_length = norm(point - int_point);
    
    end

end