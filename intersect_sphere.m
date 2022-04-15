function [bool,int_point,path_length] = intersect_sphere(r,point,dir)

    
    % if originating point is on the destination sphere
%     if norm(point) == r
%         bool = 0;
%         int_point = NaN;
%         path_length = inf;
%     end
    
    dir = dir./(norm(dir));   
    
    a = dir(1)^2 + dir(2)^2 + dir(3)^2;
    b = 2*(point(1) * dir(1) + point(2) * dir(2) + point(3) * dir(3));
    c = point(1)^2 + point(2)^2 + point(3)^2 - r^2;
    
    discriminant = (b^2) - (4*a*c);
    
    % Ray heading towards a larger radius from origin or on same sphere
       if r >= norm(point - [0 0 0])
    
        if discriminant < 0 

            %error('Ray does not intersect sphere at all')
            
            bool = 0;
            
            int_point = NaN;
            
            path_length = inf;

        elseif discriminant == 0 

            t = -b/(2*a);

            int_point = [point(1) + t*dir(1) , point(2) + t*dir(2) , point(3) + t*dir(3)];

            path_length = norm(point - int_point);
            
            bool = 1;



        elseif discriminant > 0

            t1 = (-b + sqrt(discriminant) ) / (2*a);

            t2 = (-b - sqrt(discriminant) ) / (2*a);
            
            bool = 1;

            if t1>t2

                int_point = [point(1) + t1*dir(1) , point(2) + t1*dir(2) , point(3) + t1*dir(3)];

                path_length = norm(point - int_point);

            elseif t2>t1

                int_point = [point(1) + t2*dir(1) , point(2) + t2*dir(2) , point(3) + t2*dir(3)];

                path_length = norm(point - int_point);

            end

        end
       end
       
       % Ray heading towards a smaller radius from origin
       if r < norm(point - [0 0 0])

        if discriminant < 0 

            %error('Ray does not intersect sphere at all')
            
            bool = 0;
            
            int_point = NaN;
            
            path_length = inf;

        elseif discriminant == 0 

            t = -b/(2*a);

            int_point = [point(1) + t*dir(1) , point(2) + t*dir(2) , point(3) + t*dir(3)];

            path_length = norm(point - int_point);
            
            bool = 1;



        elseif discriminant > 0

            t1 = (-b + sqrt(discriminant) ) / (2*a);

            t2 = (-b - sqrt(discriminant) ) / (2*a);
            
            bool = 1;


            if t1<t2

                int_point = [point(1) + t1*dir(1) , point(2) + t1*dir(2) , point(3) + t1*dir(3)];

                path_length = norm(point - int_point);

            elseif t2<t1

                int_point = [point(1) + t2*dir(1) , point(2) + t2*dir(2) , point(3) + t2*dir(3)];

                path_length = norm(point - int_point);

            end
        end
    end
end