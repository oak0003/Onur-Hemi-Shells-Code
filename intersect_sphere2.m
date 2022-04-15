function [bool,int_point,path_length] = intersect_sphere2(r,point,dir)

    a = dot(dir,dir);
    b = 2* dot(dir,point);
    c = dot(point,point) - r^2;
    
    discriminant = b^2 - 4*a*c;
    
    % No intersection
    if discriminant < 0
        bool = 0;
        int_point = [NaN NaN NaN];
        path_length = inf;
    
    % Only a singular intersection (tangent)
    elseif discriminant == 0
        
        t = -b / (2*a);
        
        % Intersection occurs in front
        if t>0
            bool = 1;
            int_point = point + t.*(dir);
            path_length = norm(point - int_point);
            
            % Intersection occurs below the substrate
            if int_point(3) < 0 && bool == 1
                bool = 0;
                int_point = [NaN NaN NaN];
                path_length = inf;
            end
            
            % Intersection either occurs at same point of origin or behind
        else
            bool = 0;
            int_point = [NaN NaN NaN];
            path_length = inf;
        
        end
     
      % Two points possible for intersection
     elseif discriminant > 0
         
         t1 = (-b + sqrt(discriminant))/(2*a);
         
         t2 = (-b - sqrt(discriminant))/(2*a);
         
         % Set either parameter to 0 if it's close enough
         if abs(t1) < .0000000000001
             t1 = 0;
         elseif abs(t2) < .0000000000001
             t2 = 0;
         end
         
         % If both points occur in front of ray, choose the closer one
         if t1 > 0 && t2 > 0
             
             t = min(t1,t2);
             bool = 1;
             int_point = point + t.*(dir);
             path_length = norm(point - int_point);
             
             % If closer one is below the substrate, toss it and return no
             % intersection
             if int_point(3) < 0 && bool == 1
                bool = 0;
                int_point = [NaN NaN NaN];
                path_length = inf;
             end
         
         % Both points behind ray, no intersection
         elseif t1 < 0 && t2 < 0
             
             bool = 0;
             int_point = [NaN NaN NaN];
             path_length = inf;
         
         else
             
             %Other cases, such as one parameter being pos and the other
             %neg, choose the larger one
             t = max(t1,t2);
             
             %If one is zero and the other negative, no intersection
             if (abs(t1) <.0000000001  || abs(t2) <.0000000001) && ((t1+t2) < 0)
                 bool = 0;
                 int_point = [NaN NaN NaN];
                 path_length = inf;
             else
                 bool = 1;
                 int_point = point + t.*(dir);
                 path_length = norm(point - int_point);
                 
                 %toss if below substrate
                 if int_point(3) < 0 && bool == 1
                    bool = 0;
                    int_point = [NaN NaN NaN];
                    path_length = inf;
                 end
             end
         end
                 
    end
     
    

end