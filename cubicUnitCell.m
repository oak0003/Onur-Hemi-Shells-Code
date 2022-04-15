%o:origin of ray
%n:direction vector of ray
%s:side length of cubic unit cell

% x and y ranges from -s/2 to s/2, z ranges from 0 to h.
function [bool,int_point,path_length,normal,intersection_scenario] = cubicUnitCell(o,n,s,h)
% Draw the cubic unit cell (no hemisphere) for visualization for debug
%dbg = true;

if exist('dbg','var')&&dbg
    cube_vertices = [s/2 s/2 0; s/2 -s/2 0;-s/2 s/2 0; -s/2 -s/2 0; s/2 s/2 h; s/2 -s/2 h;-s/2 s/2 h; -s/2 -s/2 h];

    cube_faces = [1 2 4 3; 1 3 7 5; 2 4 8 6; 1 2 6 5; 3 4 8 7; 5 6 8 7];

    patch('Vertices', cube_vertices, 'Faces', cube_faces, 'FaceColor', 'g','facealpha',0.1,'edgealpha',0.1);
    view(3)
    hold on
    % Draw Circle at the base
    % number of points
    numPoints = 1000;

    %// running variable
    t = linspace(0,2*pi,numPoints);
    
    %Radius of circumscribing circle
    r = s* (sqrt(2)/2);
    %Center of circle
    c = [0 0];

    x = c(1) + r*sin(t);
    y = c(2) + r*cos(t);

    %// draw line
    line(x,y);

    %// or draw polygon if you want to fill it with color
    %// fill(x,y,[1,1,1])
%     axis equal
end % if we're debugging, visualize circle intersect

%First check to see if beam leaves unit cell from the top plane
plane_normal = [0 0 -1];
plane_point = [0 0 h];
[I,rc] = line_plane_intersection(n, o, plane_normal, plane_point);

ix = I(1); iy = I(2); iz = I(3);

%Check to see if the beam is leaving through the top

%To see if it's at the edges
cornerXCheck = -.0000000001 < (abs(ix) -(s/2)) && (abs(ix) -(s/2)) < .0000000001;

cornerYCheck = -.0000000001 < (abs(iy) -(s/2)) && (abs(iy) -(s/2)) < .0000000001;

heightCheck = -.0000000001 < (iz -h) && (iz -h) < .0000000001;

IequalO = ix == o(1) && iy == o(2) && iz == o(3);

if (rc == 1 && abs(ix) < (s/2) && abs(iy) < (s/2) || (rc == 1 && (cornerXCheck || cornerYCheck))) && heightCheck && ~IequalO
    int_point = I;
    normal = [0 0 -1];
    path_length = norm(o - int_point);
    bool = 1;
    intersection_scenario = 6;
    
    if exist('dbg','var')&&dbg
        plot3([o(1) ix],[o(2),iy],[o(3),iz]);
    end % if we're debugging, visualize circle intersect
    
%If not leaving through the top, which side plane does the beam leave
%through?
else
    %Radius of circumscribing circle
    r = s* (sqrt(2)/2);
    ox = o(1); oy = o(2); oz = o(3);

    nx = n(1); ny = n(2); nz = n(3);
    
    lx   = -ox;
    ly   = -oy;
    l2   = lx^2+ly^2; % distance from ray origin to circle center
    nxy2 = sqrt(nx^2+ny^2);
    nxx  = nx/nxy2;   % scale it for a unit normal in 2d
    nyy  = ny/nxy2;
    loc  = lx*nxx+ly*nyy; % distance from origin to closest approach
    lca2 = l2-loc^2;      % closest approach squared
    lct2 = r^2-lca2;      % distance from closest approach to circle intersection
    s4xy = loc+sqrt(lct2);% scale from origin to interception

    % which means the location of intersection with the circle is
    tx   = ox+s4xy*nxx;
    ty   = oy+s4xy*nyy;
    tz   = oz+s4xy/nxy2*nz;
    
     if exist('dbg','var')&&dbg
            plot3([ox tx],[oy,ty],[oz,tz]);
     end % if we're debugging, visualize circle intersect
    
    % Check to see if the beam is incident on a vertical corner. If it is incident
    % on a corner than the circle intersect is also the square intersect
    % for the circle and square are coincident at corners
    if -.0000000001 < (abs(tx) - abs(ty)) &&  (abs(tx) - abs(ty)) < .0000000001
        
        %Top Right Corner
        if tx > 0 && ty > 0
            wx = 1;wy = 1;
        
        %Top Left Corner
        elseif tx < 0 && ty > 0
            wx = -1;wy = 1;
        
        %Bottom Right
        elseif tx > 0 && ty < 0
            wx = 1;wy = -1;
            
        
        %Bottom Left    
        elseif tx < 0 && ty < 0
            wx = -1;wy = -1;
        
        %Error    
        else
            error('Something went wrong with corners');
        end
        
    % If it's not a corner
    else
        
        wx=sign(tx)*(abs(tx)>abs(ty));
        wy=sign(ty)*(abs(ty)>abs(tx));
        
        s4   = ((s/2)-ox*wx-oy*wy)/(nx*wx+ny*wy);
            
        tx   = ox+s4*nx;
        ty   = oy+s4*ny;
        tz   = oz+s4*nz;
    end
    
    
    
    if exist('dbg','var')&&dbg
        plot3([ox tx],[oy,ty],[oz,tz]);
    end % if we're debugging, visualize circle intersect
    
%     aa = (dx^2+dy^2);
%     bb = (2*ox*dx+2*oy*dy);
%     cc = (ox^2+oy^2-r^2);
%     ss = (-bb+sqrt(bb^2-4*aa*cc))/(2*aa);
%    
    
%     tx = ox+dx*ss; % x intercept of circle
%     ty = oy+dy*ss; % y intercept of circle
    
    %THIS IS MORE GENERAL, GOT CONFUSING, MAY USE FOR HEX UNIT CELL
%     % so the angle from the x axis to our circle intersection is
%     qi   = (sign(ty)*((tx<0)*pi+sign(tx)*atan(abs(ty/tx))))
%     
%     % now to find the angle of intersection with the wall of interest
%     % rather than the circle
%     qw = (floor(qi + (pi/4))/(pi/2))*(pi/2)
    
%     wx=sign(tx)*(abs(tx)>abs(ty));
%     wy=sign(ty)*(abs(ty)>abs(tx));

%     s4   = ((s/2)-ox*wx-oy*wy)/(nx*wx+ny*wy);
%             
%     tx   = ox+s4*nx;
%     ty   = oy+s4*ny;
%     tz   = oz+s4*nz;
    
%     if exist('dbg','var')&&dbg
%         plot3([ox tx],[oy,ty],[oz,tz])
%     end % if we're debugging, visualize wall intersect


    int_point = [tx,ty,tz];

    normal = [-wx -wy 0];
    
    path_length = norm(o - int_point);
    
    bool = 1;

    intersection_scenario = 7;
    
    % If the intersection with the side walls has a negative z value, then
    % substrate is hit first. Output no intersection data
    if int_point(3) < 0
        
        int_point = [NaN NaN NaN];
        
        bool = 0;
        
        intersection_scenario = NaN;
        
        normal = NaN;
        
        path_length = inf;
    end
    
end


end