%% Inputs
%point: origin point of ray
%vector:direction of ray
%r1: inner hemi radius
%r2: outer hemi radius
%intensity_in: intensity of beam
%i_cas : where the beam is now
%kappa: absorption coefficient
%s: length of the cubic unit cell side on footprint
%h: height of unit cell
% i_cas: where is the bundle originating from

% i_cas possibilities:

% 1: outside shell
% 2: outer hemi going out
% 3: outer hemi going in
% 4: inner hemi going in
% 5: inner hemi going out
% 6: substrate within shell
% 7: substrate under shell

% intersection_scenario:

% 1: outside hemi moving out
% 2: outer hemi moving in
% 3: inner hemi moving out
% 4: inner hemi moving in
% 5: planar substrate
% 6: out the top
% 7: out the sides

% [bool1,int_point1,path_length1] = intersect_sphere2(r1,point,vector);
% 
% [bool2,int_point2,path_length2] = intersect_sphere2(r2,point,vector);
% 
% [bool3,int_point3,path_length3] = intersect_plane(point,vector);
% 
% [bool4,int_point4,path_length4,normalUC,intersection_scenarioUC] = cubicUnitCell(point,vector,s,h);


function [int_point,n,intersection_scenario,intensity_out, Aglass,path_length] = ray_trace3(point,vector,r1,r2,intensity_in,i_cas,kappa,s,h)

    % Normalization and necessary vector
    vector = vector ./(norm(vector));
    substrate_normal = [0 0 1];
    
    switch i_cas
        
        case 1 %outside shell
         [bool1,int_point1,path_length1] = intersect_sphere2(r2,point,vector);
         [bool2,int_point2,path_length2] = intersect_plane(point,vector);
         [bool3,int_point3,path_length3,normalUC,intersection_scenarioUC] = cubicUnitCell(point,vector,s,h);

         int_cell = [bool1,int_point1,path_length1,2;bool2,int_point2,path_length2,5;bool3,int_point3,path_length3,intersection_scenarioUC];

         int_cell = sortrows(int_cell,[5,1]);

         int_point = int_cell(1,2:4);

         path_length = int_cell(1,5);

         intersection_scenario = int_cell(1,6);

         intensity_out = intensity_in;

         Aglass = 0;

         %Normal Vectors
         switch intersection_scenario
             case 2 %outer hemi moving in
                 n = int_point./(norm(int_point));
             case 5 %substrate
                 n = [0 0 1];
             case 6 %out of the top
                 n = normalUC;
             case 7 %hitting a side wall
                 n = normalUC;
             otherwise
                 error('You have goofed somewhere');
         end
            
        case 2 %outer hemi going out
            [bool1,int_point1,path_length1] = intersect_plane(point,vector);
            [bool2,int_point2,path_length2,normalUC,intersection_scenarioUC] = cubicUnitCell(point,vector,s,h);
            
            int_cell = [bool1,int_point1,path_length1,5;bool2,int_point2,path_length2,intersection_scenarioUC];
            
            int_cell = sortrows(int_cell,[5,1]);

            int_point = int_cell(1,2:4);

            path_length = int_cell(1,5);

            intersection_scenario = int_cell(1,6);

            intensity_out = intensity_in;

            Aglass = 0;
            
            switch intersection_scenario
                
                case 5 %substrate
                    n = [0 0 1];
                    
                case 6 %out the top
                    n = normalUC;
                
                case 7 %out the side
                    n = normalUC;
                otherwise 
                   error('You have goofed somewhere')
            end
     
        case 3 %outer hemi going in
         
         [bool1,int_point1,path_length1] = intersect_sphere2(r1,point,vector);
         [bool2,int_point2,path_length2] = intersect_sphere2(r2,point,vector);
         [bool3,int_point3,path_length3] = intersect_plane(point,vector);
         
         int_cell = [bool1,int_point1,path_length1,4;bool2,int_point2,path_length2,1;bool3,int_point3,path_length3,5];

         int_cell = sortrows(int_cell,[5,1]);

         int_point = int_cell(1,2:4);

         path_length = int_cell(1,5);

         intersection_scenario = int_cell(1,6);
         
         switch intersection_scenario
             case 1 %outer hemi moving out
                 n = -1.*int_point;
                 n = n./norm(n);
             case 4 %inner hemi moving in
                 n = int_point./(norm(int_point));
             case 5 %substrate
                 n = [0 0 1];
             otherwise
                 error('You have goofed somewhere');
         end
         
            % Attenuation in the glass 
            intensity_out = intensity_in * exp(-kappa*path_length); %shell absorbing
            Aglass = intensity_in - intensity_out;

            
        case 4 %inner hemi going in
            [bool1,int_point1,path_length1] = intersect_sphere2(r1,point,vector);
            [bool2,int_point2,path_length2] = intersect_plane(point,vector);
            
            int_cell = [bool1,int_point1,path_length1,3;bool2,int_point2,path_length2,5];
            
            int_cell = sortrows(int_cell,[5,1]);

            int_point = int_cell(1,2:4);

            path_length = int_cell(1,5);

            intersection_scenario = int_cell(1,6);

            intensity_out = intensity_in;

            Aglass = 0;
            
            switch intersection_scenario
                
                case 3 %inner hemi moving out
                    n = -1.*int_point;
                    n = n./norm(n);
                
                case 5 %substrate
                    n = [0 0 1];             
                
                otherwise
                    error('You have goofed somewhere');         
            end
            
        case 5 %inner hemi going out
            [bool1,int_point1,path_length1] = intersect_sphere2(r2,point,vector);
            [bool2,int_point2,path_length2] = intersect_plane(point,vector);  
            
            int_cell = [bool1,int_point1,path_length1,1;bool2,int_point2,path_length2,5];
            
            int_cell = sortrows(int_cell,[5,1]);

            int_point = int_cell(1,2:4);

            path_length = int_cell(1,5);

            intersection_scenario = int_cell(1,6);
            
            switch intersection_scenario
                
                case 1 %outer hemi moving out
                    n = -1.*int_point;
                    n = n./norm(n);
                
                case 5 %substrate
                    n = [0 0 1];
                
                otherwise
                    error('You have goofed somewhere')
                    
            end
            
            % Attenuation in the glass 
            intensity_out = intensity_in * exp(-kappa*path_length); %shell absorbing
            Aglass = intensity_in - intensity_out;
        
        case 6 %substrate under shell
            [bool1,int_point,path_length] = intersect_sphere2(r1,point,vector);
            intersection_scenario = 3;
            intensity_out = intensity_in;
            Aglass = 0;
            n = -1.*int_point;
            n = n./norm(n);                   
        
        case 7 %substrate within shell
            [bool1,int_point1,path_length1] = intersect_sphere2(r1,point,vector);
            [bool2,int_point2,path_length2] = intersect_sphere2(r2,point,vector);
            
            int_cell = [bool1,int_point1,path_length1,4;bool2,int_point2,path_length2,1];
            
            int_cell = sortrows(int_cell,[5,1]);

            int_point = int_cell(1,2:4);

            path_length = int_cell(1,5);

            intersection_scenario = int_cell(1,6);
            
            switch intersection_scenario
                
                case 1 %outer hemi moving out
                    n = -1.*int_point;
                    n = n./norm(n);                    
                
                case 4 %inner hemi moving in
                    n = int_point./(norm(int_point));                    
                
                otherwise
                    error('You have goofed somewhere');
            end

            % Attenuation in the glass 
            intensity_out = intensity_in * exp(-kappa*path_length); %shell absorbing
            Aglass = intensity_in - intensity_out;
            

        
        otherwise
            error('i_cas integer may not be correct. Make sure beam location is correctly identified')      
    end
    
    
    
    
    
    
    
    
    
% Visualize if debugging    
%dbg = true;

if exist('dbg','var')&&dbg
    %plot the cube first
    cube_vertices = [s/2 s/2 0; s/2 -s/2 0;-s/2 s/2 0; -s/2 -s/2 0; s/2 s/2 s; s/2 -s/2 s;-s/2 s/2 s; -s/2 -s/2 s];

    cube_faces = [1 2 4 3; 1 3 7 5; 2 4 8 6; 1 2 6 5; 3 4 8 7; 5 6 8 7];

    cube = patch('Vertices', cube_vertices, 'Faces', cube_faces, 'FaceColor', 'g','facealpha',0.1,'edgealpha',0.1);
    view(3)
    
    %plot the hemispherical shell
    hold on
    C = zeros(101,101);

    R2 = r2;
    [X2,Y2] = meshgrid(-.01:(.01/50):.01);
    Z2 = sqrt(R2.^2 - X2.^2 - Y2.^2);
    Z2(imag(Z2) ~= 0) = 0;
    outerHemi = mesh(X2,Y2,Z2,C,'facealpha',0,'edgealpha',0.1);

    hold on

    R1 = R2*NST;
    [X1,Y1] = meshgrid(-.01:(.01/50):.01);
    Z1 = sqrt(R1.^2 - X1.^2 - Y1.^2);
    Z1(imag(Z1) ~= 0) = 0;
    innerHemi = mesh(X1,Y1,Z1,C,'facealpha',0,'edgealpha',0.1);
    
    axis equal

    hold on

        
end





end