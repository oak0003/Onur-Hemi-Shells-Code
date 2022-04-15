% Non Black Substrate, Variable Angle of Incidence,
% Hemispherical Shell with Cubic Unit Cell

% SINGLE WAVELENGTH

% clc;clear all;clf

function [AglassTotFinal, AsubTotFinal,RTotFinal,check2,Scattered,N,corner1,tangent,lostInten,segments,discardedInt,count] = CubicCellHemi(r2,NST,s,wavelength,m2,m4,incident_polar,incident_azimuth,res)

%% Setting the scene. Geometric parameters and optical parameters

% Geometric Parameters
% r2 = .005; % meters
% NST = .75; % normalized shell thickness
r1 = r2*NST;
% s = 2*r2+.001; % side length, meters
h = r2+.00001; % height, meters

% Optical Parameters
n1 = 1; % unity
% m2 = [1.532, 1.925E-7]; % complex index of refraction for absorbing shell material m = n + ik
n3 = 1; % unity
% m4 = [4.084907216, 0.036329897]; % complex index of refraction for absorbing substrate material m = n + ik
thetaCrit = asind(n3/m2(1)); % critical angle for checking for TIR
% wavelength = .55E-6; % meters, incident wavelength
kappa = (4*pi*m2(2)) / wavelength; %absorption coefficient

%% If we're debugging, visualize the hemispherical shell and cubic unit cell

dbg = true;

if exist('dbg','var')&&dbg
    %plot the cube first
    cube_vertices = [s/2 s/2 0; s/2 -s/2 0;-s/2 s/2 0; -s/2 -s/2 0; s/2 s/2 h; s/2 -s/2 h;-s/2 s/2 h; -s/2 -s/2 h];

    cube_faces = [1 2 4 3; 1 3 7 5; 2 4 8 6; 1 2 6 5; 3 4 8 7; 5 6 8 7];
    
    figure

    cube = patch('Vertices', cube_vertices, 'Faces', cube_faces, 'FaceColor', 'g','facealpha',0.1,'edgealpha',0.1);
    view(3)
    
    %plot the hemispherical shell
    hold on
    C = zeros(101,101);

    R2 = r2;
    [X2,Y2] = meshgrid(-(s/2):((s/2)/50):(s/2));
    Z2 = sqrt(R2.^2 - X2.^2 - Y2.^2);
    Z2(imag(Z2) ~= 0) = 0;
    outerHemi = mesh(X2,Y2,Z2,C,'facealpha',0,'edgealpha',0.1);

    hold on

    R1 = R2*NST;
    [X1,Y1] = meshgrid(-(s/2):((s/2)/50):(s/2));
    Z1 = sqrt(R1.^2 - X1.^2 - Y1.^2);
    Z1(imag(Z1) ~= 0) = 0;
    innerHemi = mesh(X1,Y1,Z1,C,'facealpha',0,'edgealpha',0.1);
    
    axis equal

    hold on

        
end

%% SOURCE PLANE POINT GENERATION

% Incoming linearly collimated radiation's direction
% incident_polar = 30; %incident angle of incoming radiation wrt substrate normal (pos z axis)
% incident_azimuth = 20; %incident angle of incoming radiation wrt pos x axis

% Gridding source plane to get emmission points
% res = 30; %determines resolution of the grid. Grid will be res*2 by res*2 (elements wise)
dims = s;%(r2*1.5);
ss = (s/2) - .0001;
x = -(ss):dims/res:(ss);
y = -1.*(-(ss):dims/res:(ss));

[X,Y] = meshgrid(x,y);

[dim1,dim2] = size(X);

points = [];%zeros(dim1^2,2);


for i = 1:dim1%length(X)
    for j = 1:dim2%length(Y)
       points = [points;X(i,j),Y(i,j)];
    end
end





% Number of emission points on source plane
N = length(points);

% All points are on source plane and have a z value of s
z_vals = (h*ones(N,1));
points = [points,z_vals];

% Direction vector of beams from source plane
dir = [sind(incident_polar+180)*cosd(incident_azimuth), sind(incident_polar+180)*sind(incident_azimuth),cosd(incident_polar+180)];

%% If debugging the source points, visualize
%dbgSource = true;

if exist('dbgSource','var')&&dbgSource

    scatter3(points(:,1), points(:,2),points(:,3));

end

%% Initialize counter variables

tangent = 0; % tangent to either hemisphere
corner1 = 0; % incident on seam of either hemisphere and substrate
corner2 = 0; % incident on corner substrate and side walls
error_tracker = 0; % tracks number of iterations in loops
lostInten = 0;
scatInt = 0;
count = 0;
count2 = 0;
count3 = 0;

thresh = .00075; % attenuation threshold

%% Loop establish the for loop that iterates through all of the beams from the source plane

% our rays will look like this:
% ray(ii,:)={r_org(1x3 double),d_ray(1x3 double),e_act(double),i_cas(uint8)}
% e_act: energy still in this bundle
% e_gls: energy absorbed by glass
% e_sub: energy absorbed by substrate
% e_out: energy that left the cell
% r_org: position of bundle (x,y,z)
% d    : direction of bundle (x,y,z)
% i_cas: where is the bundle originating from

% i_cas possibilities:

% 1: outside shell
% 2: outer hemi going out
% 3: outer hemi going in
% 4: inner hemi going in
% 5: inner hemi going out
% 6: substrate under shell
% 7: substrate within shell

% i_scenario will match i_cas. i_scenario is what surface we hit and will
% determine what happens to our ray at the intersection (splitting,
% transmitting, etc)
AglassTot1 = zeros(1,N);
AsubTot1 = zeros(1,N);
RTot1 = zeros(1,N);
Scattered = [];
segments = 0;
discardedInt = 0;
ii = 1;

% for ii = 1:N
while ii <= N    
    
    
        %% Begin algorithm
    % This sets up the ray cell, initially containing only the first emitted ray
    point = points(ii,:);
    rays = [point,dir,1,1];

    % Define numBeams variable to determine if there are any rays left
    [numBeams,~] = size(rays);

    % Initialize the new ray cell
   %new_rays = [];
    j = 1; % counts number of times the external while loop runs

    % Initialize the property counters
    Aglass = [];
    Asub = [];
    R = [];

    while numBeams ~= 0
            test = 0;
            new_rays = [];
            Aglass = [];
            Asub = [];
            R = [];

            % Loops through all rays that are produced by the path of the
            % initial ray
            for k = 1:numBeams

                % Check to see if that ray's intensity is above the attenuation
                % threshold
                if rays(k,7) > thresh

                    % If the ray is above the threshold then figure out where
                    % it's next intersection is
                    [int_point,n,scenario,intensity,AglassNew] = ray_trace3(rays(k,1:3),rays(k,4:6),r1,r2,rays(k,7),rays(k,8),kappa,s,h);
                    %Keep track of rays that are leaving out of the top
                    %(point, direction, intensity)
                   if scenario == 6
                    new_scattered = [int_point,rays(k,4:6),intensity];
                    Scattered = [Scattered;new_scattered];
                    intensity;
                    scatInt = scatInt+intensity;
                                       
                   end
                    % Plot ray segment if debugging  %NEEDS TO BE
                    % DE-CELLIFIED
                      if exist('dbg','var') && dbg 
                            hold on
%                             pause
        %                   delete(rayPlot);
                            rayPlot = plot3([rays{k,1}(1), int_point(1)],[rays{k,1}(2), int_point(2)],[rays{k,1}(3), int_point(3)]);
        %                     orig = rays{k,1}
        %                     dir = rays{k,2}
        %                     int_point
        %                     n  
        %                     scenario
        %                     intensity
        %                     pause

                        end

                    % Eliminate ray segments incident on seam where hemis meet
                    % substrate
                    atSub = abs(int_point(3)) < .0000000000001;

                    atInnerHemi = abs(norm(int_point - [0 0 0]) - r1) <.00000000001;

                    atOuterHemi = abs(norm(int_point - [0 0 0]) - r2) <.00000000001;

                    if atSub && (atInnerHemi || atOuterHemi) 
                        corner1 = corner1 + 1;
                         Aglass(k) =  AglassNew;
                        R(k) =  0;
                         Asub(k) =  0;
                        lostInten = lostInten+ intensity;
                        continue

                    end

                    % Eliminate ray segments that tangentially intersect either
                    % of the hemispheres
                    atSphere = scenario == 2 || scenario == 4;
                    if abs(dot(rays(k,4:6),n)) > -.0000000001 && abs(dot(rays(k,4:6),n)) < .0000000001 && atSphere
                        tangent = tangent + 1;
                        continue
                    end


                    [cast_rays,RNew,AsubNew] = ray_cast3(int_point,rays(k,4:6),n,scenario,intensity,n1,m2(1),m2(2),n3,m4(1),m4(2),thetaCrit,s,r1,r2);

                    % Add the new rays generated if any to the end of the array
                    new_rays = [new_rays;cast_rays];

                    % Gather any of the property contributions
                    Aglass(k) =  AglassNew;
                    R(k) =  RNew;
                     Asub(k) =  AsubNew;
                    
                else
                    discardedInt = discardedInt + rays(k,7);

                end

    %             Aglass = sum(Aglass);
    %             R = sum(R);
    %             Asub = sum(Asub);

            end
                        Aglass = sum(Aglass);
                R = sum(R);
                Asub = sum(Asub);

           rays = new_rays;
           [numBeams,~] = size(rays);
           segments = segments + numBeams;


            AglassTot(j) = Aglass;
           RTot(j) = R;
            AsubTot(j) = Asub; 
           j = j+1; 

    end
     AglassTot = sum(AglassTot);
    RTot = sum(RTot);% + RFirstHit;
     AsubTot = sum(AsubTot);

     check1 = AglassTot + RTot + AsubTot;
    
    % Adaptive error control
    % If more than one percent of intensity lost, run the beam on a lower
    % atennuation threshold, if less, error is acceptable and move onto
    % next beam
%     if check1 < .95 && check1 > .94
%         thresh = thresh*(.5);
%         count = count+1;
% %         ii = ii + 1;
%     elseif check1 < .94 && check1 >.93
%           thresh = thresh*(.4);
%           count2 = count2 +1;
%     elseif check1 <.93
%         thresh = thresh*(.25);
%           count3 = count3 +1;
%     else
%         thresh = .00075;
        AglassTot1(ii) = AglassTot;
        RTot1(ii) = RTot;
        AsubTot1(ii) = AsubTot;
        ii = ii+1;

%      end
%     
%     AglassTot1(ii) = AglassTot;
%     RTot1(ii) = RTot;
%     AsubTot1(ii) = AsubTot;
    

 
end

 AglassTotFinal = sum(AglassTot1)/N;
RTotFinal = sum(RTot1)/N;
 AsubTotFinal = sum(AsubTot1)/N;

 check2 = AglassTotFinal + RTotFinal + AsubTotFinal;
 check3 = AglassTotFinal + scatInt/N + AsubTotFinal;
count = [count count2 count3];


end
