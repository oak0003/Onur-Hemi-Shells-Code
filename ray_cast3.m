% i_cas possibilities:

% 1: outside shell
% 2: outer hemi going out
% 3: outer hemi going in
% 4: inner hemi going in
% 5: inner hemi going out
% 6: substrate under shell
% 7: substrate within shell

function [cast_rays,R,Asub,TIRbool] = ray_cast3(int_point,dir,n,scenario,intensity,n1,n2,k2,n3,n4,k4,thetaCrit,span,r1,r2)
    
    % Scenario 1: Ray hits inside of outer hemi moving out
    % Reflected ray cast, refracted ray sent out of shell
    % possibility of TIR
    if scenario == 1
        %CONDITIONAL HERE TO CHECK FOR TIR, NO REFLECTANCE
        [reflected_beam,theta_1] = reflect(dir,n);
        
        if theta_1>thetaCrit
            rho = 1;
            R= 0;
            Asub =0;
            intensity_flect = intensity * rho;
            i_cas = 3; % reflected beam is at outer hemi going in
            cast_rays = [int_point,reflected_beam,intensity_flect,i_cas];
            TIRbool = 1;
        else
            [refracted_beam,theta_1,theta_2] = refract(n2,n1,0,dir,n);
            rho = reduced_fresnel(theta_1, theta_2);
            intensity_flect = intensity * rho;
            intensity_fract = intensity * (1-rho);
            %Reflected ray is at outer hemi going in, refracted beam is at
            %outer hemi going out
            cast_rays = [int_point,reflected_beam,intensity_flect,3;int_point,refracted_beam,intensity_fract,2];% {int_point,reflected_beam,intensity_flect};
            R = 0;% intensity_fract;
            Asub =0;
            TIRbool = 0;
            
        end

    
    % Scenario 4: Beam hits outside of inner hemi moving in
    % Cast reflected and refracted ray, possibility of TIR (in which case
    % no refracted ray
    elseif scenario == 4
        [reflected_beam,theta_1] = reflect(dir,n);
        
        %CONDITIONAL HERE TO CHECK FOR TIR
        if theta_1>thetaCrit
           rho = 1; %WILL BE REPLACED WITH ACTUAL RHO FORMULATION
           intensity_flect = intensity * rho;
           cast_rays = [int_point,reflected_beam,intensity_flect,5];
%            if intensity_flect > .01
           TIRbool = 1;          
%            else
%                TIRbool = 0;
%            end
            
        else
            [refracted_beam,theta_1,theta_2] = refract(n2,n3,0,dir,n);
            rho = reduced_fresnel(theta_1, theta_2);
%             rho = 0.5; %WILL BE REPLACED WITH ACTUAL RHO FORMULATION
            intensity_flect = intensity * rho;
            intensity_fract = intensity * (1-rho);
            %reflected beam is at inner hemi going out, refracted beam is
            %at inner hemi going in
            cast_rays = [int_point,reflected_beam,intensity_flect,5;int_point,refracted_beam,intensity_fract,4];
            TIRbool = 0;
            
        end
        
        R = 0;
        Asub = 0;
        
     
    % Scenario 3: Beam hits inside of inner hemi moving out
    % Cast reflected and refracted ray
    elseif scenario == 3
        [reflected_beam,theta_1] = reflect(dir,n);
        [refracted_beam,theta_1,theta_2] = refract(n3,n2,k2,dir,n);
        rho = reflectivity(theta_1,n3,n2,k2);
%         rho = 0.5; %WILL BE REPLACED WITH ACTUAL RHO FORMULATION
        intensity_flect = intensity * rho;
        intensity_fract = intensity * (1-rho);
        
        %reflected beam is inner hemi going in, refracted inner hemi moving
        %out
        cast_rays = [int_point,reflected_beam,intensity_flect,4;int_point,refracted_beam,intensity_fract,5];
        
        R = 0;
        Asub = 0;
        TIRbool = 0;
    % Scenario 5: Beam hits the planar substrate
    % Cast reflected ray, refracted ray's intensity is absorbed by
    % substrate
    elseif scenario == 5
        [reflected_beam,theta_1] = reflect(dir,n);
        rho =  reflectivity(theta_1,n3,n4,k4);
        intensity_flect = intensity *rho;
        
        
        %CODE NEEDED TO DETERMINE WHERE ON THE SUBSTRATE,1,6,or 7
%         -.000000001 < (norm([ix,iy]) - r2) &&  (norm([ix,iy]) - r2) < .000000001
        ix = int_point(1); iy = int_point(2); iz = int_point(3);
        if iz~= 0
            error('Something is wrong and the beam is not at the substrate');
        elseif norm([ix,iy]) > r2
            i_cas = 1;
        elseif norm([ix,iy]) < r1
            i_cas = 6;
        else
            i_cas = 7;
        end
       
        
        cast_rays = [int_point,reflected_beam,intensity_flect,i_cas];
        Asub = intensity - intensity_flect;
        R = 0;
        TIRbool = 0;
    % Scenario 2: Beam hits outside of outer hemi moving in, a
    % refracted ray and reflected ray are cast with no TIR possibility
    elseif scenario == 2
        [reflected_beam,theta_1] = reflect(dir,n);        
        [refracted_beam,theta_1,theta_2] = refract(n1,n2,k2,dir,n);
        rho =  reflectivity(theta_1,n1,n2,k2);
%         rho = reduced_fresnel(theta_1, theta_2);
        intensity_flect = intensity * rho;
        intensity_fract = intensity * (1-rho);
        cast_rays = [int_point,reflected_beam,intensity_flect,2;int_point,refracted_beam,intensity_fract,3];% {int_point,reflected_beam,intensity_flect};
        R = 0;% intensity_fract;
        Asub =0;
        TIRbool = 0;
    % Scenario 6: Beam hits the top plane and leaves the unit cell, no new
    % ray is cast
    elseif scenario == 6
        R = intensity;
        Asub = 0;
        cast_rays = [];
        TIRbool = 0;
    
    % Scenario 7: Beam hits a side plane of the unit cell and the
    % periodicity effects occur and a single translated ray is cast.
    elseif scenario ==7
        R = 0; Asub = 0; TIRbool = 0;
        
        % Translate the beam, x and y change but z does not
        ix = int_point(1);iy = int_point(2);iz = int_point(3);
        nx = n(1); ny = n(2); nz = n(3);
        ox = ix -(-nx)*span; oy = iy - (-ny)*span;
        
        orig_point = [ox,oy,iz];
        
        % If we periodize by translating the beam, the beam is still
        % outside the shell
        cast_rays =[orig_point,dir,intensity,1];
    
    end
    
   

end