clc; clear all

r2 = [.005];%[.0005 .001 .002 .004 .008];
wavelength =.6E-6;%.55E-6; %m
% del = .25; sig = .99;
% tau = [.1 .2 .3 .4];
%  m4 = [4.084907216, 0.036329897];
% NST = 1-del;
% s = ((2*r2)/sig)
% incident_polar =[0  20 30 40 50 60 65 70 75 80 88];
% incident_azimuth = 0;
% res = 45;
% 
% for i = 1:4
%     k2 = wavelength*tau(i)/(4*pi*r2*del);
%     m2 = [1.535 k2]
%     for j = 1:11
%      [AglassTotFinal(i,j), AsubTotFinal(i,j),RTotFinal(i,j),check2(i,j),N,corner1,tangent,lostInten,segments,discardedInt,count,sideHit(i,j),TIRCountInner(i,j),path_lengthTot(i,j),neighHit(i,j),~,~,~,~,subHit(i,j)] = CubicCellHemiArrayProp(r2,NST,s,wavelength,m2,m4,incident_polar(j),incident_azimuth,res);
%     end
% end
%  m2 = [1.525,0.3E-8];  %FOR OPTICALLY THIN TAU = .00001
%    m2 = [1.525,3.5E-7]; %FOR OPTICALLY THIN TAU = .001
% m2 = [1.525,3.5E-5]; %FOR MODERATE TAU = .1
% m2 = [1.525,3.5E-5]
%     m2 = [1.525,3.8E-3]; %FOR THICK TAU = 10
%   m4 = [4.084907216, 0.036329897];
  m4 = [3.947296 .025826];
% m4 = [3.4261 1.99E-7];%[3.947296, .025826]
tau = .001;%[.001]% 0.2 0.3 0.4 0.5];%.001;%[.001 .5 10];


% wavelength = 3.5E-6;
% m2 = [1.474 1.065E-4];
% m4 = [3.432570007	1.00984E-08];

% wavelength = 5E-6;
% m2 = [1.397	3.000E-03]
% 
% m4 = [3.4261	0.000000199]


del = [.01 .25 .5 .75 .99];%[.01 .1 .2  .3  .4 .5 .6 .7 .8 .9 .99];%[.4 .42 .44 .45 .46 .48 .5];[0.99]%;
sig = .99;


% r2 = tau./(kappa.*del)
s = ((2.*r2)./sig)
r1 = -1.*(del-1).*r2

t = r2-r1;


NST = 1-del;

incident_polar =[0 30 40 50 60 65 70 72.5 75 77.5 80 82.5 85 88 89 89.5 89.7 89.9];%[0  20 30 40 50 60 65
incident_azimuth = [0];%[0 10 20 30 40 45 ];%[45];
res = 45;

%  [AglassTotFinal, AsubTotFinal,RTotFinal,check2,N,corner1,tangent,lostInten,segments,discardedInt,count,sideHit,TIRCountInner,path_lengthTot,neighHit,AglassPB,AsubPB,RPB,points] = CubicCellHemiArrayProp(r2,NST,s,wavelength,m2,m4,incident_polar,incident_azimuth,res);
len1 = length(del); len2 = length(del); len3 = length(incident_polar);
% m2 = [1.397, .003];
m2 = [1.523, 4.548E-7];
for i = 1:len3
        
%         planarGlass(i) = reflectivity(incident_polar(i),1,m2(1),0)
            planarSil(i) = reflectivity(incident_polar(i),1,m4(1),m4(2))
            planarAlR(i) = reflectivity(incident_polar(i),1,1.2,7.26);%8.67, 48.6);%1.2769,7.0328)
            planarGlassR(i) = reflectivity(incident_polar(i),1,m2(1), m2(2));
        
        for k = 1:len1
            tic
            k2 = (tau.*wavelength)./(4*pi*r2*del(k));
%             m2 = [1.397, .003];%[1.523 k2];
%             m2 = [1.523, 4.548E-7];
                m2 = [1.523 k2];
          
            
%               figure
            [AglassTotFinal(i,k), AsubTotFinal(i,k),RTotFinal(i,k),check2(i,k),N,corner1,tangent,lostInten,segments,discardedInt,count,sideHit(i,k),TIRCountInner(i,k),path_lengthTot(i,k),neighHit(i,k),~,~,~,~,subHit(i,k)] = CubicCellHemiArrayProp(r2,NST(k),s,wavelength,m2,m4,incident_polar(i),incident_azimuth,res);
%             view(0,0)
%             check2(i,j,k)
%             sideHit(i,j,k)
%             TIRCountInner(i,j,k)
            toc
        end
        

    i
end


% %%Constant del,varied tau
% del = .5;
% NST = 1-del;
% sig = .99;
% tau = [.001 .002 .004 .008 .016];
% r2 = tau./(kappa*del);
% s = ((2.*r2)./sig);
% incident_polar = [0 20 40 50 60 70 75 80 88];
% incident_azimuth = 0;
% res=45;
% len1 = length(tau); len3 = length(incident_polar);
% 
% for i = 1:len3
%     
%     for k = 1:len1
%         
%         tic
%        [AglassTotFinal(i,k), AsubTotFinal(i,k),RTotFinal(i,k),check2(i,k),N,corner1,tangent,lostInten,segments,discardedInt,count,sideHit(i,k),TIRCountInner(i,k)] = CubicCellHemiArrayProp(r2(k),NST,s(k),wavelength,m2,m4,incident_polar(i),incident_azimuth,res);
% 
%        toc
%         
%     end
%     
%     i
% end