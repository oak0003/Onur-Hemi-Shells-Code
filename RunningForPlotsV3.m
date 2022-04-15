clc; clear all

r2 = [.002];%[.0005 .001 .002 .004 .008];
wavelength = .55E-6; %m
  m2 = [1.525,3.5E-7]; %FOR OPTICALLY THIN TAU = .001
%  m2 = [1.525,3.5E-5]; %FOR MODERATE TAU = .1
 %m2 = [1.525,3.8E-3]; %FOR MODERATE TAU = 10
 m4 = [4.084907216, 0.036329897];



% wavelength = 3.5E-6;
% m2 = [1.474 1.065E-4];
% m4 = [3.432570007	1.00984E-08];

% wavelength = 5E-6;
% m2 = [1.397	3.000E-03]
% 
% m4 = [3.4261	0.000000199]


kappa = (4*pi*m2(2))/wavelength;

del =[.5];%[.01 .1 .2  .3  .4 .5 .6 .7 .8 .9 .99];
sig = .99;
tau = .001;

r2 = tau./(kappa.*del)
s = ((2.*r2)./sig)
r1 = -1.*(del-1).*r2

t = r2-r1;


NST = 1-del;

incident_polar = [80];%[40 50 70 75 85 88 89];% [0  20 30 40 50 60 65 70 75 80 88];%[0  20 30 40 50 60 65
incident_azimuth = 0;
res = 7;
len1 = length(del); len2 = length(del); len3 = length(incident_polar);

for i = 1:len3
        
        planarGlass(i) = reflectivity(incident_polar(i),1,m2(1),m2(2))
        
        for k = 1:len1
            tic
            figure
            [AglassTotFinal(i,k), AsubTotFinal(i,k),RTotFinal(i,k),check2(i,k),N,corner1,tangent,lostInten,segments,discardedInt,count,sideHit(i,k),TIRCountInner(i,k)] = CubicCellHemiArrayPropSingleBeam(r2(k),NST(k),s(k),wavelength,m2,m4,incident_polar(i),incident_azimuth,res);
%             view(0,0)
%             check2(i,j,k)
%             sideHit(i,j,k)
%             TIRCountInner(i,j,k)
%  view(0,0)
            toc
            hold on
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