clc, clear all
%% FIRST PLOT %%%%%%%%%%%
%Input parameters
wavelength = .55E-6; %m
m2 = [1.532, 1.925E-7]; %1.532%1.925E-7
m4 = [4.084907216, 0.036329897];%4.084907216

NST = .9;


r2 = .005; %m

s =(1.5 *(r2*2));

h = r2+.0001;%m

incident_polar = 30;

incident_azimuth =0;


%Call the simulation
[AGlass, Asub, RTotFinal,check2,Scattered,N,~,~,~,~,discardedInt] = CubicCellHemiArray(r2,NST,s,wavelength,m2,m4,incident_polar,incident_azimuth,7);
check2
%Obtain the direction of all the beams leaving the system through top of
%unit cell
n_cart = zeros(3,length(Scattered));
i_beam = zeros(1,length(Scattered));
for ii = 1:length(Scattered)
    n_cart(:,ii) = Scattered(ii,4:6);
    i_beam(:,ii) = Scattered(ii,7);
end

n_polr(2,:) = acos(n_cart(3,:));
n_polr(1,:) = atan2(n_cart(2,:),n_cart(1,:));%acos(n_cart(1,:)./sin(n_polr(2,:)));

n_bins = 30;
d = pi/n_bins;
% Turn negative angles into positive ones
for k = 1:length(n_polr)
    if n_polr(1,k) < 0
            n_polr(1,k) = n_polr(1,k) + (2*pi);
    end
    
    if n_polr(1,k) >=0 && n_polr(1,k) < d
        
            n_polr(1,k) = n_polr(1,k) + (2*pi);
    end
end


%% Characterize the bins

dpsi = (2*pi)/n_bins;

omega = (2*pi)/((n_bins*n_bins)- n_bins + 1); %sr in each bin

centrPolr1 = acos(1-(omega/(2*pi))); %polar angle in central bin



a_bins1 = linspace(0+d,2*pi + d,n_bins+1);

% for k = 1:length(a_bins1)
%     
%     if a_bins1(k) > (2*pi)
%         
%         a_bins1(k) = a_bins1(k) - (2*pi);
%     end
% end

% a_bins1 = linspace(0,2*pi,n_bins+1);

p_bins1 = zeros(1,n_bins);

p_bins1(1) = 0;

p_bins1(2) = centrPolr1;%acos(cos(centrPolr1) - (omega/dpsi));

for i = 3:n_bins+1
    
   p_bins1(i) = acos(cos(p_bins1(i-1)) - (omega/dpsi));
   
end

% To check the bin boundaries in degrees
rad2deg(p_bins1);
rad2deg(a_bins1);

% To check if the bins contain the same solid angle
elementSolAng = dpsi * (cos(p_bins1(4)) - cos(p_bins1(5))); %sr

centerSolAng = (2*pi) * (cos(0) - cos(centrPolr1)); %sr

% Generate the polar plot
ray = ones(size(a_bins1));
figure

% subplot(1,2,1)
polar(a_bins1,p_bins1(end)*ray)
hold on

%% Plot the central element
% i_ptch = zeros(1,n_bins*n_bins);
ii = 1;
x_ptchCenter = zeros(n_bins,1);
y_ptchCenter = zeros(n_bins,1);

for kk = 1:n_bins
    %sums up intensity in all azimuthal bins that fall within the first
    %polar bin
    i_this = sum(i_beam(n_polr(2,:) >= p_bins1(ii) &...
    n_polr(2,:) < p_bins1(ii+1)));

    
    x_ptchCenter(kk) = cos(a_bins1(kk  ))*centrPolr1;
    y_ptchCenter(kk) = sin(a_bins1(kk  ))*centrPolr1;
          
end

%Divide by total number of beams
cen_int = i_this/N;

%Divide by solid angle in each bin
cen_intPerSA = cen_int/omega;
hold on
p = patch(x_ptchCenter,y_ptchCenter,cen_intPerSA);
p.LineStyle = 'none';



%% Populate remaining bins in polar plot
nn = 0;
for jj = 1:n_bins
    for ii = 2:n_bins

        
        i_this = sum(i_beam(n_polr(2,:)>=p_bins1(ii)  &...
                            n_polr(2,:)<p_bins1(ii+1)&...
                            n_polr(1,:)>=a_bins1(jj)  &...
                            n_polr(1,:)<a_bins1(jj+1)));
         i_bins(ii,jj) =  i_this;
        
        nn = nn+1;
        
        x1 = cos(a_bins1(jj  ))*p_bins1(ii  );
        x2 = cos(a_bins1(jj+1))*p_bins1(ii  );
        x3 = cos(a_bins1(jj+1))*p_bins1(ii+1);
        x4 = cos(a_bins1(jj  ))*p_bins1(ii+1);
        y1 = sin(a_bins1(jj  ))*p_bins1(ii  );
        y2 = sin(a_bins1(jj+1))*p_bins1(ii  );
        y3 = sin(a_bins1(jj+1))*p_bins1(ii+1);
        y4 = sin(a_bins1(jj  ))*p_bins1(ii+1);
        
        x_ptch(:,nn) = [x1 x2 x3 x4];
        y_ptch(:,nn) = [y1 y2 y3 y4];
        i_ptch(nn) = i_this;
        
    end
end

% Divide by total number of beams
i_ptch = i_ptch./N;

% Divide by solid angle in each bin
i_ptchPerSA = i_ptch/omega;

%% Plot the remaining bins

p=patch(x_ptch,y_ptch,i_ptchPerSA);

% cb1 = colorbar('Location','westoutside');
set(gca,'ColorScale','log');
% t = findall(gcf,'type','text');
% % delete the text objects
% delete(t);
p.LineStyle = 'none';

b1 = min(i_ptchPerSA); b2 = max(i_ptchPerSA);


% caxis([10^-1 ,  max(i_ptchPerSA)]);


check = RTotFinal - sum(i_ptch) - cen_int

max([i_ptchPerSA,cen_intPerSA]);


% %% SECOND PLOT %%%%%%%%%%%
% n_polr = [];
% %Input parameters
% wavelength = .55E-6; %m
% m2 = [1.532, 1000]; %1.532%1.925E-7
% m4 = [4.084907216, 0.036329897];%4.084907216
% 
% NST = .8;
% 
% 
% r2 = .005; %m
% 
% s =(1.01 *(r2*2));
% 
% h = r2+.0001;%m
% 
% incident_polar = 0;
% 
% incident_azimuth =0;
% 
% 
% %Call the simulation
% [AglassTotFinal, AsubTotFinal,RTotFinal,check2,Scattered,N] = CubicCellHemi(r2,NST,s,wavelength,m2,m4,incident_polar,incident_azimuth,125);
% 
% %Obtain the direction of all the beams leaving the system through top of
% %unit cell
% n_cart = zeros(3,length(Scattered));
% i_beam = zeros(1,length(Scattered));
% for ii = 1:length(Scattered)
%     n_cart(:,ii) = Scattered{ii,2};
%     i_beam(:,ii) = Scattered{ii,3};
% end
% 
% n_polr(2,:) = acos(n_cart(3,:));
% n_polr(1,:) = atan2(n_cart(2,:),n_cart(1,:));%acos(n_cart(1,:)./sin(n_polr(2,:)));
% 
% %n_bins = 50;
% d = pi/n_bins;
% % Turn negative angles into positive ones
% for k = 1:length(n_polr)
%     if n_polr(1,k) < 0
%             n_polr(1,k) = n_polr(1,k) + (2*pi);
%     end
%     
%     if n_polr(1,k) >=0 && n_polr(1,k) < d
%         
%             n_polr(1,k) = n_polr(1,k) + (2*pi);
%     end
% end
% 
% 
% %% Characterize the bins
% 
% dpsi = (2*pi)/n_bins;
% 
% omega = (2*pi)/((n_bins*n_bins)- n_bins + 1); %sr in each bin
% 
% centrPolr1 = acos(1-(omega/(2*pi))); %polar angle in central bin
% 
% 
% 
% a_bins1 = linspace(0+d,2*pi + d,n_bins+1);
% 
% % for k = 1:length(a_bins1)
% %     
% %     if a_bins1(k) > (2*pi)
% %         
% %         a_bins1(k) = a_bins1(k) - (2*pi);
% %     end
% % end
% 
% % a_bins1 = linspace(0,2*pi,n_bins+1);
% 
% p_bins1 = zeros(1,n_bins);
% 
% p_bins1(1) = 0;
% 
% p_bins1(2) = centrPolr1;%acos(cos(centrPolr1) - (omega/dpsi));
% 
% for i = 3:n_bins+1
%     
%    p_bins1(i) = acos(cos(p_bins1(i-1)) - (omega/dpsi));
%    
% end
% 
% % To check the bin boundaries in degrees
% rad2deg(p_bins1);
% rad2deg(a_bins1);
% 
% % To check if the bins contain the same solid angle
% elementSolAng = dpsi * (cos(p_bins1(4)) - cos(p_bins1(5))); %sr
% 
% centerSolAng = (2*pi) * (cos(0) - cos(centrPolr1)); %sr
% 
% % Generate the polar plot
% ray = ones(size(a_bins1));
% % figure
% 
% subplot(1,2,2)
% polar(a_bins1,p_bins1(end)*ray)
% hold on
% 
% %% Plot the central element
% % i_ptch = zeros(1,n_bins*n_bins);
% ii = 1;
% x_ptchCenter = zeros(n_bins,1);
% y_ptchCenter = zeros(n_bins,1);
% 
% for kk = 1:n_bins
%     %sums up intensity in all azimuthal bins that fall within the first
%     %polar bin
%     i_this = sum(i_beam(n_polr(2,:) >= p_bins1(ii) &...
%     n_polr(2,:) < p_bins1(ii+1)));
% 
%     
%     x_ptchCenter(kk) = cos(a_bins1(kk  ))*centrPolr1;
%     y_ptchCenter(kk) = sin(a_bins1(kk  ))*centrPolr1;
%           
% end
% 
% %Divide by total number of beams
% cen_int = i_this/N;
% 
% %Divide by solid angle in each bin
% cen_intPerSA = cen_int/omega;
% hold on
% p1 = patch(x_ptchCenter,y_ptchCenter,cen_intPerSA);
% p1.LineStyle = 'none';
% 
% 
% 
% %% Populate remaining bins in polar plot
% nn = 0;
% for jj = 1:n_bins
%     for ii = 2:n_bins
% 
%         
%         i_this = sum(i_beam(n_polr(2,:)>=p_bins1(ii)  &...
%                             n_polr(2,:)<p_bins1(ii+1)&...
%                             n_polr(1,:)>=a_bins1(jj)  &...
%                             n_polr(1,:)<a_bins1(jj+1)));
%          i_bins(ii,jj) =  i_this;
%         
%         nn = nn+1;
%         
%         x1 = cos(a_bins1(jj  ))*p_bins1(ii  );
%         x2 = cos(a_bins1(jj+1))*p_bins1(ii  );
%         x3 = cos(a_bins1(jj+1))*p_bins1(ii+1);
%         x4 = cos(a_bins1(jj  ))*p_bins1(ii+1);
%         y1 = sin(a_bins1(jj  ))*p_bins1(ii  );
%         y2 = sin(a_bins1(jj+1))*p_bins1(ii  );
%         y3 = sin(a_bins1(jj+1))*p_bins1(ii+1);
%         y4 = sin(a_bins1(jj  ))*p_bins1(ii+1);
%         
%         x_ptch(:,nn) = [x1 x2 x3 x4];
%         y_ptch(:,nn) = [y1 y2 y3 y4];
%         i_ptch(nn) = i_this;
%         
%     end
% end
% 
% % Divide by total number of beams
% i_ptch = i_ptch./N;
% 
% % Divide by solid angle in each bin
% i_ptchPerSA = i_ptch/omega;
% 
% %% Plot the remaining bins
% 
% p1=patch(x_ptch,y_ptch,i_ptchPerSA);
% 
% cb2 = colorbar('Location','westoutside');
% set(gca,'ColorScale','log');
% % t = findall(gcf,'type','text');
% % % delete the text objects
% % delete(t);
% p1.LineStyle = 'none';
% 
% b3 = min(i_ptchPerSA); b4 = max(i_ptchPerSA);
% 
% if b1<=b3
%     c_min = b1;
% else
%     c_min = b3;
% end
% 
% if b2 >= b4
%     c_max = b2;
% else
%     c_max = b4;
% end
% 
% caxis([10^-3,c_max]);
% 
% subplot(1,2,1)
% 
% cb1 = colorbar('Location','eastoutside');
% 
% caxis([10^-3,c_max]);
% 
% set(gca,'ColorScale','log');
% 
% cb1.Visible = 'off';
% 
% % %b1 = min(i_ptchPerSA); b2 = max(i_ptchPerSA);
% % b = max(i_ptchPerSA);
% % 
% % if b > b2
% %     b2 = b;
% %     
% % end
% % 
% % caxis([0 ,  b2]);
% 
% 
% 
% 
% check = RTotFinal - sum(i_ptch) - cen_int
% 
% max([i_ptchPerSA,cen_intPerSA])
