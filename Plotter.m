% %Silicon and Glass for Solar
figure
plot(incident_polar, planarAlR,'LineWidth',2)
hold on
plot(incident_polar, planarSil,'LineWidth',2)
plot(incident_polar, RTotFinal(:,5),'LineWidth',2)
x = xlabel('\theta_{polar}','fontweight','bold','fontsize',14)
y = ylabel('Reflectance','fontweight','bold','fontsize',14)
legend('Planar Aluminum ','Planar Silicon','Modified Silicon \delta = 0.99','FontWeight','bold');

figure
plot(incident_polar, 1-planarAlR,'LineWidth',2)
hold on
plot(incident_polar, 1-planarSil,'LineWidth',2)
plot(incident_polar, AsubTotFinal(:,5),'LineWidth',2)
x = xlabel('\theta_{polar}','fontweight','bold','fontsize',14)
y = ylabel('Absorptance','fontweight','bold','fontsize',14)
legend('Planar Aluminum ','Planar Silicon','Modified Silicon (A_{sub} Only)\delta = 0.99','FontWeight','bold');

% %Glass and Glass for Windows
% figure
% plot(incident_polar, planarGlassR,'LineWidth',2)
% hold on
% plot(incident_polar, planarAlR,'LineWidth',2)
% plot(incident_polar, RTotFinal(:,5),'LineWidth',2)
% x = xlabel('\theta_{polar}','fontweight','bold','fontsize',14)
% y = ylabel('Reflectance','fontweight','bold','fontsize',14)
% legend('Planar Glass ','Planar Aluminum','Modified Glass \delta = 0.99','FontWeight','bold');
% 
% figure
% plot(incident_polar, 1-planarAlR,'LineWidth',2)
% hold on
% plot(incident_polar, 1-planarGlassR,'LineWidth',2)
% plot(incident_polar, AsubTotFinal(:,5),'LineWidth',2)
% x = xlabel('\theta_{polar}','fontweight','bold','fontsize',14)
% y = ylabel('Absorptance','fontweight','bold','fontsize',14)
% legend('Planar Aluminum ','Planar Glass','Modified Glass (A_{sub} Only) \delta = 0.99','FontWeight','bold');
