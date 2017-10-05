clear all; close all;clc;
fig_loc = 'C:\Users\Matt\Desktop\wc171002\figures\';

set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesfontsize',12);
fig_nompos = [100 100 550 400];
%% vg+
fig_name = [fig_loc,'vg_pl'];
N = 800;
Plin = linspace(-0.8,2.5,N);
Qlin = linspace(-1.3,1.3,N);

[P,Q] = meshgrid(Plin,Qlin);

V0 = 1.1;

STB_det = ((V0.^4)./4) + ((V0.^2).*P) - (Q.^2);

Vg2pl = ((V0.^2)/2) + P + sqrt(STB_det);

DNE = find(STB_det<=0); %does not exist
Vg2pl(DNE) = NaN;

fig = figure('Color','White','Position',fig_nompos);
contourf(P,Q,sqrt(Vg2pl),(0.5:0.1:2.4)); hold on;
contour(P,Q,sqrt(Vg2pl),(1.1:0.001:1.101),'r','Linewidth',2);
contour(P,Q,sqrt(Vg2pl),(0.9:0.002:0.902),'g','Linewidth',2);


Q_bdry = sqrt(((V0.^4)./4) + ((V0.^2).*Plin));
stb_det = ((V0.^4)./4) + ((V0.^2).*Plin);
dne = find(stb_det<=0);
Q_bdry(dne) = NaN;

plot(Plin,Q_bdry,'k')
plot(Plin,-Q_bdry,'k')

xlabel('$\tilde{P_{g}}$')
ylabel('$\tilde{Q_{g}}$')
axis('equal');
axis([min(Plin) max(Plin) min(Qlin) max(Qlin)]);

lgnd = legend('$|V_{g}|_{+}$','$V_{+}$','$V_{-}$','Location','NorthWest');
set(lgnd,'FontSize',12,'Interpreter','Latex');

% export_fig(fig,fig_name);
% export_fig(fig,strcat(fig_name,'.pdf'),'-dpdf');
%% vg+/-
fig_name = [fig_loc,'vg_pl_mn'];
Vg2mn = ((V0.^2)/2) + P - sqrt(STB_det);
Vg2mn(DNE) = NaN;
fig = figure('Color','White','Position',fig_nompos);
[C2,H2] = contourf(P,Q,sqrt(Vg2mn),(0:0.1:1)); hold on;
[C1,H1] = contour(P,Q,sqrt(Vg2pl),(0.5:0.1:2.4)); hold on;
contour(P,Q,sqrt(Vg2pl),(1.1:0.001:1.101),'r','Linewidth',2);
contour(P,Q,sqrt(Vg2pl),(0.9:0.002:0.902),'g','Linewidth',2);
contour(P,Q,sqrt(Vg2mn),(1.1:0.001:1.101),'r','Linewidth',2);
contour(P,Q,sqrt(Vg2mn),(0.9:0.002:0.902),'g','Linewidth',2);

plot(Plin,Q_bdry,'k')
plot(Plin,-Q_bdry,'k')

xlabel('$\tilde{P_{g}}$')
ylabel('$\tilde{Q_{g}}$')
axis('equal');
axis([min(Plin) max(Plin) min(Qlin) max(Qlin)]);

lgnd = legend('$|V_{g}|_{+}$','$|V_{g}|_{-}$','$V_{+}$','$V_{-}$','Location','NorthWest');
set(lgnd,'FontSize',12,'Interpreter','Latex');

% export_fig(gcf,fig_name);
% export_fig(gcf,strcat(fig_name,'.pdf'),'-dpdf');
%% vg_slt
% fig_name = strcat(fig_loc,'/vg_slt');
% N = 200;
% Plin = linspace(-0.8,1.1,N);
% Qlin = linspace(-1.0,1.0,N);
% 
% [P,Q] = meshgrid(Plin,Qlin);
% 
% V0 = 1.1;
% 
% STB_det = ((V0.^4)./4) + ((V0.^2).*P) - (Q.^2);
% 
% Vg2 = ((V0.^2)/2) + P + sqrt(STB_det);
% Sl = ((V0.^2)/2) + P - sqrt(STB_det);
% 
% 
% DNE = find(STB_det<=0); %does not exist
% Vg2(DNE) = NaN;
% Sl(DNE) = NaN;
% 
% 
% fig = figure('Color','White','Position',fig_nompos);
% [C2,H2] = contourf(P,Q,Sl,(0:0.1:1)); hold on;
% [C1,H1] = contour(P,Q,sqrt(Vg2),(0.5:0.1:1.8),'LineWidth',2); hold on;
% %contour(P,Q,sqrt(Vg2),(0.9:0.2:1.1),'r');
% 
% Q_bdry = sqrt(((V0.^4)./4) + ((V0.^2).*Plin));
% stb_det = ((V0.^4)./4) + ((V0.^2).*Plin);
% dne = find(stb_det<=0);
% Q_bdry(dne) = NaN;
% 
% plot(Plin,Q_bdry,'k')
% plot(Plin,-Q_bdry,'k')
% 
% xlabel('$\tilde{P_{g}}$')
% ylabel('$\tilde{Q_{g}}$')
% title('$|V_{g}|$ and $\tilde{S}_{l}$ against ($\tilde{P_{g}}$,$\tilde{Q_{g}}$)');
% axis('equal'); grid on;
% axis([min(Plin) max(Plin) min(Qlin) max(Qlin)]);
% 
% lgnd = legend([H1,H2],'$|V_{g}|$','$\tilde{S}_{l}$','Location','NorthWest');
% set(lgnd,'FontSize',12,'Interpreter','Latex');
% 
% % export_fig(gcf,fig_name);
% % export_fig(gcf,strcat(fig_name,'.pdf'),'-dpdf');
