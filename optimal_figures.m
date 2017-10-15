% This script creates the `measures of efficiency' in Figure 3 - thermal
% efficiency, and the power factors at generator and feeder head.

clear all; close all; 
fig_loc = 'figures';

set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',12);
fig_nompos = [100 100 550 320];
%% Calculate the power factors and 
Vg = 1.06;
v0_lin = linspace(0.9,1.2,300);
theta = linspace(0.003,(pi/2) - 0.0006,1000);

lz_lin = cot(theta); % lambda
lv_lin = Vg./v0_lin;
nu_lin = 1./lv_lin;

[lz,nu] = meshgrid(lz_lin,nu_lin);
lv = 1./nu;
V0 = Vg./lv;

lr = lz./sqrt(1 + (lz.^2));
lx = 1./sqrt(1 + (lz.^2));
Z = lr + 1i*lx;

Pgt = (V0.^2).*lv.*(lv - lr); % we know it is the two negative solutions from the derivation
Qgt = -(V0.^2).*lv.*lx;

P0pr = (V0.^2).*lr.*((lv.^2) - 1) - (lx.*(lz.*Pgt + Qgt)); %NB |Z| = 1!
Q0pr = (V0.^2).*lx.*(lv.^2 - (2.*lr.*lv));
Pgpr = lr.*Pgt - lx.*Qgt;
Qgpr = lx.*Pgt + lr.*Qgt;

eta = (P0pr./Pgpr).*(P0pr>=0) ... 
                                    + (P0pr<0)*-0.001; %efficiency
PF0pr = abs(P0pr)./sqrt(P0pr.^2 + Q0pr.^2);
PFgpr = abs(Pgpr)./sqrt(Pgpr.^2 + Qgpr.^2);


%% plot eta
figure('Color','White','Position',fig_nompos);
fig_name = [fig_loc,'\mgtt_eta'];
[c_con,~] = contour(lz,V0,eta,(-1:0.1:1)); % ,(0:0.1:1)
set(gca,'xscale','log');
clabel(c_con,[0 0.2 0.4 0.6 0.8 0.9]); grid on;
xlabel('$\lambda = R/X$ '); ylabel('$V_{0}$ (p.u.)');

lgnd = legend('$\frac{P_{0}''}{P_{g}''}$ ','Location','NorthWest');
set(lgnd,'FontSize',18,'Interpreter','Latex');

% export_fig(gcf,fig_name);
% export_fig(gcf,strcat(fig_name,'.pdf'),'-dpdf');
%% plot pf_g
figure('Color','White','Position',fig_nompos);
fig_name = [fig_loc,'\mgtt_pf_g'];
[c_con,~] = contour(lz,V0,PFgpr);
set(gca,'xscale','log');
clabel(c_con,[0.1 0.3 0.5 0.7 0.8 0.9]); grid on;

xlabel('$\lambda = R/X$ '); ylabel('$V_{0}$ (p.u.)');
lgnd = legend('$\frac{P_{g}''}{|S_{g}''|}$');
set(lgnd,'FontSize',18,'Interpreter','Latex','Location','NorthWest');

% export_fig(gcf,fig_name);
% export_fig(gcf,strcat(fig_name,'.pdf'),'-dpdf');
%% plot pf_0
figure('Color','White','Position',fig_nompos);
fig_name = [fig_loc,'\mgtt_pf_0'];
[c_con,~] = contour(lz,V0,PF0pr);
set(gca,'xscale','log');
clabel(c_con,[0.1 0.3 0.5 0.7 0.9]); grid on;

xlabel('$\lambda = R/X$ '); ylabel('$V_{0}$ (p.u.)');
lgnd = legend('$\frac{P_{0}''}{|S_{0}''|}$');
set(lgnd,'FontSize',18,'Interpreter','Latex','Location','NorthWest');

% export_fig(gcf,fig_name);
% export_fig(gcf,strcat(fig_name,'.pdf'),'-dpdf');