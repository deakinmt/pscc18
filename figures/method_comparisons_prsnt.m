% method_compatisons.m is a script that plots a comparison between the
% marginal loss-induced power transmission limit and solution boundary and
% unity power factor control (Figure 2)


close all; clear all; clc;
fig_loc = 'C:\Users\Matt\Documents\DPhil\pscc18\pscc18_paper\figures';

set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',12);
fig_nompos = [100 100 550 320];
%%
V0 = 1;
Vg = 1.06;

theta_lo = 0.01;
theta_hi = pi/2 - 0.01;
theta = linspace(theta_lo,theta_hi,10000);
% theta = linspace(theta_lo,theta_hi,300);
lz = cot(theta);
RX = lz;
Z = cos(theta) + 1i*sin(theta);

lr = lz./sqrt(1 + (lz.^2));
lx = 1./sqrt(1 + (lz.^2));
lv = Vg./V0;

Pgt = (V0.^2).*lv.*(lv - lr); % we know it is the two negative solutions from the derivation
Qgt = -(V0.^2).*lv.*lx;

Sl = (Pgt.^2 + Qgt.^2)./(Vg.*(lr - 1i*lx));

P0pr = (V0.^2).*lr.*((lv.^2) - 1) - (lx.*(lz.*Pgt + Qgt)); %NB |Z| = 1!
Pgpr = lr.*Pgt - lx.*Qgt;
Qgpr = lx.*Pgt + lr.*Qgt;


% stability values:
P0stb = (V0.^2).*( (-0.5*lr) + (lx.*sqrt((lv.^2) - 0.25)) );
Pgstb = P0stb + (((V0.*lv).^2).*lr);
Slstb = (Vg.^2)./conj(Z);

% stab/volt values:
lz_prpr = 1./sqrt((2*lv)^2 - 1);

% Unity PF values:
lz_0 = sqrt((lv^2) - 1);

Pgupf = (V0.^2).^lv.*( (lr.*lv) - sqrt(1 - ((lx.*lv).^2) ) );
P0upf = Pgupf.*(1 - ((lr.*Pgupf)./(Vg.^2)));

Pgupf = Pgupf./(1-(lz<lz_0)); %NaN if doesn't exist
P0upf = P0upf./(1-(lz<lz_0));

%% mgtt_1dcomparison_2
fig_name = [fig_loc,'\mgtt_1dcomparison_2'];
fig = figure('Color','White','Position',fig_nompos); 


semilogx(lz,real(Pgupf),'b--'); hold on;
semilogx(lz,real(P0upf),'b');
semilogx(lz,Pgstb,'r--'); 
semilogx(lz,P0stb,'r');
semilogx(lz,Pgpr,'k--');
semilogx(lz,P0pr,'k');

grid on;
xlabel('$\lambda = R/X$');
ylabel('$P$ (p.u.)');

lgnd = legend('$P_{g}$\textsuperscript{UPF}','$P_{0}$\textsuperscript{UPF}',...
              '$P_{g}$\textsuperscript{Stb.}','$P_{0}$\textsuperscript{Stb.}',...
              '$P_{g}''$','$P_{0}''$');
set(lgnd,'Interpreter','latex','Location','SouthWest','FontSize',14);

% export_fig(gcf,fig_name);
% export_fig(gcf,[fig_name,'.pdf'],'-dpdf');
%% mgtt_1dcomparison_1
fig_name = [fig_loc,'\mgtt_1dcomparison_1'];
fig = figure('Color','White','Position',fig_nompos); 


semilogx(lz,real(Pgupf),'b--'); hold on;
semilogx(lz,real(P0upf),'b');
semilogx(lz,Pgstb,'r--'); 
semilogx(lz,P0stb,'r');
% semilogx(lz,Pgpr,'k--');
% semilogx(lz,P0pr,'k');

grid on;
xlabel('$\lambda = R/X$');
ylabel('$P$ (p.u.)');

lgnd = legend('$P_{g}$\textsuperscript{UPF}','$P_{0}$\textsuperscript{UPF}',...
              '$P_{g}$\textsuperscript{Stb.}','$P_{0}$\textsuperscript{Stb.}');
%               '$P_{g}$\textsuperscript{Mrg.}','$P_{0}$\textsuperscript{Mrg.}');
set(lgnd,'Interpreter','latex','Location','SouthWest','FontSize',14);
xs = axis;
% export_fig(gcf,fig_name);
% export_fig(gcf,[fig_name,'.pdf'],'-dpdf');
%% mgtt_1dcomparison_0
fig_name = [fig_loc,'\mgtt_1dcomparison_0'];
fig = figure('Color','White','Position',fig_nompos); 


semilogx(lz,real(Pgupf),'b--'); hold on;
semilogx(lz,real(P0upf),'b');
% semilogx(lz,Pgstb,'r--'); 
% semilogx(lz,P0stb,'r');
% semilogx(lz,Pgpr,'k--');
% semilogx(lz,P0pr,'k');

grid on;
xlabel('$\lambda = R/X$');
ylabel('$P$ (p.u.)');
axis(xs);
lgnd = legend('$P_{g}$\textsuperscript{UPF}','$P_{0}$\textsuperscript{UPF}');
%               '$P_{g}$\textsuperscript{Stb.}','$P_{0}$\textsuperscript{Stb.}',...
%               '$P_{g}$\textsuperscript{Mrg.}','$P_{0}$\textsuperscript{Mrg.}');
set(lgnd,'Interpreter','latex','Location','SouthWest','FontSize',14);

% export_fig(gcf,fig_name);
% export_fig(gcf,[fig_name,'.pdf'],'-dpdf');