close all; clear all;

set(0,'defaulttextinterpreter','latex');
set(0,'defaultaxesfontsize',14);

fig_pos = [300,200,400,340];

pg = linspace(-0.8,1.2,3e2);
qg = linspace(-1,1,3e2);
v0 = 1;

[Pg,Qg] = meshgrid(pg,qg);


dscr = ((v0^4)/4) + ((v0^2)*Pg) - (Qg.^2);

outs = 0./(dscr>0);

Vg = sqrt(Pg + (v0^2)/2 + sqrt(dscr)) + outs*(1 + 1i);
Sl = Pg + (v0^2)/2 - sqrt(dscr) + outs*(1 + 1i);
bdry = ((qg.^2)-((v0^4)/4))/(v0^2);
%
figure('Color','White','Position',fig_pos);

[~,plts{1}] = contourf(Pg,Qg,real(Sl)); hold on;
[~,plts{2}] = contour(Pg,Qg,real(Vg),'Linewidth',2);
plot(bdry,qg,'k');

axis equal;
xlabel('$\tilde{P}_{g}$');
ylabel('$\tilde{Q}_{g}$');
lgnd = legend([plts{2},plts{1}],'$V_{g}$\textsuperscript{[+ve]}','$\tilde{S}_{l}$\,\textsuperscript{[-ve]}');
set(lgnd,'Interpreter','Latex','Fontsize',14,'Location','NorthWest');

export_fig(gcf,'sg_slt_v1','-r300');