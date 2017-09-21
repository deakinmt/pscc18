close all; clear all; clc;
cd('C:\Users\chri3793\Documents\MATLAB\DPhil\pscc18_mtlb')
addpath('pscc_fncs');
fig_loc = [pwd,'\figures\'];
%%
cd('C:\Users\chri3793\Documents\MATLAB\DPhil\pscc18_mtlb')
% INPUTS TO DETERMINE BEHAVIOUR -------------------
F.Vp = 1.06; %pu
F.Ip = 180;
% F.pg_ssc = linspace(-0.001,0.2,20);
% F.qg_ssc = linspace(-0.4,0.03,20);
F.pg_ssc = linspace(-0.001,0.2,50);
F.qg_ssc = linspace(-0.4,0.03,100);

% remain the same:
F.SRC = 'SOURCEBUS';

F.filename = '\opendss_models\34Bus\ieee34Mod1_fxd';
F.elmpwr_fnm = [pwd,'\opendss_models\34Bus\ieee34-1_EXP_ElemPowers.csv'];
F.feeder = '34bus'; % for calc_ww and T_ratios
F.tr_buses = {{'888','832'};{'832','852'};{'850','814'};{'800','SOURCEBUS'}};
F.fig_nompos = [200 150 450 340];

%%
F.NUT = '834';
pl_options={'imax'};
[ RR,figs ] = run_pscc_feeder( F,pl_options );

%%
NUTs = {'812','828','834','848'};
pl_options={'imax'};

RR = cell(size(NUTs)); figs = RR;
for i = 1:numel(NUTs)
    F.NUT = NUTs{i};
    [RR{i},figs{i}] = run_pscc_feeder(F,pl_options);
end

%%
fig = figure('Color','White');
p1 = zeros(size(RR)); p2 = zeros(size(RR)); 
p3 = zeros(size(RR)); p4 = zeros(size(RR));
for i = 1:numel(RR)
    R = RR{i};
    
    subplot(121)
    p1(i) = plot(R.S0(1),'o'); hold on; grid on;
    p2(i) = plot(R.S0(2:end),'x-');
    xlabel('Re(S0)'); ylabel('Im(S0)');
    axis equal;
    
    subplot(122)
    p3(i) = plot(R.Sgen(1),'o'); hold on; grid on;
    p4(i) = plot(R.Sgen(2:end),'x-');
    xlabel('Re(Sgen)'); ylabel('Im(Sgen)');
    axis equal;
    
end
subplot(121); 
legend([p1,p2],'812','828','834','848','812','828','834','848');
title('Ieee 34 nom tap, S0 pred v. meas.');
subplot(122); legend([p3,p4],'812','828','834','848','812','828','834','848');
title('Ieee 34 nom tap, Sg pred v. meas.');

figname = [fig_loc,'ieee34_4tests_nom'];
% export_fig(fig,figname);
% export_fig(fig,[figname,'.pdf'],'-dpdf');



