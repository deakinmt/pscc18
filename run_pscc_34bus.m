close all; clear all; clc;
cd('C:\Users\chri3793\Documents\MATLAB\DPhil\pscc18_mtlb')
% cd('C:\Users\Matt\Documents\MATLAB\DPhil\pscc18_mtlb');
addpath('pscc_fncs');
% fig_loc = [pwd,'\figures\'];
fig_loc = 'C:\Users\chri3793\Documents\DPhil\pscc18\pscc18_paper\figures\';
% fig_loc = 'C:\Users\Matt\Documents\DPhil\pscc18\pscc18_paper\figures\';
%%
cd('C:\Users\chri3793\Documents\MATLAB\DPhil\pscc18_mtlb')
% cd('C:\Users\Matt\Documents\MATLAB\DPhil\pscc18_mtlb');
% INPUTS TO DETERMINE BEHAVIOUR -------------------
F.Vp = 1.06; %pu
F.Ip = 180; %A
F.Pp = 1; %pu
F.pg_ssc = linspace(-0.001,0.2,20);
F.qg_ssc = linspace(-0.4,0.03,20);
% F.pg_ssc = linspace(-0.001,0.2,50);
% F.qg_ssc = linspace(-0.4,0.03,100);

% remain the same:
F.SRC = 'SOURCEBUS';

F.filename = '\opendss_models\34Bus\ieee34Mod1_fxd';
F.feeder = '34bus'; % for calc_ww and T_ratios
F.tr_buses = {{'888','832'};{'832','852'};{'850','814'};{'800','SOURCEBUS'}};
F.fig_nompos = [200 150 450 340];

% %% First results figure: 834 bus feeder voltage
% figname1 = [fig_loc,'34bus_pgp0'];
% F.NUT = '834';
% F.pg_ssc = linspace(0.001,0.2,20);
% F.qg_ssc = linspace(-0.4,0.03,20);
% pl_options={'pgp0'};
% [ ~,figs ] = run_pscc_feeder( F,pl_options );
% 
% % export_fig(figs,figname);
% % export_fig(figs,[figname,'.pdf']);
% %% First results figure: 834 bus currents
% figname2 = [fig_loc,'34bus_imax'];
% F.NUT = '834';
% F.pg_ssc = linspace(1e-4,0.2,20);
% F.qg_ssc = linspace(-0.4,0.03,20);
% pl_options={'imax'};
% [ ~,figs ] = run_pscc_feeder( F,pl_options );
% % export_fig(figs,figname);
% % export_fig(figs,[figname,'.pdf']);
% %%
% figname3 = [fig_loc,'34bus_pgqgq0'];
% F.NUT = '834';
% F.pg_ssc = linspace(1e-4,0.2,20);
% F.qg_ssc = linspace(-0.4,0.03,20);
% pl_options={'pgqgq0'};
% [ ~,figs ] = run_pscc_feeder( F,pl_options );
% % export_fig(figs,figname);
% % export_fig(figs,[figname,'.pdf']);

%%
F.NUT = '834';
F.pg_ssc = linspace(1e-4,0.2,80);
F.qg_ssc = linspace(-0.4,0.03,400);
pl_options={'pgp0','pgqgq0','imax'}; %NB the order DOES matter!
[ ~,figs ] = run_pscc_feeder( F,pl_options );
% for i = 1:numel(pl_options)
%     figname = [fig_loc,'34bus_',pl_options{i}];
%     export_fig(figs(i),figname);
%     export_fig(figs(i),[figname,'.pdf']);
% end

%%
NUTs = {'812','828','834','848'};
F.pg_ssc = linspace(1e-4,0.2,80);
F.qg_ssc = linspace(-0.4,0.03,400);
pl_options=NaN;

RR = cell(size(NUTs)); figs = RR;
%%
for i = 4:numel(NUTs)
    F.NUT = NUTs{i};
    [RR{i},figs{i}] = run_pscc_feeder(F,pl_options);
end
%%
fig = figure('Color','White','Position',[100 150 550 450],'defaultaxesfontsize',12,'defaulttextinterpreter','latex');
figname = [fig_loc,'34bus_accuracy'];

for i=1:4
    rslt = [RR{i}.Pgen_prm RR{i}.Pgen_hat RR{i}.Psub_prm RR{i}.Psub_hat ]';

    subplot(2,2,i)
    bar(rslt);
    title(['Bus ',NUTs{i}]); grid on; 
    lgnd = legend('Msrd.','Estd.');
    set(lgnd,'Interpreter','Latex');


    
    Yticks = yticks; ytcks = cell(numel(yticks),1);
    for j = 1:numel(Yticks)
        ytcks{j} = num2str(Yticks(j));
    end
    xtcks = {'$P_{gen}''$','$\hat{P}_{gen}$','$P^{Sub \, \prime}_{0} $ ','$\hat{P}^{Sub}_{0}$'};
    format_ticks(gca,xtcks,[],[],[],[],[],[]);
    ylabel('$P$ (pu)');
    axis([0.5 4.5 0 round(max(max(rslt)))+0.5]);    
end

export_fig(fig,figname);
export_fig(fig,[figname,'.pdf']);

%% ERROR TABLE
clear input; clc;
ERR = zeros(numel(NUTs),5);
for i = 1:4
        RES = [RR{i}.Pgen_prm RR{i}.Pgen_hat RR{i}.Psub_prm RR{i}.Psub_hat ];
        ERR(i,1) = str2double(NUTs{i});
        ERR(i,2:end) = (RES(1,:) - RES(2,:));
end

input.data = ERR;
input.tableColLabels = {'Bus','$\epsilon (P_{gen}'')$','$\epsilon (\hat{P}_{gen})$',...
                        '$\epsilon (P_{0}^{Sub \, \prime})$ ','$\epsilon (\hat{P}_{0}^{Sub})$  '    };
input.dataFormat = {'%.0f',1,'%.2f',4}; % three digits precision for first two columns, one digit for the last                    
input.tableColumnAlignment = 'l';
input.booktabs = 1;
input.tableCaption = 'Marginal and thermal power limit error (pu)';
input.tableLabel = 'res_error';

err_table = latexTable(input);