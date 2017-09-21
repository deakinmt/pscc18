function [ RR,figs ] = run_pscc_feeder( FF,pl_options )

% INPUTS TO DETERMINE BEHAVIOUR -------------------
Vp = FF.Vp; % we set for this network.
NUT = FF.NUT;
SRC = FF.SRC;
pg_ssc = FF.pg_ssc;
qg_ssc = FF.qg_ssc;

filename = FF.filename;
elmpwr_fnm = FF.elmpwr_fnm;
fig_nompos = FF.fig_nompos;
%% ------------------------------------------------
% Run the DSS
[~, DSSObj, DSSText] = DSSStartup;
DSSCircuit=DSSObj.ActiveCircuit; DSSSolution=DSSCircuit.Solution;
DSSText.command=['Compile (',pwd,filename,'.dss)'];
DSSSolution.Solve;

% DSSCircuit.SetActiveClass('PDElements');

Zbus = create_zbus(DSSCircuit);

T_ratio = 1;
V0 = DSSCircuit.Vsource.pu*T_ratio;

DSSCircuit.SetActiveBus(NUT);
vbN = sqrt(3)*DSSCircuit.ActiveBus.kVBase; % LL (@ NUT)
DSSCircuit.SetActiveBus(SRC);
vbS = sqrt(3)*DSSCircuit.ActiveBus.kVbase; % LL (@ SRC)
nTr0 = vbN/vbS;
DSSCircuit.Transformers.First;
sb = DSSCircuit.Transformers.kVA/10; %NB - this is wrong in ieee34mod1 script.
zbN = 1e3*(vbN^2)/sb; %ohms
ibN = sb*sqrt(1/3)/vbN; %A, per phase


YZNodeOrder = DSSCircuit.YNodeOrder;

Z_ohm  = find_node_Z_pscc( YZNodeOrder,Zbus,SRC,NUT,nTr0 );
Z = Z_ohm/zbN; %pu

Ssc = ((vbN^2)/abs(Z*zbN))*1e3; %kW

DSSObj.Text.Command='export elempowers';
CapQ = calc_capQ(elmpwr_fnm)/sb; %pu 
LdNom = (calc_load(elmpwr_fnm)/sb); %pu
Sload = CapQ + LdNom; %pu

% create generator at given node:
DSSText.Command=['new generator.gen phases=3 bus1=',NUT,'.1.2.3 kw=0 pf=0'];
% DSSText.Command='new monitor.genvi element=generator.gen terminal=1 mode=32 VIPolar=yes';
DSSText.Command='new monitor.genpq element=generator.gen terminal=1 mode=65 PPolar=no';

%% Run the DSS
[Pg, Qg] = meshgrid(Ssc*pg_ssc,Ssc*qg_ssc);
Pf = sign(Qg).*Pg./sqrt(Pg.^2 + Qg.^2);

Pl_mes_lin = zeros(size(Pg));
totlss = zeros(numel(Pg),2); totpwr = zeros(numel(Pg),2);
VmaxMes = zeros(numel(Pg),1); VminMes = zeros(numel(Pg),1);
ImaxMes = zeros(numel(Pg),1);

DSSCircuit.Monitors.ResetAll;
tic
for i = 1:numel(Pg)
    DSSText.Command=strcat('edit generator.gen kw=',num2str(Pg(i)));
    DSSText.Command=strcat('edit generator.gen pf=',num2str(Pf(i)));
    DSSSolution.Solve;
    DSSCircuit.Sample;
    
    [ I,~ ] = meas_pde_i( DSSCircuit );
    ImaxMes(i) = max(max(abs(I)));
    AllBusVmagPu = DSSCircuit.AllBusVmagPu; 
    VmaxMes(i) = max(AllBusVmagPu);
    VminMes(i) = min(AllBusVmagPu);
    totlss(i,:) = 1e-3*DSSObj.ActiveCircuit.Losses/sb; %pu (output in W)
    totpwr(i,:) = DSSObj.ActiveCircuit.TotalPower/sb; %pu
end
toc
% withdraw results
% DSSMon=DSSCircuit.Monitors; DSSMon.name='genvi';
% DSSText.Command='show monitor genvi'; % this has to be in for some reason?
% Vgen_lin = ExtractMonitorData(DSSMon,1,vbN*1e3/sqrt(3)); %phase 1
% Vgen_lin2 = ExtractMonitorData(DSSMon,2,vbN*1e3/sqrt(3)); %phase 2
% Vgen_lin3 = ExtractMonitorData(DSSMon,3,vbN*1e3/sqrt(3)); %phase 3
% Vgen1=reshape(Vgen_lin,size(Pg));
% Vgen2=reshape(Vgen_lin2,size(Pg));
% Vgen3=reshape(Vgen_lin3,size(Pg));

ImaxMat = reshape(ImaxMes,size(Pg));
VmaxMat = reshape(VmaxMes,size(Pg));
VminMat = reshape(VminMes,size(Pg));
% Vgen = max(cat(3,Vgen1,Vgen2,Vgen3),[],3); %assume voltage limit refers to a maximum phase LN voltage magnitude
% Vgen = reshape(Vgen_lin,size(Pg));
Vgen = VmaxMat;

DSSMon=DSSCircuit.Monitors; DSSMon.name='genpq';
DSSText.Command='show monitor genpq';
Pgen_lin = -ExtractMonitorData(DSSMon,1,1)/sb; %pu, +ve is generation
Qgen_lin = -ExtractMonitorData(DSSMon,2,1)/sb; %pu, +ve is generation

TotLss = reshape(totlss(:,1),size(Pg)) + 1i*reshape(totlss(:,2),size(Pg));% +ve implies losses
TotPwr = reshape(totpwr(:,1),size(Pg)) + 1i*reshape(totpwr(:,2),size(Pg)); %-ve implies load

% calculate maximum power export s.t. voltage constraint:
% NaN_outs = 0./(Vgen<Vp);
% Vg_in = Vgen + NaN_outs; % remove numbers where the voltage is too high
% Vg_in_mn = (Vg_in==max(Vg_in));
NaN_outs = 0./(VmaxMat<Vp);
Vg_in = VmaxMat + NaN_outs; % remove numbers where the voltage is too high
Vg_in_mn = (Vg_in==max(Vg_in));
max_pwr = max(real(TotPwr(Vg_in_mn)),[],2); % for each row, find the max power
Vmn = find(real(TotPwr)==repmat(max_pwr',numel(TotPwr)/numel(max_pwr),1)); %find the index of these points

Pgenmat=reshape(Pgen_lin,size(Pg));% - real(TotLd);
Qgenmat=reshape(Qgen_lin,size(Pg));% - imag(TotLd);
Sgenmat = Pgenmat+1i*Qgenmat;

Vg2_l= V2_Sl_calc( Sgenmat - Sload,Z,V0,'p' );
Vg_l = sqrt(Vg2_l);

Pgen = Pgenmat(Vmn);

[S0,Sg,Sl,~,Ipr,P0,Q0] = pred_S0_pscc( Pgen, Sload, Z, V0, Vp  );
%%
% Results table
TotPwr_max = max(max_pwr);
P0_l_mx = max(P0);
RR.P0 = [TotPwr_max;P0_l_mx];

Pgen_meas = Pgenmat(real(TotPwr)==TotPwr_max);
Pgen_l_mx = real(Sg(P0_l_mx==P0) + Sload);
RR.Pgen = [Pgen_meas;Pgen_l_mx];

Qgen_meas = Qgenmat(real(TotPwr)==TotPwr_max);
Qgen_l_mx = imag(Sg(P0_l_mx==P0) + Sload);
RR.Qgen = [Qgen_meas;Qgen_l_mx];

Q0_meas = imag(TotPwr(real(TotPwr)==TotPwr_max));
Q0_l_mx = Q0(P0_l_mx==P0);
RR.Q0 = [Q0_meas;Q0_l_mx];

Imax_mes = ImaxMat(real(TotPwr)==TotPwr_max);
Imax_prd = Ipr(P0_l_mx==P0)*ibN;
RR.Imax = [Imax_mes;Imax_prd];

RR.Sgen = RR.Pgen + 1i*RR.Qgen;
RR.Sg = Sg(P0_l_mx==P0);
RR.S0 = RR.P0 + 1i*RR.Q0;

RR.FF = FF;

RR.T_ratio = T_ratio;
RR.sbase = sb;
RR.Z = Z;
RR.V0 = V0;
RR.Vp = Vp;
RR.Ssc = Ssc; %kW
RR.Sload = Sload; %in pu

Vl_kst = real(Z*Sg(P0_l_mx==P0));
% RR.dV_ldc = [Vg_max - V0; Vn_ldc - 1; Vl_ldc - 1; Vw_ldc - 1];
RR.dV_kst = [Vp - V0; Vl_kst - 1];
%% pl_options - to look through the data and compare methods etc.
figs = [];
% possible: check v_comparison pgp0 pgqgq0 plosses qlosses p_d_losses
%           q_d_losses pgp0s0 vmaxmin v_comparison_lin

if sum(ismember(pl_options,'check')) || sum(ismember(pl_options,'0'))
    figs = [figs;figure('Color','White')];
    % Seems to work! :DDDDD
    cc=contourf(Pgenmat,Qgenmat,real(TotPwr)); hold on;
    clabel(cc); grid on;
    contour(Pgenmat,Qgenmat,Vgen*0.01,0.01*(1.0:0.06:1.06),'LineColor','r'); hold on;
    axis equal;
end
if sum(ismember(pl_options,'v_comparison')) || sum(ismember(pl_options,'0'))
    figs = [figs;figure('Color','White')];
    subplot(221)
    [cc,~]=contourf(Pgenmat,Qgenmat,Vgen);
    clabel(cc); axis equal;
    subplot(222)
    [cc,~]=contourf(Pgenmat,Qgenmat,Vg_l);
    clabel(cc); axis equal;
    
    subplot(224)
    [cc,~]=contourf(Pgenmat,Qgenmat,100*(Vg_l-Vgen));
    clabel(cc); %axis equal;
end
if sum(ismember(pl_options,'pgp0'))  || sum(ismember(pl_options,'0'))
    figs = [figs;figure('Color','White','Position',fig_nompos)];%,'defaultaxesfontsize',12,'defaulttextinterpreter','latex',fig_nompos);
    
    set(figs,'defaulttextinterpreter','latex');
    set(figs,'defaultaxesfontsize',12);

    plot(Pgenmat(Vmn),max_pwr,'x-'); grid on; hold on;
    plot(Pgen,P0);

    xlabel('$P_{g}$ (pu)');
    ylabel('$P_{0}$ (pu)');
%     lgnd=legend('Measured','No Load','Full Load','Weighted Load','Location','SouthEast');
%     set(lgnd,'interpreter','latex');
    axis equal; %axis([-0.002 0.050 -0.020 0.007]);
end
if sum(ismember(pl_options,'pgqgq0')) || sum(ismember(pl_options,'0'))
    figs = [figs;figure('Color','White')];
    plot(Pgenmat(Vmn),Qgenmat(Vmn),'k'); hold on; grid on;
    plot(Pgen,imag(Sg + Sload),'k:');
    plot(Pgenmat(Vmn),imag(TotPwr(Vmn)),'r'); hold on;
    plot(Pgen,Q0,'r:');

%     legend('Qg meas','Qgw pred','Qgl pred','Q0 meas','Q0w pred','Q0l pred');
    xlabel('Pg (p.u)'); ylabel('Q (p.u.)');
    axis equal;
end
if sum(ismember(pl_options,'plosses')) || sum(ismember(pl_options,'0'))
    figs = [figs;figure('Color','White','Position',[300 200 600 600])];
    plot(Pgenmat(Vmn),real(TotLss(Vmn)),'k','LineWidth',2); hold on;
    plot(Pgen,real(Sl));
%     legend('Measured losses','N losses','L losses','W losses','Location','NorthWest');
    xlabel('Pg (pu)');
    ylabel('Pl (pu)');
    axis equal; grid on;
end
if sum(ismember(pl_options,'qlosses')) || sum(ismember(pl_options,'0'))
    figs = [figs;figure('Color','White','Position',[300 200 600 600])];
    plot(Pgenmat(Vmn),imag(TotLss(Vmn)),'k','LineWidth',2); hold on;
    plot(Pgen,imag(Sl));
%     legend('N losses','L losses','W losses','Location','NorthWest');
    xlabel('Pg (pu)');
    ylabel('Ql (pu)');
    axis equal; grid on;
end
if sum(ismember(pl_options,'p_d_losses')) || sum(ismember(pl_options,'0'))
    figs = [figs;figure('Color','White','Position',[300 200 600 600])];
    TotLssMes = TotLss(Vmn);
    plot(Pgen,real(TotLssMes-Sl));
%     legend('N losses','L losses','W losses','Location','NorthWest');
    xlabel('Pg (pu)');
    ylabel('Plmes - Plpred (pu)');
    axis equal; grid on;
end
if sum(ismember(pl_options,'q_d_losses')) || sum(ismember(pl_options,'0'))
    figs = [figs;figure('Color','White','Position',[300 200 600 600])];
    TotLssMes = TotLss(Vmn);
    plot(Pgen,imag(TotLssMes-Sl));
%     legend('N losses','L losses','W losses','Location','NorthWest');
    xlabel('Pg (pu)');
    ylabel('Qlmes - Qlpred (pu)');
    axis equal; grid on;
end
if sum(ismember(pl_options,'pgp0s0')) || sum(ismember(pl_options,'0'))
    figs = [figs;figure('Color','White','Position',fig_nompos+[0 0 -100 100])];
    set(figs,'defaulttextinterpreter','latex');
    set(figs,'defaultaxesfontsize',12);

    p1 = plot(Pgenmat(Vmn),(max_pwr+real(Sload)),'x-'); grid on;  hold on;
    p2 = plot(Pgen,(P0+real(Sload))); plot(0,0);plot(0,0);
    p3 = plot(Pgenmat(Vmn),abs(max_pwr+imag(TotPwr(Vmn))),'x-');
    p4 = plot(Pgen,abs(S0));
    p5 = plot([-0.7 9.4]/(sb*1e-3),[2.5 2.5]/(sb*1e-3),'k--');
    
    axis equal; axis([-1.5 9.9 -1 14]/(sb*1e-3));
    xlabel('$P_{g}$, p.u.'); ylabel('$S_{(\cdot)}$, p.u.'); 
%     lgnd = legend([p1 p2 p3 p4 p5],'$P_{g}-P_{l}$, measrd.', '$P_{g}-P_{l}$, wghtd.',...
%             '$|S_{0}|$, measrd.','$|S_{0}|$, wghtd.','$|S_{0}|$ Sbstn. limit','Location','NorthWest');
%     set(lgnd,'Interpreter','Latex');
end
if sum(ismember(pl_options,'vmaxmin')) || sum(ismember(pl_options,'0'))
    figs = [figs;figure('Color','White','Position',fig_nompos)];%,'defaultaxesfontsize',12,'defaulttextinterpreter','latex',fig_nompos);
    plot(Pgenmat(Vmn),VmaxMat(Vmn)); grid on; hold on;
    plot(Pgenmat(Vmn),Vgen(Vmn));
    plot(Pgenmat(Vmn),VminMat(Vmn));
    plot([min(Pgenmat(Vmn)) max(Pgenmat(Vmn))],[1.06 1.06],'k--');
    plot([min(Pgenmat(Vmn)) max(Pgenmat(Vmn))],[0.94 0.94],'k--');
    legend('Max. Circuit Voltage','|Vg|','Min. Circuit Voltage',...
                    'MV Voltage Limits','Location','SouthWest');
    xlabel('Pg (pu)');
    ylabel('|V| (pu)');
end
if sum(ismember(pl_options,'imax')) || sum(ismember(pl_options,'0'))
    figs = [figs;figure('Color','White','Position',fig_nompos)];%,'defaultaxesfontsize',12,'defaulttextinterpreter','latex',fig_nompos);
    plot(Pgenmat(Vmn),ImaxMat(Vmn)); grid on; hold on;
    plot(Pgenmat(Vmn),ones(nnz(Vmn),1)*FF.Ip);
    legend('Max. Circuit Current','I+','Location','SouthWest');
    xlabel('Pgen (pu)');
    ylabel('|I| (pu)');
end









end

