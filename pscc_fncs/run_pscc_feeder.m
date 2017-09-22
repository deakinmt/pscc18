function [ RR,figs ] = run_pscc_feeder( FF,pl_options )

% INPUTS TO DETERMINE BEHAVIOUR -------------------
Vp = FF.Vp; % we set for this network.
Ip = FF.Ip;
Pp = FF.Pp;
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

V0 = DSSCircuit.Vsource.pu;

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

% DSSObj.Text.Command='export elempowers';
% CapQ = calc_capQ(elmpwr_fnm)/sb; %pu 
% LdNom = (calc_load(elmpwr_fnm)/sb); %pu
% Sload = CapQ + LdNom; %pu

Sload = (1769 + 294*1i)/sb;

% create generator at given node:
DSSText.Command=['new generator.gen phases=3 bus1=',NUT,'.1.2.3 kw=0 pf=0'];
DSSText.Command='new monitor.genpq element=generator.gen terminal=1 mode=65 PPolar=no';

%% Run the DSS
[Pg, Qg] = meshgrid(Ssc*pg_ssc,Ssc*qg_ssc);
Pf = sign(Qg).*Pg./sqrt(Pg.^2 + Qg.^2);

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

ImaxMat = reshape(ImaxMes,size(Pg));
VmaxMat = reshape(VmaxMes,size(Pg));
VminMat = reshape(VminMes,size(Pg));

DSSMon=DSSCircuit.Monitors; DSSMon.name='genpq';
DSSText.Command='show monitor genpq';
Pgen_lin = -ExtractMonitorData(DSSMon,1,1)/sb; %pu, +ve is generation
Qgen_lin = -ExtractMonitorData(DSSMon,2,1)/sb; %pu, +ve is generation

TotLss = reshape(totlss(:,1),size(Pg)) + 1i*reshape(totlss(:,2),size(Pg));% +ve implies losses
TotPwr = reshape(totpwr(:,1),size(Pg)) + 1i*reshape(totpwr(:,2),size(Pg)); %-ve implies load

Pgenmat=reshape(Pgen_lin,size(Pg));% - real(TotLd);
Qgenmat=reshape(Qgen_lin,size(Pg));% - imag(TotLd);
Sgenmat = Pgenmat+1i*Qgenmat;

Vg2= V2_Sl_calc(Sgenmat - Sload,Z,V0,'p' );
Vg = sqrt(Vg2);
%%
% calculate maximum power export s.t. voltage constraint:
VNaN_outs = 0./(VmaxMat<Vp);
Vg_inV = VmaxMat + VNaN_outs; % remove numbers where the voltage is too high
Vmn = (Vg_inV==max(Vg_inV));
Vmax_pwr = max(real(TotPwr(Vmn)),[],2); % for each row, find the max power
PgenV = Pgenmat(Vmn);

% Vmax_pwr = max(real(TotPwr(Vg_in_mnV)),[],2); % for each row, find the max power
% Vmn = find(real(TotPwr)==repmat(Vmax_pwr',numel(TotPwr)/numel(Vmax_pwr),1)); %find the index of these points
% PgenV = Pgenmat(Vmn);

% calculate maximum power export s.t. current:
INaN_outs = 0./((ImaxMat<Ip).*(real(TotPwr)<Pp).*(VmaxMat<Vp).*(VmaxMat>1.05)) ;
Vg_inI = VmaxMat + INaN_outs + VNaN_outs; % remove high voltages and currents

Imn = (Vg_inI==max(Vg_inI));
Imax_pwr = max(real(TotPwr(Imn)),[],2); % for each row, find the max power
% PgenI = Pgenmat(Imn);


[S0,Sg,Sl,~,Iest_pu,P0,Q0] = pred_S0_pscc( PgenV, Sload, Z, V0, Vp  );
Iest = Iest_pu*ibN;
Iest_in = ((Iest<Ip).*(P0<Pp));
%%
% Results table
TotPwr_maxV = max(Vmax_pwr);
P0_mxV = max(P0);
RR.P0V = [TotPwr_maxV;P0_mxV];
Pgen_measV = Pgenmat(real(TotPwr)==TotPwr_maxV);
Pgen_mxV = real(Sg(P0_mxV==P0) + Sload);
RR.PgenV = [Pgen_measV;Pgen_mxV];
Qgen_measV = Qgenmat(real(TotPwr)==TotPwr_maxV);
Qgen_mxV = imag(Sg(P0_mxV==P0) + Sload);
RR.QgenV = [Qgen_measV;Qgen_mxV];
Q0_measV = imag(TotPwr(real(TotPwr)==TotPwr_maxV));
Q0_mxV = Q0(P0_mxV==P0);
RR.Q0V = [Q0_measV;Q0_mxV];

TotPwr_maxI = Imax_pwr(end);
% P0_I = P0(Iest_in); %remove overcurrents
P0_mxI = P0(find(Iest_in,1,'last'));
RR.P0I = [TotPwr_maxI;P0_mxI];
Pgen_measI = Pgenmat(real(TotPwr)==TotPwr_maxI);
Pgen_mxI = real(Sg(P0_mxI==P0) + Sload);
RR.PgenI = [Pgen_measI;Pgen_mxI];


% RR.SgenV = RR.PgenV + 1i*RR.QgenV;
% RR.SgV = Sg(P0_mxV==P0);
% RR.S0V = RR.P0V + 1i*RR.Q0V;

Imax_mes = ImaxMat(real(TotPwr)==TotPwr_maxV);
Imax_prd = Iest(P0_mxV==P0)*ibN;
RR.Imax = [Imax_mes;Imax_prd];

RR.FF = FF;

RR.sbase = sb;
RR.Z = Z;
RR.Zabs = abs(Z);
RR.lz = lambdas(Z);
RR.V0 = V0;
RR.Ssc = Ssc; %kW
RR.Sload = Sload; %in pu

%% pl_options - to look through the data and compare methods etc.
figs = [];
% possible: check v_comparison pgp0 pgqgq0 plosses qlosses p_d_losses
%           q_d_losses pgp0s0 vmaxmin v_comparison_lin

if sum(ismember(pl_options,'check')) || sum(ismember(pl_options,'0'))
    figs = [figs;figure('Color','White')];
    cc=contourf(Pgenmat,Qgenmat,real(TotPwr)); hold on;
    clabel(cc); grid on;
    contour(Pgenmat,Qgenmat,VmaxMat*0.01,0.01*(1.0:0.06:1.06),'LineColor','r'); hold on;
    axis equal;
end
if sum(ismember(pl_options,'v_comparison')) || sum(ismember(pl_options,'0'))
    figs = [figs;figure('Color','White')];
    subplot(221)
    [cc,~]=contourf(Pgenmat,Qgenmat,VmaxMat);
    clabel(cc); axis equal;
    subplot(222)
    [cc,~]=contourf(Pgenmat,Qgenmat,Vg);
    clabel(cc); axis equal;
    
    subplot(224)
    [cc,~]=contourf(Pgenmat,Qgenmat,100*(Vg-VmaxMat));
    clabel(cc); %axis equal;
end
if sum(ismember(pl_options,'pgp0'))  || sum(ismember(pl_options,'0'))
    figs = [figs;figure('Color','White','Position',fig_nompos,'defaultaxesfontsize',12,'defaulttextinterpreter','latex')];
    
    set(figs,'defaulttextinterpreter','latex');
    set(figs,'defaultaxesfontsize',12);

    plot(Pgenmat(Vmn),Vmax_pwr,'x-'); grid on; hold on;
    plot(PgenV,P0);

    xlabel('$P_{g}$ (pu)');
    ylabel('$P_{0}$ (pu)');
    lgnd=legend('Measured','2 bus P-V curve','Location','SouthEast');
    set(lgnd,'interpreter','latex');
    axis equal;
end
if sum(ismember(pl_options,'pgqgq0')) || sum(ismember(pl_options,'0'))
    figs = [figs;figure('Color','White','Position',fig_nompos,'defaultaxesfontsize',12,'defaulttextinterpreter','latex')];
    plot(Pgenmat(Vmn),Qgenmat(Vmn),'x-','Color',[0.3 0.3 0.3]); hold on; grid on;
    plot(PgenV,imag(Sg + Sload),'k');
    plot(Pgenmat(Vmn),imag(TotPwr(Vmn)),'x-','Color',[1.0 0.5 0.5]); hold on;
    plot(PgenV,Q0,'r');

    lgnd = legend('$Q_{gen}$ msrd','$Q_{gen}$ estd','$Q_{comp}$ msrd','$Q_{comp}$ estd');
    set(lgnd,'Interpreter','Latex');
    xlabel('$P_{gen}$ (pu)'); ylabel('$Q$ (pu)'); axis equal;
end
if sum(ismember(pl_options,'plosses')) || sum(ismember(pl_options,'0'))
    figs = [figs;figure('Color','White','Position',[300 200 600 600])];
    plot(Pgenmat(Vmn),real(TotLss(Vmn)),'k','LineWidth',2); hold on;
    plot(PgenV,real(Sl));
    legend('Measured P losses','Pred P losses','Location','NorthWest');
    xlabel('Pg (pu)');
    ylabel('Pl (pu)');
    axis equal; grid on;
end
if sum(ismember(pl_options,'qlosses')) || sum(ismember(pl_options,'0'))
    figs = [figs;figure('Color','White','Position',[300 200 600 600])];
    plot(Pgenmat(Vmn),imag(TotLss(Vmn)),'k','LineWidth',2); hold on;
    plot(PgenV,imag(Sl));
    legend('Meas Q losses','Pred Q losses','Location','NorthWest');
    xlabel('Pg (pu)');
    ylabel('Ql (pu)');
    axis equal; grid on;
end
if sum(ismember(pl_options,'pgp0s0')) || sum(ismember(pl_options,'0'))
    figs = [figs;figure('Color','White','Position',fig_nompos+[0 0 -100 100])];
    set(figs,'defaulttextinterpreter','latex');
    set(figs,'defaultaxesfontsize',12);

    p1 = plot(Pgenmat(Vmn),(Vmax_pwr+real(Sload)),'x-'); grid on;  hold on;
    p2 = plot(PgenV,(P0+real(Sload))); plot(0,0);plot(0,0);
    p3 = plot(Pgenmat(Vmn),abs(Vmax_pwr+imag(TotPwr(Vmn))),'x-');
    p4 = plot(PgenV,abs(S0));
    p5 = plot([-0.7 9.4]/(sb*1e-3),[2.5 2.5]/(sb*1e-3),'k--');
    
    axis equal; axis([-1.5 9.9 -1 14]/(sb*1e-3));
    xlabel('$P_{g}$, p.u.'); ylabel('$S_{(\cdot)}$, p.u.'); 
    lgnd = legend([p1 p2 p3 p4 p5],'$P_{g}-P_{l}$, measrd.', '$P_{g}-P_{l}$, wghtd.',...
            '$|S_{0}|$, measrd.','$|S_{0}|$, wghtd.','$|S_{0}|$ Sbstn. limit','Location','NorthWest');
    set(lgnd,'Interpreter','Latex');
end
if sum(ismember(pl_options,'vmaxmin')) || sum(ismember(pl_options,'0'))
    figs = [figs;figure('Color','White','Position',fig_nompos)];%,'defaultaxesfontsize',12,'defaulttextinterpreter','latex',fig_nompos);
    plot(Pgenmat(Vmn),VmaxMat(Vmn)); grid on; hold on;
    plot(Pgenmat(Vmn),VmaxMat(Vmn));
    plot(Pgenmat(Vmn),VminMat(Vmn));
    plot([min(Pgenmat(Vmn)) max(Pgenmat(Vmn))],[1 1]*Vp,'k--');
    plot([min(Pgenmat(Vmn)) max(Pgenmat(Vmn))],[0.94 0.94],'k--');
    legend('Max. Circuit Voltage','|Vg|','Min. Circuit Voltage',...
                    'MV Voltage Limits','Location','SouthWest');
    xlabel('Pg (pu)');
    ylabel('|V| (pu)');
end
if sum(ismember(pl_options,'imax')) || sum(ismember(pl_options,'0'))
    figs = [figs;figure('Color','White','Position',fig_nompos,'defaultaxesfontsize',12,'defaulttextinterpreter','latex')];%
    plot(Pgenmat(Vmn),ImaxMat(Vmn),'x-'); grid on; hold on;
    plot(PgenV,Iest); grid on; hold on;
    xs = axis;
    plot(xs(1:2),[1 1]*FF.Ip,'k--');
    plot([1 1]*Pgen_measV,xs(3:4),'k:');
    lgnd = legend('Max $I$, meas.','$I_{g}$, estimated','$I_{+}$','$P_{g}''$, meas','Location','NorthWest');
    set(lgnd,'Interpreter','Latex');
    xlabel('$P_{gen}$ (pu)'); ylabel('$|I|$ (A)');
end









end

