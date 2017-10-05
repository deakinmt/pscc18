function [ RR,figs ] = run_pscc_feeder( FF,pl_options )

Vp = FF.Vp; Ip = FF.Ip; Pp = FF.Pp; 
fig_nompos = FF.fig_nompos;

%% Find nominal network parameters from opendss

[~, DSSObj, DSSText] = DSSStartup;
DSSCircuit=DSSObj.ActiveCircuit; DSSSolution=DSSCircuit.Solution;
DSSText.command=['Compile (',pwd,FF.filename,'.dss)'];
DSSSolution.Solve;

Zbus = create_zbus(DSSCircuit);
V0 = DSSCircuit.Vsource.pu;

DSSCircuit.SetActiveBus(FF.NUT);
vbN = sqrt(3)*DSSCircuit.ActiveBus.kVBase; % LL (@ NUT)
DSSCircuit.SetActiveBus(FF.SRC);
vbS = sqrt(3)*DSSCircuit.ActiveBus.kVbase; % LL (@ SRC)
nTr0 = vbN/vbS;
DSSCircuit.Transformers.First;
sb = DSSCircuit.Transformers.kVA/10; %ieee34mod1 script is factor of 10 out
zbN = 1e3*(vbN^2)/sb; %ohms
ibN = sb*sqrt(1/3)/vbN; %A, per phase

Z_ohm  = find_node_Z_pscc( DSSCircuit.YNodeOrder,Zbus,FF.SRC,FF.NUT,nTr0 );
Z = Z_ohm/zbN; %pu
Ssc = ((vbN^2)/abs(Z*zbN))*1e3; %kW
Sload = (1769 + 294*1i)/sb; %feeder nominal load power

% create generator at given node:
DSSText.Command=['new generator.gen phases=3 bus1=',FF.NUT,'.1.2.3 kw=0 pf=0'];
DSSText.Command='new monitor.genpq element=generator.gen terminal=1 mode=65 PPolar=no';
%% Run the DSS
[Pg, Qg] = meshgrid(Ssc*FF.pg_ssc,Ssc*FF.qg_ssc);
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
P0_mxPrm = max(P0);
RR.Psub_prm = [TotPwr_maxV;P0_mxPrm];

Pgen_prm_meas = Pgenmat(real(TotPwr)==TotPwr_maxV);
Pgen_mxV = real(Sg(P0_mxPrm==P0) + Sload);
RR.Pgen_prm = [Pgen_prm_meas;Pgen_mxV];

TotPwr_maxI = Imax_pwr(end);
Psub_hat_est = P0(find(Iest_in,1,'last'));
RR.Psub_hat = [TotPwr_maxI;Psub_hat_est];

Pgen_hat_meas = Pgenmat(real(TotPwr)==TotPwr_maxI);
Pghat_est = lemma_1( Vp, Ip/ibN, V0, Z );
Pgen_hat_est = Pghat_est + real(Sload);
% Pgen_mxI = real(Sg(P0_mxI==P0) + Sload);
RR.Pgen_hat = [Pgen_hat_meas;Pgen_hat_est];


Imax_mes = ImaxMat(real(TotPwr)==TotPwr_maxV);
Imax_prd = Iest(P0_mxPrm==P0)*ibN;
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
DAFS = 14; %default axis font size

if sum(ismember(pl_options,'pgp0'))  || sum(ismember(pl_options,'0'))
    figs = [figs;figure('Color','White','Position',fig_nompos,'defaultaxesfontsize',DAFS,'defaulttextinterpreter','latex')];
    plot(Pgenmat(Vmn),Vmax_pwr,'x-'); grid on; hold on;
    plot(PgenV,P0);
    axis equal;
    axis([0 4 -1 1.5]);
    xs = axis;
    plot([1 1]*Pgen_prm_meas,xs(3:4),'k:');

    xlabel('$P_{gen}$ (pu)');
    ylabel('$P_{0}^{Sub}$ (pu)');
    lgnd=legend('Msrd.','Estd., (5)','$P_{g}''$, msrd.','Location','NorthWest');
    set(lgnd,'interpreter','latex');

end
if sum(ismember(pl_options,'pgqgq0')) || sum(ismember(pl_options,'0'))
    figs = [figs;figure('Color','White','Position',fig_nompos,'defaultaxesfontsize',DAFS,'defaulttextinterpreter','latex')];
    plot(Pgenmat(Vmn),Qgenmat(Vmn),'x-','Color',[0.3 0.3 0.3]); hold on; grid on;
    plot(PgenV,imag(Sg + Sload),'k');
    plot(Pgenmat(Vmn),imag(TotPwr(Vmn)),'x-','Color',[1.0 0.5 0.5]); hold on;
    plot(PgenV,Q0,'r');

    lgnd = legend('$Q_{gen}$ msrd.','$Q_{gen}$ estd.','$Q_{comp}$ msrd.','$Q_{comp}$ estd.');
    set(lgnd,'Interpreter','Latex');
    xlabel('$P_{gen}$ (pu)'); ylabel('$Q$ (pu)'); axis equal;
end
if sum(ismember(pl_options,'imax')) || sum(ismember(pl_options,'0'))
    figs = [figs;figure('Color','White','Position',fig_nompos,'defaultaxesfontsize',DAFS,'defaulttextinterpreter','latex')];%
    plot(Pgenmat(Vmn),ImaxMat(Vmn),'x-'); grid on; hold on;
    plot(PgenV,Iest); grid on; hold on;
    
    xs = axis;
    plot(xs(1:2),[1 1]*FF.Ip,'k--');
    plot([1 1]*Pgen_prm_meas,xs(3:4),'k:');
    plot([1 1]*Pgen_hat_meas,xs(3:4),'k-.');
    lgnd = legend('Max $I$, msrd.','$I_{g}$, estd.','$I_{+}$','$P_{g}''$, msrd.','$\hat{P}_{g}''$, msrd.','Location','NorthWest');
    set(lgnd,'Interpreter','Latex');
    xlabel('$P_{gen}$ (pu)'); ylabel('$|I|$ (A)');
end








end

