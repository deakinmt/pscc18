

% for I ratings, see Kersting's 2001 'Radial Distribution Test Feeders'
% paper.

fig_loc = 'C:\Users\chri3793\Documents\DPhil\malcolm_updates\wc170918\figures';

% line types:
v_max = 1.06; %pu;
n_ph = 3; % no. phases

types = {'aa1000','acsr5565','aa500','acsr3364','aa250','acrs4/0','aa2/0','acrs1/0',...
            'aa1/0','aa2','acsr2','acsr4','acsr4','cu10','cu12','cu14'};
i_ratings = [698; 730; 783; 530; 329; 340; 230; 230; 310; 156; 180; 140; 80; 75; 20];

% I ratings of the IEEE 34 bus feeder
% 34bus_linetypes = {'acsr1/0','acsr2'}; % double check for 301 (the second one)?????

ohl_i_rating34 = i_ratings([8,11,12]);
sb34 = 2500;
vb34 = 24.9/sqrt(3);

s_rating34 = v_max*ohl_i_rating34*vb34*n_ph;

fig = figure('Color','White');

figname = [fig_loc,'\psjul_figs_line_ratings'];
plot(s_rating34/sb34,'x-');
xlabel('Conductor type'); ylabel('OHL rating (Sline/Sbase)');
legend('IEEE13','IEEE34');
title('Thermal line ratings of OHLs c.f. transformers');

% export_fig(fig,figname);
% export_fig(fig,[figname,'.pdf'],'-dpdf');