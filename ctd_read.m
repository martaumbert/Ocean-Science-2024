% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load CTD data
% Figure vertical profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('/Users/MUmbert/Downloads/LSSL_ctd2018/dcruise_201881.mat') % '27-Sep-2018'

data = dcruise.data;

z = data.depth;

sal = data.salinity;
lat = dcruise.latitude;
lon = dcruise.longitude;

% TOPAZ 2018/09/15
t_sal = ncread('20180915_mm-12km-NERSC-MODEL-TOPAZ4B-ARC-RAN.fv2.0_regrided_025.nc','so');
t_lat = ncread('20180915_mm-12km-NERSC-MODEL-TOPAZ4B-ARC-RAN.fv2.0_regrided_025.nc','lat');
t_lon = ncread('20180915_mm-12km-NERSC-MODEL-TOPAZ4B-ARC-RAN.fv2.0_regrided_025.nc','lon');
aa = find(t_lat>=lat & t_lat<lat+0.25);
bb = find(t_lon>=lon & t_lon<lon+0.25);
topaz_sal = t_sal(bb,aa,:);
%% Plot
figure;set(gcf,'position',[400 1200 800 800],'color','w');
set(gca,'fontsize',32,'boxStyle','full')
plot (sal,z,'LineWidth',3);
ax(1) = gca; set(ax(1),'FontSize',24,'LineWidth',4)
xlabel('Salinity [psu] ','Fontsize',28);
ylabel('depth [m]','Fontsize',29);
set(gca, 'YDir', 'reverse')
set(gca, 'XAxisLocation', 'top')
legend ('in-situ CTD','TOPAZ4','SMOS + TOPAZ4')