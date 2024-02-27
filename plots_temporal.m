%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% 1. Read fwc outputs
% 2. Read original monthly smos
% 3. Read in-situ estimates
% 4. Create temporal evolution plot
% -------------------------------------------------------------------------
% Author    : M. Umbert 
% Email     : mumbert@icm.csic.es 
% Creation  : 20/02/2023
% Modified  : 25/04/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all
%%
% area of each pixel
load('area.mat');
area = permute(area,[2 1]);

% Beaufort Gyre limits (bx by)
load poligon_bg.mat

%% 1. Read FWC outputs
% READ The Original NETCDF File and generate a .mat
addpath('/Users/mumbert/Library/Mobile Documents/com~apple~CloudDocs/Analysis/Curie/fwc/outputs_smos_topaz/');

filename = 'topaz4b_monthly_fullArctic_fwc2011_2019_latlon_vSmean.nc';
filename16 = 'smos_topaz4b_monthly_fullArctic_fwc2011_2019_16m_vSmean.nc';
filename22 = 'smos_topaz4b_monthly_fullArctic_fwc2011_2019_22m_vSmean.nc';
filename25 = 'smos_topaz4b_monthly_fullArctic_fwc2011_2019_25m_vSmean.nc';
filename29 = 'smos_topaz4b_monthly_fullArctic_fwc2011_2019_29m_vSmean.nc';

lat       = ncread(filename,'latitude');
lon       = ncread(filename,'longitude');
time       = ncread(filename,'time');
sss      = ncread(filename,'sss');
sst     = ncread(filename,'sst');
fwc    = ncread(filename,'fwc');
fwc_16 = ncread(filename16,'fwc');
fwc_22 = ncread(filename22,'fwc');
fwc_25 = ncread(filename25,'fwc');
fwc_29 = ncread(filename29,'fwc');


[yy2,xx2] = meshgrid(lat,lon);
X = time/(60*60*24)+datenum(1970,1,1);
X2 = datestr(X);X3 = datetime(X2);

%% 2. OPEN FWC from INSITU DATA
% All data are in meters. Uncertainties for each grid cell are determined using 
% the optimal interpolation technique previously described in Proshutinsky et al. (2009)

% Execute /fwc/scripts/open_bgfwconint.m
% Execute /fwc/scripts/open_bgfwconint_error.m

% Mean temporal evolution over the region
FWC = [fwc_2003,fwc_2004,fwc_2005,fwc_2006,fwc_2007,fwc_2008,fwc_2009,fwc_2010,fwc_2011,fwc_2012,fwc_2013,fwc_2014,fwc_2015,fwc_2016,fwc_2017,fwc_2018,fwc_2019,fwc_2020];

years = 2003:1:2020;
for i = 1:length(years)
    mean_fwc(i) = nanmean(FWC(:,i));   
end

FWC_err = [err_2003,err_2004,err_2005,err_2006,err_2007,err_2008,err_2009,err_2010,err_2011,err_2012,err_2013,err_2014,err_2015,err_2016,err_2017,err_2018];

years_e = 2003:1:2018;
for i = 1:length(years_e)
    mean_err(i) = nanmean(FWC_err(:,i));   
end

%% Fill a vector mean_fwc_ex and mean_err_ex of same time lenght as monthly smos time period with in-situ
x = 1:108; x = x';
for i = 1:108
    if (i == 1)
        mean_fwc_ex(i) = mean_fwc(9);
        mean_err_ex(i) = mean_err(9);
    elseif (i == 12*1)
        mean_fwc_ex(i) = mean_fwc(10);
        mean_err_ex(i) = mean_err(10);
    elseif (i == 12*2)
        mean_fwc_ex(i) = mean_fwc(11);
        mean_err_ex(i) = mean_err(11);
    elseif (i == 12*3)
        mean_fwc_ex(i) = mean_fwc(12);
        mean_err_ex(i) = mean_err(12);
    elseif (i == 12*4)
        mean_fwc_ex(i) = mean_fwc(13);
        mean_err_ex(i) = mean_err(13);
    elseif (i == 12*5)
        mean_fwc_ex(i) = mean_fwc(14);
        mean_err_ex(i) = mean_err(14);
    elseif (i == 12*6)
        mean_fwc_ex(i) = mean_fwc(15);
        mean_err_ex(i) = mean_err(15);
    elseif (i == 12*7)
        mean_fwc_ex(i) = mean_fwc(16);
        mean_err_ex(i) = mean_err(16);
    elseif (i == 12*8)
        mean_fwc_ex(i) = mean_fwc(17);
        mean_err_ex(i) = mean_err(16);
    elseif (i == 12*9)
        mean_fwc_ex(i) = mean_fwc(18);
        mean_err_ex(i) = mean_err(16);
    else
        mean_fwc_ex(i) = nan;
        mean_err_ex(i) = nan;
    end
end
figure
set(gcf,'position',[249 608 1200 600],'color','w');
errorbar(x,mean_fwc_ex,mean_err_ex,'Linewidth',2,'Marker','diamond','color',rgb('pink'))
set(gca,'fontsize',26)
grid on; 
ylabel('fwc [m]')


%% Load SMOS monthly data

% 3. Load SMOS SSS (monthly data) ----------------------------------------- 

smosfile = 'sss_smos_bec_2011_2019_regrided_025_monthly_data_40N.nc';

slat     = ncread(smosfile,'lat');
slon     = ncread(smosfile,'lon');
stime    = ncread(smosfile,'time');
ssss     = ncread(smosfile,'sss');


%% Poligon shape
% pgon = polyshape(bx,by); % bx and by have been loaded in 'poligon_bg.mat'
% 
% TFin=zeros(size(xx2));
% TFon=zeros(size(xx2));
% 
% for i=1:length(xx2(:,1))
%     [TFin(i,:), TFon(i,:)] = isinterior(pgon,xx2(i,:),yy2(i,:));%en TFin guardamos los puntos dentro del pol??gono (incluyendo la frontera), en TFon identificamos s??lo los puntos en la frontera
% end

%%
load('pgon.mat')
load('TFon.mat')
load('TFin.mat')

% BG Mask for data selection
mask_bg = TFin;
mask_bg(mask_bg ==0) = nan;
lon_bg=xx2.*mask_bg;
lat_bg=yy2.*mask_bg;

% topaz bg selection

for i=1:108
    
    BG_topaz_fwc(:,:,i) = (fwc(:,:,i)).*mask_bg;
    BG_topaz_fwc16(:,:,i) = (fwc_16(:,:,i)).*mask_bg;
    BG_topaz_fwc22(:,:,i) = (fwc_22(:,:,i)).*mask_bg;
    BG_topaz_fwc25(:,:,i) = (fwc_25(:,:,i)).*mask_bg;
    BG_topaz_fwc29(:,:,i) = (fwc_29(:,:,i)).*mask_bg;
    
    BG_topaz_sss(:,:,i) = (sss(:,:,i)).*mask_bg;
    
    % smos bg selection
    BG_smos_sss(:,:,i) = (ssss(:,:,i)).*mask_bg;    
end

% % %PLOT: time evolution of BG sea ice extension, fwc and sss 
mean_mens_fwc = squeeze(nanmean(nanmean(BG_topaz_fwc)));
mean_mens_fwc16 = squeeze(nanmean(nanmean(BG_topaz_fwc16)));
mean_mens_fwc22 = squeeze(nanmean(nanmean(BG_topaz_fwc22)));
mean_mens_fwc25 = squeeze(nanmean(nanmean(BG_topaz_fwc25)));
mean_mens_fwc29 = squeeze(nanmean(nanmean(BG_topaz_fwc29)));

%% PLOT  
figure
plot(x,mean_mens_fwc16,'color',rgb('slategrey'),'LineWidth',2)
hold on
%plot(x,mean_mens_fwc22,'color',rgb('green'),'LineWidth',2)
plot(x,mean_mens_fwc25,'color',rgb('forestgreen'),'LineWidth',2)
plot(x,mean_mens_fwc29,'color',rgb('gold'),'LineWidth',2)
plot(x,mean_mens_fwc,'color',rgb('black'),'LineWidth',2)
errorbar(x,mean_fwc_ex,mean_err_ex,'Linewidth',2,'Marker','diamond','color',rgb('darkviolet'))
ylabel('fwc (m)');
xlabel(' '); xticks(1:12:109)
xticklabels({'2011', '2012','2013','2014','2015','2016','2017','2018','2019','2020'})
set(gcf,'color','w','position',[249 608 1200 500]);
set(gca,'Fontsize',24)
grid on
legend('SMOS 16m','SMOS 25m','SMOS 29m','TOPAZ','In situ','Fontsize',18,'Location','best')


%% Compute Percentage increase = [(new value - initial value) / initial value] x 100%

Pinc = ((mean_mens_fwc29 - mean_mens_fwc ) / mean_mens_fwc) * 100;
figure
plot(x,Pinc,'color',rgb('slategrey'),'LineWidth',2)


%% SMOS-Mask for data selection 
mask_smos = double(~isnan(BG_smos_sss));
mask_smos(mask_smos == 0) = NaN;

BG_topaz_sss_smosmask = BG_topaz_sss.*mask_smos;
BG_topaz_fwc_smosmask = BG_topaz_fwc.*mask_smos;
BG_topaz16_fwc_smosmask = BG_topaz_fwc16.*mask_smos;
BG_topaz22_fwc_smosmask = BG_topaz_fwc22 .*mask_smos;
BG_topaz25_fwc_smosmask = BG_topaz_fwc25 .*mask_smos;
BG_topaz29_fwc_smosmask = BG_topaz_fwc29 .*mask_smos;

% % %PLOT: time evolution of BG sea ice extension, fwc and sss 
mean_mens_fwc_smosmask = squeeze(nanmean(nanmean(BG_topaz_fwc_smosmask)));
mean_mens_fwc16_smosmask = squeeze(nanmean(nanmean(BG_topaz16_fwc_smosmask)));
mean_mens_fwc22_smosmask = squeeze(nanmean(nanmean(BG_topaz22_fwc_smosmask)));
mean_mens_fwc25_smosmask = squeeze(nanmean(nanmean(BG_topaz25_fwc_smosmask)));
mean_mens_fwc29_smosmask = squeeze(nanmean(nanmean(BG_topaz29_fwc_smosmask)));

% PLOT 

i = 1;
for tt = 2011:2019
    kk = find(year(X3)==tt & (month(X3)==7 | month(X3)==8 | month(X3)==9 | month(X3)==10)); % from July to October to emulate the means from in-situ
    fwc_year_topaz(:,:,i) = nanmean(BG_topaz_fwc(:,:,kk),3);
    fwc_year_topaz16(:,:,i) = nanmean(BG_topaz_fwc16(:,:,kk),3);
    fwc_year_topaz22(:,:,i) = nanmean(BG_topaz_fwc22(:,:,kk),3);
    fwc_year_topaz25(:,:,i) = nanmean(BG_topaz_fwc25(:,:,kk),3);
    fwc_year_topaz29(:,:,i) = nanmean(BG_topaz_fwc29(:,:,kk),3);
    fwc_year_topaz_smosmask(:,:,i) = nanmean(BG_topaz_fwc_smosmask(:,:,kk),3);
    fwc_year_topaz16_smosmask(:,:,i) = nanmean(BG_topaz16_fwc_smosmask(:,:,kk),3);
    fwc_year_topaz22_smosmask(:,:,i) = nanmean(BG_topaz22_fwc_smosmask(:,:,kk),3);
    fwc_year_topaz25_smosmask(:,:,i) = nanmean(BG_topaz25_fwc_smosmask(:,:,kk),3);
    fwc_year_topaz29_smosmask(:,:,i) = nanmean(BG_topaz29_fwc_smosmask(:,:,kk),3);
    i = i +1;
end

% Plot yearly evolution of fwc and in-situ only where smos and mean only where smos
FWC_situ = [fwc_2011,fwc_2012,fwc_2013,fwc_2014,fwc_2015,fwc_2016,fwc_2017,fwc_2018,fwc_2019];
FWC_error_situ = [err_2011,err_2012,err_2013,err_2014,err_2015,err_2016,err_2017,err_2018,err_2018];

x = 1:9; x = x';
figure
errorbar(x,squeeze(nanmean(nanmean(fwc_year_topaz16))),squeeze(nanstd(nanstd(fwc_year_topaz16))),'Linewidth',1,'Marker','diamond','color',rgb('slategrey'))
hold on
%errorbar(x,squeeze(nanmean(nanmean(fwc_year_topaz22))),squeeze(nanstd(nanstd(fwc_year_topaz22))),'Linewidth',1,'Marker','diamond')
errorbar(x,squeeze(nanmean(nanmean(fwc_year_topaz25))),squeeze(nanstd(nanstd(fwc_year_topaz25))),'Linewidth',1,'Marker','diamond','color',rgb('forestgreen'))
errorbar(x,squeeze(nanmean(nanmean(fwc_year_topaz29))),squeeze(nanstd(nanstd(fwc_year_topaz29))),'Linewidth',1,'Marker','diamond','color',rgb('gold'))
errorbar(x,squeeze(nanmean(nanmean(fwc_year_topaz))),squeeze(nanstd(nanstd(fwc_year_topaz))),'Linewidth',1,'Marker','diamond','color',rgb('black'))
errorbar(x,nanmean(FWC_situ),nanmean(FWC_error_situ),'Linewidth',2,'Marker','diamond','color',rgb('darkviolet'))
ylabel('fwc (m)');
%xlabel('time'); 
xticks(1:1:9)
xticklabels({'2011', '2012','2013','2014','2015','2016','2017','2018','2019','2020'})
set(gcf,'color','w','position',[249 608 1200 500]);
set(gca,'Fontsize',22)
grid on
legend('SMOS 16m','SMOS 25m','SMOS 29m','TOPAZ','In situ','Fontsize',18,'Location','best')

% Plot yearly evolution of fwc and in-situ only where smos and mean only where smos

x = 1:9; x = x';
figure
errorbar(x,squeeze(nanmean(nanmean(fwc_year_topaz16_smosmask))),squeeze(nanstd(nanstd(fwc_year_topaz16_smosmask))),'Linewidth',1,'Marker','diamond','color',rgb('slategrey'))
hold on
errorbar(x,squeeze(nanmean(nanmean(fwc_year_topaz22_smosmask))),squeeze(nanstd(nanstd(fwc_year_topaz22_smosmask))),'Linewidth',1,'Marker','diamond')
errorbar(x,squeeze(nanmean(nanmean(fwc_year_topaz25_smosmask))),squeeze(nanstd(nanstd(fwc_year_topaz25_smosmask))),'Linewidth',1,'Marker','diamond','color',rgb('forestgreen'))
errorbar(x,squeeze(nanmean(nanmean(fwc_year_topaz29_smosmask))),squeeze(nanstd(nanstd(fwc_year_topaz29_smosmask))),'Linewidth',1,'Marker','diamond','color',rgb('gold'))
errorbar(x,squeeze(nanmean(nanmean(fwc_year_topaz_smosmask))),squeeze(nanstd(nanstd(fwc_year_topaz_smosmask))),'Linewidth',1,'Marker','diamond','color',rgb('black'))
errorbar(x,nanmean(FWC_situ),nanmean(FWC_error_situ),'Linewidth',2,'Marker','diamond','color',rgb('darkviolet'))
ylabel('fwc (m)');
%xlabel('time'); 
xticks(1:1:9)
xticklabels({'2011', '2012','2013','2014','2015','2016','2017','2018','2019','2020'})
set(gcf,'color','w','position',[249 608 1200 500]);
set(gca,'Fontsize',22)
grid on
legend('SMOS 16m','SMOS 22m','SMOS 25m','SMOS 29m','TOPAZ','In situ','Fontsize',18,'Location','best')


%% FWC estimates in km3 and plots to compare with Protushinsky et al

for i=1:108
    %fwc_km3(:,:,i) =(fwc(:,:,i)*0.001).*area;% fwc in km3/pixel
     fwc_km3(:,:,i) =(BG_topaz_fwc(:,:,i)*0.001).*area;
end

i = 1;
for tt = 2011:2019
    kk = find(year(X3)==tt & (month(X3)==7 | month(X3)==8 | month(X3)==9 | month(X3)==10)); % from July to October to emulate the means from in-situ
    fwc_km3_year(:,:,i) = nanmean(fwc_km3(:,:,kk),3);
    i = i +1;
end

x = 1:9; x = x';
figure
%plot(x,squeeze(nansum(nansum(fwc_km3))),'color',rgb('black'),'LineWidth',2); hold on
bar(x,squeeze(nansum(nansum(fwc_km3_year))))
hold on
%plot(x,mean_mens_fwc40,'LineWidth',2)
%errorbar(x,fwc_km3,mean_err_ex,'Linewidth',2,'Marker','diamond','color',rgb('darkviolet'))
ylabel('FWC (km^3)');
%xlabel('time'); 
xticks(1:1:9)
xticklabels({'2011', '2012','2013','2014','2015','2016','2017','2018','2019','2020'})
set(gcf,'color','w','position',[249 608 1200 500]);
set(gca,'Fontsize',22)
grid on
legend('TOPAZ','Fontsize',18,'Location','best')


% New plot: year estimates of in-situ against TOPAZ yearly estimates

x = 1:9; x = x';

coefficients = polyfit(x, nanmean(FWC_situ), 1);  % Fit a first-degree polynomial (linear fit)
slope = coefficients(1);
intercept = coefficients(2);
xfit = min(x):0.1:max(x);
yfit = polyval(coefficients, xfit);

figure
errorbar(x,nanmean(FWC_situ),nanmean(FWC_error_situ),'Linewidth',2,'Marker','diamond','color',rgb('darkviolet'));hold on
plot(xfit, yfit, '--r', 'LineWidth', 2);

coefficients = polyfit(x,squeeze(nanmean(nanmean(fwc_year_topaz))), 1);  % Fit a first-degree polynomial (linear fit)
slope = coefficients(1);
intercept = coefficients(2);
xfit = min(x):0.1:max(x);
yfit = polyval(coefficients, xfit);

errorbar(x,squeeze(nanmean(nanmean(fwc_year_topaz))),squeeze(nanstd(nanstd(fwc_year_topaz))),'Linewidth',2,'Marker','diamond','color',rgb('black'))
plot(xfit, yfit, '--k', 'LineWidth', 2);

ylabel('fwc (m)');
%xlabel('time'); 
xticks(1:1:9)
xticklabels({'2011', '2012','2013','2014','2015','2016','2017','2018','2019','2020'})
set(gcf,'color','w','position',[249 608 1200 500]);
set(gca,'Fontsize',22)
grid on
legend('In-situ','Trend = ','TOPAZ','Trend = ','Fontsize',18,'Location','best')
%saveas(gcf,'FWC_km3_temporal','png')


%% Scatter plots In-situ FWC / TOPAZ+SMOS FWC
% months 8, 9, 10 maximum coverage of smos data in the region

FWC_situ = [fwc_2011,fwc_2012,fwc_2013,fwc_2014,fwc_2015,fwc_2016,fwc_2017,fwc_2018,fwc_2019];

% Interpolate TOPAZ+SMOS estimates to in-situ yearly (jul-oct) and 50*50 km sapling
i = 1;
for tt = 2011:2019
    kk = find(year(X3)==tt & (month(X3)==8 | month(X3)==9 | month(X3)==10));
    BG_topaz_fwc_mean(:,:,i) = nanmean(BG_topaz_fwc(:,:,kk),3);
    BG_smos16_fwc_mean(:,:,i) = nanmean(BG_topaz_fwc16(:,:,kk),3);
    BG_smos25_fwc_mean(:,:,i) = nanmean(BG_topaz_fwc25(:,:,kk),3);
    BG_smos29_fwc_mean(:,:,i) = nanmean(BG_topaz_fwc29(:,:,kk),3);

    FWC_topaz_int(:,i) = interpn(xx2,yy2,squeeze(BG_topaz_fwc_mean(:,:,i)),Longitude, Latitude);
    FWC_smos16_int(:,i) = interpn(xx2,yy2,squeeze(BG_smos16_fwc_mean(:,:,i)),Longitude, Latitude);
    FWC_smos25_int(:,i) = interpn(xx2,yy2,squeeze(BG_smos25_fwc_mean(:,:,i)),Longitude, Latitude);
    FWC_smos29_int(:,i) = interpn(xx2,yy2,squeeze(BG_smos29_fwc_mean(:,:,i)),Longitude, Latitude);

    i = i +1;
end



%% Scatterplots

% TOPAZ Only
%X_vec = FWC_situ(:); Y_vec = FWC_topaz_int(:); % tots els anys
for i = 1:9
    X_vec = FWC_situ(:,i); X_vec = X_vec(:); Y_vec = FWC_topaz_int(:,i);Y_vec = Y_vec(:);
    nn = find(~isnan(X_vec) & ~isnan(Y_vec));
    X_vec = X_vec(nn);Y_vec = Y_vec(nn);

    figure;set(gcf,'position',[400 1200 800 800],'color','w');
    set(gca,'fontsize',16,'boxStyle','full')
    scatter(X_vec(:),Y_vec(:),26,[0.5 0.5 0.5],'filled');hold on;
    xlim([10 35]);ylim([10 35])
    plot([10 35],[10 35],'--r','LineWidth',3)
    hold on
    [stat.curve, stat.gof,stat.output] = fit(double(X_vec(:)), double(Y_vec(:)), 'poly1' );
    ci = predint(stat.curve,Longitude,0.99);
    p1=plot(stat.curve,'--k');set(p1,'LineWidth',3);legend off;grid on
    [R p] = corrcoef(X_vec(:),Y_vec(:));
    str1(1) = {['n            : ' num2str(stat.output.numobs,'%6.0f')]};
    str1(2) = {['Slope      : ',num2str(stat.curve.p1,'%6.2f')]};
    str1(3) = {['Y-intercept  : ',num2str(stat.curve.p2,'%6.2f')]};
    str1(4) = {['R^2        : ',num2str(stat.gof.rsquare,'%6.2f')]};
    str1(5) = {['bias        : ',num2str((nanmean(X_vec(:)-Y_vec(:))),'%6.2f')]};
    str1(6) = {['std        : ',num2str((nanstd(X_vec(:)-Y_vec(:))),'%6.2f')]};
    ht      = text(11,31,str1);
    set(ht,'fontsize',28,'FontName','times','fontweight','b')
    ax(1) = gca; set(ax(1),'FontSize',24,'LineWidth',4)
    xlabel('fwc [m] In situ ','Fontsize',28);
    ylabel('fwc [m] TOPAZ ','Fontsize',28);
end

% TOPAZ+SMOS 16m
%X_vec = FWC_situ(:); Y_vec = FWC_smos16_int(:);
for i = 1:9
    X_vec = FWC_situ(:,i); X_vec = X_vec(:); Y_vec = FWC_smos16_int(:,i);Y_vec = Y_vec(:);

    nn = find(~isnan(X_vec) & ~isnan(Y_vec));
    X_vec = X_vec(nn);Y_vec = Y_vec(nn);

    figure;set(gcf,'position',[400 1200 800 800],'color','w');
    set(gca,'fontsize',16,'boxStyle','full')
    scatter(X_vec(:),Y_vec(:),26,[0.5 0.5 0.5],'filled');hold on;
    xlim([10 35]);ylim([10 35])
    plot([10 35],[10 35],'--r','LineWidth',3)
    hold on
    [stat.curve, stat.gof,stat.output] = fit(double(X_vec(:)), double(Y_vec(:)), 'poly1' );
    ci = predint(stat.curve,Longitude,0.99);
    p1=plot(stat.curve,'--k');set(p1,'LineWidth',3);legend off;grid on
    [R p] = corrcoef(X_vec(:),Y_vec(:));
    str1(1) = {['n            : ' num2str(stat.output.numobs,'%6.0f')]};
    str1(2) = {['Slope      : ',num2str(stat.curve.p1,'%6.2f')]};
    str1(3) = {['Y-intercept  : ',num2str(stat.curve.p2,'%6.2f')]};
    str1(4) = {['R^2        : ',num2str(stat.gof.rsquare,'%6.2f')]};
    str1(5) = {['bias        : ',num2str((nanmean(X_vec(:)-Y_vec(:))),'%6.2f')]};
    str1(6) = {['std        : ',num2str((nanstd(X_vec(:)-Y_vec(:))),'%6.2f')]};
    ht      = text(11,31,str1);
    set(ht,'fontsize',28,'FontName','times','fontweight','b')
    ax(1) = gca; set(ax(1),'FontSize',24,'LineWidth',4)
    xlabel('fwc [m] In situ ','Fontsize',28);
    ylabel('fwc [m] SMOS 16m TOPAZ ','Fontsize',28);
end

% TOPAZ+SMOS 25m

%X_vec = FWC_situ(:); Y_vec = FWC_smos25_int(:);
for i = 1:9
    X_vec = FWC_situ(:,i); X_vec = X_vec(:); Y_vec = FWC_smos25_int(:,i);Y_vec = Y_vec(:);

    nn = find(~isnan(X_vec) & ~isnan(Y_vec));
    X_vec = X_vec(nn);Y_vec = Y_vec(nn);

    figure;set(gcf,'position',[400 1200 800 800],'color','w');
    set(gca,'fontsize',16,'boxStyle','full')
    scatter(X_vec(:),Y_vec(:),26,[0.5 0.5 0.5],'filled');hold on;
    xlim([10 35]);ylim([10 35])
    plot([10 35],[10 35],'--r','LineWidth',3)
    hold on
    [stat.curve, stat.gof,stat.output] = fit(double(X_vec(:)), double(Y_vec(:)), 'poly1' );
    ci = predint(stat.curve,Longitude,0.99);
    p1=plot(stat.curve,'--k');set(p1,'LineWidth',3);legend off;grid on
    [R p] = corrcoef(X_vec(:),Y_vec(:));
    str1(1) = {['n            : ' num2str(stat.output.numobs,'%6.0f')]};
    str1(2) = {['Slope      : ',num2str(stat.curve.p1,'%6.2f')]};
    str1(3) = {['Y-intercept  : ',num2str(stat.curve.p2,'%6.2f')]};
    str1(4) = {['R^2        : ',num2str(stat.gof.rsquare,'%6.2f')]};
    str1(5) = {['bias        : ',num2str((nanmean(X_vec(:)-Y_vec(:))),'%6.2f')]};
    str1(6) = {['std        : ',num2str((nanstd(X_vec(:)-Y_vec(:))),'%6.2f')]};
    ht      = text(11,31,str1);
    set(ht,'fontsize',28,'FontName','times','fontweight','b')
    ax(1) = gca; set(ax(1),'FontSize',24,'LineWidth',4)
    xlabel('fwc [m] In situ ','Fontsize',28);
    ylabel('fwc [m] SMOS 25m TOPAZ ','Fontsize',28);
end

% TOPAZ+SMOS 29m
%X_vec = FWC_situ(:); Y_vec = FWC_smos29_int(:);
for i = 1:9
    X_vec = FWC_situ(:,i); X_vec = X_vec(:); Y_vec = FWC_smos29_int(:,i);Y_vec = Y_vec(:);

    nn = find(~isnan(X_vec) & ~isnan(Y_vec));
    X_vec = X_vec(nn);Y_vec = Y_vec(nn);

    figure;set(gcf,'position',[400 1200 800 800],'color','w');
    set(gca,'fontsize',16,'boxStyle','full')
    scatter(X_vec(:),Y_vec(:),26,[0.5 0.5 0.5],'filled');hold on;
    xlim([10 35]);ylim([10 35])
    plot([10 35],[10 35],'--r','LineWidth',3)
    hold on
    [stat.curve, stat.gof,stat.output] = fit(double(X_vec(:)), double(Y_vec(:)), 'poly1' );
    ci = predint(stat.curve,Longitude,0.99);
    p1=plot(stat.curve,'--k');set(p1,'LineWidth',3);legend off;grid on
    [R p] = corrcoef(X_vec(:),Y_vec(:));
    str1(1) = {['n            : ' num2str(stat.output.numobs,'%6.0f')]};
    str1(2) = {['Slope      : ',num2str(stat.curve.p1,'%6.2f')]};
    str1(3) = {['Y-intercept  : ',num2str(stat.curve.p2,'%6.2f')]};
    str1(4) = {['R^2        : ',num2str(stat.gof.rsquare,'%6.2f')]};
    str1(5) = {['bias        : ',num2str((nanmean(X_vec(:)-Y_vec(:))),'%6.2f')]};
    str1(6) = {['std        : ',num2str((nanstd(X_vec(:)-Y_vec(:))),'%6.2f')]};
    ht      = text(11,31,str1);
    set(ht,'fontsize',28,'FontName','times','fontweight','b')
    ax(1) = gca; set(ax(1),'FontSize',24,'LineWidth',4)
    xlabel('fwc [m] In situ ','Fontsize',28);
    ylabel('fwc [m] SMOS 29m TOPAZ ','Fontsize',29);
end


% 4. Hacer el scatterplot con solo los puntos donde hay smos:

% Interpolate TOPAZ+SMOS estimates to in-situ yearly (jul-oct) and 50*50 km sapling
i = 1;
for tt = 2011:2019
    kk = find(year(X3)==tt & ( month(X3)==8 | month(X3)==9 | month(X3)==10));
    BG_topaz_fwc_mean(:,:,i) = nanmean(BG_topaz_fwc_smosmask(:,:,kk),3);
    BG_smos16_fwc_mean(:,:,i) = nanmean(BG_topaz16_fwc_smosmask(:,:,kk),3);
    BG_smos25_fwc_mean(:,:,i) = nanmean(BG_topaz25_fwc_smosmask(:,:,kk),3);
    BG_smos29_fwc_mean(:,:,i) = nanmean(BG_topaz29_fwc_smosmask(:,:,kk),3);

    FWC_topaz_int(:,i) = interpn(xx2,yy2,squeeze(BG_topaz_fwc_mean(:,:,i)),Longitude, Latitude);
    FWC_smos16_int(:,i) = interpn(xx2,yy2,squeeze(BG_smos16_fwc_mean(:,:,i)),Longitude, Latitude);
    FWC_smos25_int(:,i) = interpn(xx2,yy2,squeeze(BG_smos25_fwc_mean(:,:,i)),Longitude, Latitude);
    FWC_smos29_int(:,i) = interpn(xx2,yy2,squeeze(BG_smos29_fwc_mean(:,:,i)),Longitude, Latitude);

    i = i +1;
end

% TOPAZ Only SMOS MASK
X_vec = FWC_situ(:); Y_vec = FWC_topaz_int(:);
nn = find(~isnan(X_vec) & ~isnan(Y_vec)); 
X_vec = X_vec(nn);Y_vec = Y_vec(nn);

figure;set(gcf,'position',[400 1200 800 800],'color','w');
set(gca,'fontsize',16,'boxStyle','full')
scatter(X_vec(:),Y_vec(:),26,[0.5 0.5 0.5],'filled');hold on;
xlim([10 35]);ylim([10 35])
plot([10 35],[10 35],'--r','LineWidth',3)
hold on
[stat.curve, stat.gof,stat.output] = fit(double(X_vec(:)), double(Y_vec(:)), 'poly1' );
ci = predint(stat.curve,Longitude,0.99);
p1=plot(stat.curve,'--k');set(p1,'LineWidth',3);legend off;grid on
[R p] = corrcoef(X_vec(:),Y_vec(:));
str1(1) = {['n            : ' num2str(stat.output.numobs,'%6.0f')]};
str1(2) = {['Slope      : ',num2str(stat.curve.p1,'%6.2f')]};
str1(3) = {['Y-intercept  : ',num2str(stat.curve.p2,'%6.2f')]};
str1(4) = {['R^2        : ',num2str(stat.gof.rsquare,'%6.2f')]};
str1(5) = {['bias        : ',num2str((nanmean(X_vec(:)-Y_vec(:))),'%6.2f')]};
str1(6) = {['std        : ',num2str((nanstd(X_vec(:)-Y_vec(:))),'%6.2f')]};
ht      = text(11,31,str1);
set(ht,'fontsize',28,'FontName','times','fontweight','b')
ax(1) = gca; set(ax(1),'FontSize',24,'LineWidth',4)
xlabel('fwc [m] In situ ','Fontsize',28);
ylabel('fwc [m] TOPAZ (smos mask)','Fontsize',28);

% TOPAZ+SMOS 16m SMOS MASK
X_vec = FWC_situ(:); Y_vec = FWC_smos16_int(:);
nn = find(~isnan(X_vec) & ~isnan(Y_vec)); 
X_vec = X_vec(nn);Y_vec = Y_vec(nn);

figure;set(gcf,'position',[400 1200 800 800],'color','w');
set(gca,'fontsize',16,'boxStyle','full')
scatter(X_vec(:),Y_vec(:),26,[0.5 0.5 0.5],'filled');hold on;
xlim([10 35]);ylim([10 35])
plot([10 35],[10 35],'--r','LineWidth',3)
hold on
[stat.curve, stat.gof,stat.output] = fit(double(X_vec(:)), double(Y_vec(:)), 'poly1' );
ci = predint(stat.curve,Longitude,0.99);
p1=plot(stat.curve,'--k');set(p1,'LineWidth',3);legend off;grid on
[R p] = corrcoef(X_vec(:),Y_vec(:));
str1(1) = {['n            : ' num2str(stat.output.numobs,'%6.0f')]};
str1(2) = {['Slope      : ',num2str(stat.curve.p1,'%6.2f')]};
str1(3) = {['Y-intercept  : ',num2str(stat.curve.p2,'%6.2f')]};
str1(4) = {['R^2        : ',num2str(stat.gof.rsquare,'%6.2f')]};
str1(5) = {['bias        : ',num2str((nanmean(X_vec(:)-Y_vec(:))),'%6.2f')]};
str1(6) = {['std        : ',num2str((nanstd(X_vec(:)-Y_vec(:))),'%6.2f')]};
ht      = text(11,31,str1);
set(ht,'fontsize',28,'FontName','times','fontweight','b')
ax(1) = gca; set(ax(1),'FontSize',24,'LineWidth',4)
xlabel('fwc [m] In situ ','Fontsize',28);
ylabel('fwc [m] SMOS 16m TOPAZ (smos mask)','Fontsize',28);

% TOPAZ+SMOS 25m SMOS MASK

X_vec = FWC_situ(:); Y_vec = FWC_smos25_int(:);
nn = find(~isnan(X_vec) & ~isnan(Y_vec)); 
X_vec = X_vec(nn);Y_vec = Y_vec(nn);

figure;set(gcf,'position',[400 1200 800 800],'color','w');
set(gca,'fontsize',16,'boxStyle','full')
scatter(X_vec(:),Y_vec(:),26,[0.5 0.5 0.5],'filled');hold on;
xlim([10 35]);ylim([10 35])
plot([10 35],[10 35],'--r','LineWidth',3)
hold on
[stat.curve, stat.gof,stat.output] = fit(double(X_vec(:)), double(Y_vec(:)), 'poly1' );
ci = predint(stat.curve,Longitude,0.99);
p1=plot(stat.curve,'--k');set(p1,'LineWidth',3);legend off;grid on
[R p] = corrcoef(X_vec(:),Y_vec(:));
str1(1) = {['n            : ' num2str(stat.output.numobs,'%6.0f')]};
str1(2) = {['Slope      : ',num2str(stat.curve.p1,'%6.2f')]};
str1(3) = {['Y-intercept  : ',num2str(stat.curve.p2,'%6.2f')]};
str1(4) = {['R^2        : ',num2str(stat.gof.rsquare,'%6.2f')]};
str1(5) = {['bias        : ',num2str((nanmean(X_vec(:)-Y_vec(:))),'%6.2f')]};
str1(6) = {['std        : ',num2str((nanstd(X_vec(:)-Y_vec(:))),'%6.2f')]};
ht      = text(11,31,str1);
set(ht,'fontsize',28,'FontName','times','fontweight','b')
ax(1) = gca; set(ax(1),'FontSize',24,'LineWidth',4)
xlabel('fwc [m] In situ ','Fontsize',28);
ylabel('fwc [m] SMOS 25m TOPAZ (smos mask)','Fontsize',28);

% TOPAZ+SMOS 29m SMOS MASK
X_vec = FWC_situ(:); Y_vec = FWC_smos29_int(:);
nn = find(~isnan(X_vec) & ~isnan(Y_vec)); 
X_vec = X_vec(nn);Y_vec = Y_vec(nn);

figure;set(gcf,'position',[400 1200 800 800],'color','w');
set(gca,'fontsize',16,'boxStyle','full')
scatter(X_vec(:),Y_vec(:),26,[0.5 0.5 0.5],'filled');hold on;
xlim([10 35]);ylim([10 35])
plot([10 35],[10 35],'--r','LineWidth',3)
hold on
[stat.curve, stat.gof,stat.output] = fit(double(X_vec(:)), double(Y_vec(:)), 'poly1' );
ci = predint(stat.curve,Longitude,0.99);
p1=plot(stat.curve,'--k');set(p1,'LineWidth',3);legend off;grid on
[R p] = corrcoef(X_vec(:),Y_vec(:));
str1(1) = {['n            : ' num2str(stat.output.numobs,'%6.0f')]};
str1(2) = {['Slope      : ',num2str(stat.curve.p1,'%6.2f')]};
str1(3) = {['Y-intercept  : ',num2str(stat.curve.p2,'%6.2f')]};
str1(4) = {['R^2        : ',num2str(stat.gof.rsquare,'%6.2f')]};
str1(5) = {['bias        : ',num2str((nanmean(X_vec(:)-Y_vec(:))),'%6.2f')]};
str1(6) = {['std        : ',num2str((nanstd(X_vec(:)-Y_vec(:))),'%6.2f')]};
ht      = text(11,31,str1);
set(ht,'fontsize',28,'FontName','times','fontweight','b')
ax(1) = gca; set(ax(1),'FontSize',24,'LineWidth',4)
xlabel('fwc [m] In situ ','Fontsize',28);
ylabel('fwc [m] SMOS 29m TOPAZ (smos mask)','Fontsize',29);


% and 130??W - 170??W where water depths exceed 300 m 