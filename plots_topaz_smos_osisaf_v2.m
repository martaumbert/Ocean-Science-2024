%% Script for the study of Topaz-FWC, SMOS-SSS & OSISAF in the Beaufort Gyre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% 1. Load all data products & apply BG Mask
% 2. Calculate monthly means for topaz, smos, osisaf & smos-topaz variables
% 3. FWC vs SeaIce Coverage%
% 4. Temporal evolution of SMOS-SSS & Topaz-FWC
% -------------------------------------------------------------------------
% Author   : Marta Umbert from plots_topaz_smos_osisaf
% Email    : mumbert@icm.csic.es 
% Creation  : 26/04/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add folders and subfolders of the current path
%addpath('/mnt/lustre/users/BECPolar/analysis/fwc/outputs_smos_topaz/');
% study years
period=(2011:2019);

% area of each pixel
load('area.mat');

% Beaufort Gyre limits (bx by)
load poligon_bg.mat

%% 1. Load all data products
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% topaz data
filename = 'topaz4b_monthly_fullArctic_fwc2011_2019_latlon_vSmean.nc';

yy = ncread(filename, 'latitude'); %yy = yy(125:end,:); % common to all data sets
xx = ncread(filename, 'longitude');%xx = xx(125:end,:); % common to all data sets
[xx2,yy2] = meshgrid(xx,yy);
xx2 = permute(xx2,[2 1]);
yy2 = permute(yy2,[2 1]);
time = ncread(filename, 'time');

topaz.fwc  = ncread(filename, 'fwc');%topaz.fwc=topaz.fwc(125:end,:,:); % delimited to the same data set as smos & osisaf
topaz.sss  = ncread(filename, 'sss');%topaz.sss=topaz.sss(125:end,:,:);
topaz.time = ncread(filename, 'time');
% Topaz time and lat lon
X = topaz.time/(60*60*24)+datenum(1970,1,1);
X2 = datestr(X);
X3 = datetime(X2);

%clear addpath
area = permute(area,[2 1]);

for i=1:108
    topaz.fwc_km3(:,:,i) =(topaz.fwc(:,:,i)*0.001).*area;% fwc in km3/pixel
end
% smos data
%addpath('/mnt/lustre/users/BECPolar/data/regrid_cdo/artic_sss_smos_bec/')
%smos.sss = load('smos_mercator025_daily_beaufort_sss.mat');
filename2 = 'sss_smos_bec_2011_2019_regrided_025_monthly_data_40N.nc';
smos.sss = ncread(filename2, 'sss');
smos.lat = ncread(filename2, 'lat');
smos.lon = ncread(filename2, 'lon');
[yy2,xx2] = meshgrid(smos.lat,smos.lon);

filename4 = 'sss_smos_bec_2011_2019_monthly_data.nc';
smos.time = ncread(filename4, 'time');
zinterval = datetime('01-Jan-1970 00:00:00');
sss_tt = zinterval + seconds(smos.time);
%clear addpath

% osisaf data
%osisaf = load('osisaf_mercator025_daily_beaufort_fice.mat');
filename3 = 'sic_osisaf_2011_2022_regrided_025_monthly_data.nc'; 
osisaf.fice = ncread(filename3, 'sic');

for i=1:108
    osisaf.fice_km2(:,:,i) =osisaf.fice(:,:,i).*area; % iced area/pixel
end

%% Plot some maps

% SMOS - OSISAF Sep 2011

figure 
set(gcf, 'position',[149 108 900 900], 'color', 'w');
set(gca,'fontsize',28) 
m_proj('Azimuthal Equal-area','lat',74,'long',-150,'radius',12, 'rectbox','on'); 
m_pcolor(xx2,yy2,smos.sss(:,:,9));shading flat;hold on;colormap('jet');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
m_gshhs_l('color','b');
m_grid('xtick',20,'tickdir','in','ytick',[50 60 70 80],'linest','-', 'linewidth', 1); %Para comparacion 
hold on 
cc2 = [0.5,0.5]; cc3 = [0.1,0.1]; cc4 = [0.3,0.3];cc5 = [0.7,0.7];
[C, h] =m_contour(xx2,yy2,osisaf.fice(:,:,9),cc3,'color',rgb('dimgray'),'Linewidth',3);clabel(C, h, 'FontSize', 20, 'Color', 'k');
[C, h] =m_contour(xx2,yy2,osisaf.fice(:,:,9),cc4,'color',rgb('slategray'),'Linewidth',3);clabel(C, h, 'FontSize', 20, 'Color', 'k');
[C, h] =m_contour(xx2,yy2,osisaf.fice(:,:,9),cc2,'color',rgb('silver'),'Linewidth',3);clabel(C, h, 'FontSize', 20, 'Color', 'k');
[C, h] =m_contour(xx2,yy2,osisaf.fice(:,:,9),cc5,'color',rgb('lightgray'),'Linewidth',3);clabel(C, h, 'FontSize', 20, 'Color', 'k');
c = colorbar('fontsize',32,'location','southoutside');
c.Label.String = 'psu ';
caxis([22 32]);
m_plot(bx,by,'--k','Linewidth',2)

% SMOS - OSISAF Sep 2016

figure 
set(gcf, 'position',[149 108 900 900], 'color', 'w');
set(gca,'fontsize',28) 
m_proj('Azimuthal Equal-area','lat',74,'long',-150,'radius',12, 'rectbox','on'); 
m_pcolor(xx2,yy2,smos.sss(:,:,69));shading flat;hold on;colormap('jet');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
m_gshhs_l('color','b');
m_grid('xtick',20,'tickdir','in','ytick',[50 60 70 80],'linest','-', 'linewidth', 1); %Para comparacion 
hold on 
cc2 = [0.5,0.5]; cc3 = [0.1,0.1]; cc4 = [0.3,0.3];cc5 = [0.7,0.7];
[C, h] =m_contour(xx2,yy2,osisaf.fice(:,:,69),cc3,'color',rgb('dimgray'),'Linewidth',3);clabel(C, h, 'FontSize',20, 'Color', 'k');
[C, h] =m_contour(xx2,yy2,osisaf.fice(:,:,69),cc4,'color',rgb('slategray'),'Linewidth',3);clabel(C, h, 'FontSize', 20, 'Color', 'k');
[C, h] =m_contour(xx2,yy2,osisaf.fice(:,:,69),cc2,'color',rgb('silver'),'Linewidth',3);clabel(C, h, 'FontSize', 20, 'Color', 'k');
[C, h] =m_contour(xx2,yy2,osisaf.fice(:,:,69),cc5,'color',rgb('lightgray'),'Linewidth',3);clabel(C, h, 'FontSize', 20, 'Color', 'k');
c = colorbar('fontsize',32,'location','southoutside');
c.Label.String = 'psu ';
caxis([22 32]);
m_plot(bx,by,'--k','Linewidth',2)

% TOPAZ - OSISAF Sep 2011

figure 
set(gcf, 'position',[149 108 900 900], 'color', 'w');
set(gca,'fontsize',28) 
m_proj('Azimuthal Equal-area','lat',74,'long',-150,'radius',12, 'rectbox','on'); 
m_pcolor(xx2,yy2,topaz.sss(:,:,9));shading flat;hold on;colormap('jet');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
m_gshhs_l('color','b');
m_grid('xtick',20,'tickdir','in','ytick',[50 60 70 80],'linest','-', 'linewidth', 1); %Para comparacion 
hold on 
cc2 = [0.5,0.5]; cc3 = [0.1,0.1]; cc4 = [0.3,0.3];cc5 = [0.7,0.7];
[C, h] =m_contour(xx2,yy2,osisaf.fice(:,:,9),cc3,'color',rgb('dimgray'),'Linewidth',3);clabel(C, h, 'FontSize', 20, 'Color', 'k');
[C, h] =m_contour(xx2,yy2,osisaf.fice(:,:,9),cc4,'color',rgb('slategray'),'Linewidth',3);clabel(C, h, 'FontSize', 20, 'Color', 'k');
[C, h] =m_contour(xx2,yy2,osisaf.fice(:,:,9),cc2,'color',rgb('silver'),'Linewidth',3);clabel(C, h, 'FontSize',20, 'Color', 'k');
[C, h] =m_contour(xx2,yy2,osisaf.fice(:,:,9),cc5,'color',rgb('lightgray'),'Linewidth',3);clabel(C, h, 'FontSize',20, 'Color', 'k');
c = colorbar('fontsize',32,'location','southoutside');
c.Label.String = 'psu ';
caxis([22 32]);
m_plot(bx,by,'--k','Linewidth',2)

% TOPAZ - OSISAF Sep 2016

figure 
set(gcf, 'position',[149 108 900 900], 'color', 'w');
set(gca,'fontsize',28) 
m_proj('Azimuthal Equal-area','lat',74,'long',-150,'radius',12, 'rectbox','on'); 
m_pcolor(xx2,yy2,topaz.sss(:,:,69));shading flat;hold on;colormap('jet');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
m_gshhs_l('color','b');
m_grid('xtick',20,'tickdir','in','ytick',[50 60 70 80],'linest','-', 'linewidth', 1); %Para comparacion 
hold on 
cc2 = [0.5,0.5]; cc3 = [0.1,0.1]; cc4 = [0.3,0.3];cc5 = [0.7,0.7];
[C, h] =m_contour(xx2,yy2,osisaf.fice(:,:,69),cc3,'color',rgb('dimgray'),'Linewidth',3);clabel(C, h, 'FontSize', 20, 'Color', 'k');
[C, h] =m_contour(xx2,yy2,osisaf.fice(:,:,69),cc4,'color',rgb('slategray'),'Linewidth',3);clabel(C, h, 'FontSize', 20, 'Color', 'k');
[C, h] =m_contour(xx2,yy2,osisaf.fice(:,:,69),cc2,'color',rgb('silver'),'Linewidth',3);clabel(C, h, 'FontSize', 20, 'Color', 'k');
[C, h] =m_contour(xx2,yy2,osisaf.fice(:,:,69),cc5,'color',rgb('lightgray'),'Linewidth',3);clabel(C, h, 'FontSize', 20, 'Color', 'k');
c = colorbar('fontsize',32,'location','southoutside');
c.Label.String = 'psu ';
caxis([22 32]);
m_plot(bx,by,'--k','Linewidth',2)


%% 2. Monthly means
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrix to be filled
% topaz --> not necesary for loaded data, already monthly means

% smos  
%%% OJO ahora ya leemos monthly smos sss
%smos.sss_mm_std=zeros([size(smos.sss(:,:,1)),length(topaz.time)]);
%smos.sss_mm_mean=zeros([size(smos.sss(:,:,1)),length(topaz.time)]);
%smos.sss_mean = zeros(length(topaz.time),1);
%smos.sss_std = zeros(length(topaz.time),1);

% osisaf  
%%% OJO leemos monthly osisaf
%fice_free_mm=zeros([size(osisaf.fice(:,:,1)),length(topaz.time)]);
%fice_days_mm=zeros([size(osisaf.fice(:,:,1)),length(topaz.time)]);
%osisaf.fice_mm_km2=zeros([size(osisaf.fice(:,:,1)),length(topaz.time)]);
%oissaf.fice_mean

% DIF. smos-topaz
sss_dif=zeros(size(smos.sss));
sss_dif=smos.sss-topaz.sss;

%tt = 1;% general counter

%for i = 1:length(period)
    
%    for j =1:12
        % topaz 
               
        % smos
 %       mm = find((month(smos.time) == j) & (year (smos.time) == period(i)));
 %       smos.sss_mm_std(:,:,tt)=std(smos.sss(:,:,mm),0,3,'omitnan');
 %       smos.sss_mm_mean(:,:,tt)=mean(smos.sss(:,:,mm),3,'omitnan');

        % osisaf
  %      mm = find((month(osisaf.time) == j) & (year (osisaf.time) == period(i)));
  %      fice_positive = osisaf.fice(:,:,mm) > 0;
  %      fice_filtered = osisaf.fice(:,:,mm);
        
   %     fice_filtered(fice_positive) = nan;
        
   %     fice_free_mm(:,:,tt) = mean(fice_filtered,3,'omitnan');
   %     fice_days_mm(:,:,tt) = sum(fice_positive,3,'omitnan');
        
    %    osisaf.fice_mm_km2 (:,:,tt) = mean(osisaf.fice_km2(:,:,mm),3,'omitnan');   
       
        % DIF. smos-topaz
 %       sss_dif(:,:,tt)=smos.sss(:,:,tt)-topaz.sss(:,:,tt);

 %       tt = tt +1;% general counter
 %   end
    
%end

%%%%%%%%%
% PLOTs %
%%%%%%%%%
% optimizing subplot sizes
% load pos.mat
% 
% % plotting Sept months (max. data coverage)
% sept = 9;
% 
% for i=1:length(period)
%     mes = sept+((i-1)*12);
% %     figure(1);set(gcf,'position',[249 608 780 767],'color','w');subplot(2,5,i);set(subplot(2,5,i),'position',[pos(i,:)], 'color','w');
% %         plot_BG_subp(xx,yy,topaz.sss(:,:,mes),'SSS (psu)','haline',21,32,bx,by);title(['Sept ',num2str(period(i))]);
% %     figure(2);set(gcf,'position',[249 608 780 767],'color','w');subplot(2,5,i);set(subplot(2,5,i),'position',[pos(i,:)], 'color','w');
% %         plot_BG_subp(xx,yy,smos.sss_mm_mean(:,:,mes),'SSS (psu)','haline',21, 32,bx,by);title(['Sept ',num2str(period(i))])
% %     figure(3);set(gcf,'position',[249 608 780 767],'color','w');subplot(2,5,i);set(subplot(2,5,i),'position',[pos(i,:)], 'color','w');
% %         plot_BG_subp(xx,yy,fice_days_mm(:,:,mes),'Number of days covered by sea ice','ice',0,30,bx,by);title(['Sept ',num2str(period(i))])
% %     figure(4);set(gcf,'position',[249 608 780 767],'color','w');subplot(2,5,i);set(subplot(2,5,i),'position',[pos(i,:)], 'color','w');
% %         plot_BG_subp(xx,yy,sss_dif(:,:,mes),'SSS_s_m_o_s - SSS_t_o_p_a_z (psu)','oxy',-2,2,bx,by);title(['Sept ',num2str(period(i))])
% %     figure(5);set(gcf,'position',[249 608 780 767],'color','w');subplot(2,5,i);set(subplot(2,5,i),'position',[pos(i,:)], 'color','w');
% %         plot_BG_subp(xx,yy,topaz.fwc(:,:,mes)/2.5,'FWC (m)','deep',10,30,bx,by);title(['Sept ',num2str(period(i))])
% %     
% end



% % PLOT: time evolution of full arctic sea ice extension, fwc and sss 
suma_mens_fwc=squeeze(sum(topaz.fwc_km3,[1 2],'omitnan'));

suma_mens_hielo=squeeze(sum(osisaf.fice_km2,[1 2],'omitnan'));

media_mens_sss = squeeze(mean(smos.sss, [1 2],'omitnan'));
std_mens_sss = squeeze(std(smos.sss,0,[1 2],'omitnan'));

media_mens_sss_topaz=squeeze(mean(topaz.sss,1, 'omitnan'));
media_mens_sss_topaz=squeeze(mean(media_mens_sss_topaz,1, 'omitnan'));


createplot_icefwc(1:108,suma_mens_hielo(1:108),1:108,suma_mens_fwc)

figure;set(gcf,'color','w','position',[249 608 1200 500]);
plot(1:108,suma_mens_hielo,'Linewidth',2,'Color',rgb('skyblue'))
xlabel('time(months since 2011)')
ylabel('mean SIC in Arctic')
grid on;set(gca,'Fontsize',24)

figure;set(gcf,'color','w','position',[249 608 1200 500]);
plot(1:108,suma_mens_fwc,'Linewidth',2,'Color',rgb('darkviolet'))
xlabel('time(months since 2011)')
ylabel('mean FWC in Arctic')
grid on;set(gca,'Fontsize',24)

figure;set(gcf,'color','w','position',[249 608 1200 500]);
plot(1:108,media_mens_sss,'Linewidth',2,'Color',rgb('green'))
xlabel('time(months since 2011)')
ylabel('mean SSS SMOS in Arctic')
grid on;set(gca,'Fontsize',24)

figure;set(gcf,'color','w','position',[249 608 1200 500]);
plot(1:108,media_mens_sss_topaz,'Linewidth',2,'Color',rgb('darkgreen'))
hold on
xlabel('time(months since 2011)')
ylabel('mean SSS TOPAZ in Arctic')
grid on;set(gca,'Fontsize',24)


%%%% shadow for the STD
yyaxis right
hold on
x=1:108;
x2 = [x, fliplr(x)];
curve1 = media_mens_sss(1:108) + std_mens_sss(1:108);
curve2 = media_mens_sss(1:108) - std_mens_sss(1:108);
inBetween = [curve1', fliplr(curve2')];
fill(x2,inBetween,'r')

yyaxis right
hold on;
plot(x,media_mens_sss(1:108),'Color',[0.85 0.325 0.098])
% correlation
figure
scatter(suma_mens_hielo(1:100),media_mens_sss(1:100));ylabel('SSS (psu');xlabel('Sea Ice ext (km²)');

%%%%%%%%%%
% VIDEOS %
%%%%%%%%%%
% sept=9;
% year_selct=[2011 2014 2015 2016];
% 
% for i = 1:length(year_selct)
%     mm = find((month(smos.time) == sept) & (year (smos.time) == year_selct(i)));
%     sss_smos_select=smos.sss(:,:,mm);
%     figure()
%     set(gcf,'position',[249 608 780 767],'color','w');
%     plot_BG_sept_video(xx,yy,sss_smos_select,'SSS (psu)','haline',bx,by,20,32, year_selct(i));
% end

% 3. FWC vs SeaIce Coverage & smos-sss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%
% BG data selection
%%%%%%%%%%%%%%%%%%%%%%%%%
% creating the BG poligon
% 72, 80
% -127, -178

pgon = polyshape(bx,by); % bx and by have been loaded in 'poligon_bg.mat'

TFin=zeros(size(xx2));
TFon=zeros(size(xx2));

for i=1:length(xx2(:,1))
    [TFin(i,:), TFon(i,:)] = isinterior(pgon,xx2(i,:),yy2(i,:));%en TFin guardamos los puntos dentro del polígono (incluyendo la frontera), en TFon identificamos sólo los puntos en la frontera
end
% BG Mask for data selection
mask_bg = TFin;
mask_bg(mask_bg ==0) = nan;
lon_bg=xx2.*mask_bg;
lat_bg=yy2.*mask_bg;

% topaz bg selection
BG.topaz_fwc = (topaz.fwc).*mask_bg;
BG.topaz_sss = (topaz.sss).*mask_bg;

% smos bg selection
BG.smos_sss = (smos.sss).*mask_bg;
%BG.smos_sss_mm_mean = (smos.sss_mm_mean).*mask_bg;
%BG.smos_sss_mm_std = (smos.sss_mm_std).*mask_bg;

% osisaf bg selection
BG.osisaf_fice = (osisaf.fice).*mask_bg;
BG.osisaf_fice_km2 = (osisaf.fice_km2).*mask_bg;

% SMOS-Mask for data selection
mask_smos = double(~isnan(BG.smos_sss));
mask_smos(mask_smos == 0) = NaN;

BG.topaz_sss_smosmask = BG.topaz_sss(:,:,108).*mask_smos;
BG.topaz_fwc_smosmask = BG.topaz_fwc(:,:,108).*mask_smos;

% % %PLOT: time evolution of BG sea ice extension, fwc and sss 
mean_mens_fwc = squeeze(mean(BG.topaz_fwc_smosmask,[1 2],'omitnan'));
suma_mens_hielo = squeeze(sum(BG.osisaf_fice_km2, [1 2],'omitnan'));
media_mens_sss_topaz = squeeze(mean(BG.topaz_sss_smosmask,[1 2], 'omitnan'));
media_mens_sss = squeeze(mean(BG.smos_sss, [1 2],'omitnan'));
% %createplot_icefwc(1:108,suma_mens_hielo(1:108),1:108,suma_mens_fwc(1:108))
x = 1:108; x = x';
y1 = mean_mens_fwc(1:108);
y2 = media_mens_sss(1:108);
y3 = suma_mens_hielo(1:108);
y4 = media_mens_sss_topaz(1:108);
% FIGURE 3
figure
%plot(x,y1/2.5,x,y2, x,y4)
plot(x,y2, x,y4,'LineWidth',3)
%ylabel('SSS (psu)                fwc (m)');
ylabel('SSS (psu)');yaxis([20 30])
xlabel('time'); xticks(1:12:109)
yyaxis right; plot(x,y3,'LineWidth',3)
ylabel('sea ice extension (km²)')
xticklabels({'2011', '2012','2013','2014','2015','2016','2017','2018','2019','2020'})
%legend('fwc','SMOS SSS','TOPAZ SSS','SI extension','Location','South')
legend('SMOS SSS','TOPAZ SSS','Location','South')
set(gcf,'color','w','position',[249 608 1200 500]);
set(gca,'Fontsize',24)


% % % correlation omitting nan values
inc1=~isnan(y1);inc2=~isnan(y2);inc4=~isnan(y4);
[R,P] = corrcoef(y4(inc4),y1(inc1));
% % % 
figure(1);
% % TOPAZ
% subplot(1,3,1);
p = polyfit(y4(inc4), y1(inc1)/2.5, 1);
px = [min(y4(inc4)) max(y4(inc4))];
py = polyval(p, px);
scatter(y4(inc4), y1(inc1)/2.5, 'filled');ylabel('FWC (m)');xlabel('SSS_t_o_p_a_z (psu)');
hold on
plot(px, py, 'Color',[0.2 0.2 0.6],'LineWidth', 1);
hold off
% % SMOS
% subplot(1,3,2);
[R,P] = corrcoef(y2(inc2),y1(inc1))
p = polyfit(y2(inc2), y1(inc1)/2.5, 1);
px = [min(y2(inc2)) max(y2(inc2))];
py = polyval(p, px);
figure(1);hold on
scatter(y2(inc2), y1(inc1)/2.5, 'filled');ylabel('FWC (m)');xlabel('SSS_s_m_o_s (psu)');
hold on
plot(px, py, 'Color',[0.2 0.2 0.6],'LineWidth', 1);
hold off
% %
% subplot(1,3,3);
% p = polyfit(y2(inc2), y4(inc4), 1);
% px = [min(y2(inc2)) max(y2(inc2))];
% py = polyval(p, px);
% scatter(y2(inc2), y4(inc4), 'filled');ylabel('SSS_t_o_p_a_z (psu)');xlabel('SSS_s_m_o_s (psu)');
% hold on
% plot(px, py, 'Color',[0.2 0.2 0.6],'LineWidth', 1);
% hold off
% % scatter(y4(inc4),y2(inc2));ylabel('FWC (m)');xlabel('SSS_t_o_p_a_z (psu)');

%%%%%% Daily evolution of Sea Ice and SSS
sept=9;july = 7;aug=8;
year_selct=2011:2019;
%year_selct=[2011 2014 2016 2018];
clear sss iced_area


for i = 1:length(year_selct)
    nn = find((month(sss_tt) == sept |month(sss_tt) == july | month(sss_tt) == aug) & (year (sss_tt) == year_selct(i)));
    %nn = find((year (sss_tt) == year_selct(i)));

    sss(i)=mean(BG.smos_sss(:,:,nn),[1 2 3],'omitnan');
    iced_area(i)=sum(BG.osisaf_fice(:,:,nn).*area,[1 2 3],'omitnan');
end


figure;set(gcf,'color','w','position',[249 608 1200 500]);
plot(year_selct,iced_area,'Linewidth',2,'Color',rgb('darkviolet'))
xlabel('time(months since 2011)')
ylabel('mean FWC in Arctic')
grid on;set(gca,'Fontsize',24)
