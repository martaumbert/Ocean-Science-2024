%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% 1. Read in-situ data from Beaufort Gyre Exploration Project
% (https://www2.whoi.edu/site/beaufortgyre/data/freshwater-content-gridded-data/
%freshwater-content-fwc-2003-2018-gridded-data/)
% 2. Read data from AVISO and EN4 regrided and 
%    results from SATCUR code 
% 3. Matchup 
% 4. Scatterplot
% -------------------------------------------------------------------------
% Author    : M. Umbert
% Email     : mumbert@icm.csic.es
% Creation  : 17/02/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 1. OPEN FWC from INSITU DATA
% All data are in meters. Uncertainties for each grid cell are determined using 
% the optimal interpolation technique previously described in Proshutinsky et al. (2009)


% Execute /fwc/scripts/open_bgfwconint.m
% Execute /fwc/scripts/open_bgfwconint_error.m

%% 2. Plot the data into a map

figure 
set(gcf, 'position',[149 108 780 567], 'color', 'w');
set(gca,'fontsize',16) 
m_proj('Azimuthal Equal-area','lat',74,'long',-150,'radius',12, 'rectbox','on'); 
m_scatter(Longitude,Latitude,100,fwc_2011,'filled');shading flat;hold on;colormap('jet');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
m_gshhs_l('color','b');
c = colorbar('fontsize',24,'location','southoutside');
m_grid('xtick',16,'tickdir','in','ytick',[50 60 70 80],'linest','-', 'linewidth', 2); %Para comparacion 
c.Label.String = 'fwc [m]';
caxis([10 30]);

figure 
set(gcf, 'position',[149 108 780 567], 'color', 'w');
set(gca,'fontsize',16) 
m_proj('Azimuthal Equal-area','lat',74,'long',-150,'radius',12, 'rectbox','on'); 
m_scatter(Longitude,Latitude,100,err_2011,'filled');shading flat;hold on;
m_coast('patch',[.7 .7 .7],'edgecolor','k');
m_gshhs_l('color','b');
c = colorbar('fontsize',24,'location','southoutside');
m_grid('xtick',16,'tickdir','in','ytick',[50 60 70 80],'linest','-', 'linewidth', 2); %Para comparacion 
c.Label.String = 'fwc [m]';
colormap(cmocean('dense'));
caxis([0 4]);

FWC_situ = [fwc_2011,fwc_2012,fwc_2013,fwc_2014,fwc_2015,fwc_2016,fwc_2017,fwc_2018,fwc_2019];

mm = [9,21,33,45,57,69,81,93,105];
for yy = 1:9
    FWC_top_int = interpn(xx2,yy2,fwc(:,:,mm(yy)),Longitude, Latitude);
    FWC_topsmos25_int = interpn(xx2,yy2,fwc_25(:,:,mm(yy)),Longitude, Latitude);

    diff_top(:,:,yy) = FWC_top_int - FWC_situ(:,yy);
    diff_topsmos(:,:,yy) = FWC_topsmos25_int - FWC_situ(:,yy);

end

for yy = 1:9

FWC_top_int = interpn(xx2,yy2,fwc(:,:,mm(yy)),Longitude, Latitude);
FWC_topsmos25_int = interpn(xx2,yy2,fwc_25(:,:,mm(yy)),Longitude, Latitude);
end

figure 
set(gcf, 'position',[149 108 780 567], 'color', 'w');
set(gca,'fontsize',16) 
m_proj('Azimuthal Equal-area','lat',74,'long',-150,'radius',12, 'rectbox','on'); 
m_scatter(Longitude,Latitude,100,FWC_top_int,'filled');shading flat;hold on;colormap('jet');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
m_gshhs_l('color','b');
c = colorbar('fontsize',24,'location','southoutside');
m_grid('xtick',24,'tickdir','in','ytick',[50 60 70 80],'linest','-', 'linewidth', 2); %Para comparacion 
c.Label.String = 'fwc [m]';
caxis([10 30]);

figure 
set(gcf, 'position',[149 108 780 567], 'color', 'w');
set(gca,'fontsize',16) 
m_proj('Azimuthal Equal-area','lat',74,'long',-150,'radius',12, 'rectbox','on'); 
m_scatter(Longitude,Latitude,100,FWC_topsmos25_int,'filled');shading flat;hold on;colormap('jet');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
m_gshhs_l('color','b');
c = colorbar('fontsize',16,'location','southoutside');
m_grid('xtick',24,'tickdir','in','ytick',[50 60 70 80],'linest','-', 'linewidth', 2); %Para comparacion 
c.Label.String = 'fwc [m]';
caxis([10 30]);


figure 
set(gcf, 'position',[149 108 780 567], 'color', 'w');
set(gca,'fontsize',16) 
m_proj('Azimuthal Equal-area','lat',74,'long',-150,'radius',12, 'rectbox','on'); 
m_scatter(Longitude,Latitude,100,mean(diff_top,3),'filled');shading flat;hold on;colormap(cmocean('balance'));
m_coast('patch',[.7 .7 .7],'edgecolor','k');
m_gshhs_l('color','b');
c = colorbar('fontsize',24,'location','southoutside');
m_grid('xtick',24,'tickdir','in','ytick',[50 60 70 80],'linest','-', 'linewidth', 2); %Para comparacion 
c.Label.String = 'fwc [m]';
caxis([-4 4]);

figure 
set(gcf, 'position',[149 108 780 567], 'color', 'w');
set(gca,'fontsize',16) 
m_proj('Azimuthal Equal-area','lat',74,'long',-150,'radius',12, 'rectbox','on'); 
m_scatter(Longitude,Latitude,100,mean(diff_topsmos,3),'filled');shading flat;hold on;colormap(cmocean('balance'));
m_coast('patch',[.7 .7 .7],'edgecolor','k');
m_gshhs_l('color','b');
c = colorbar('fontsize',24,'location','southoutside');
m_grid('xtick',24,'tickdir','in','ytick',[50 60 70 80],'linest','-', 'linewidth', 2); %Para comparacion 
c.Label.String = 'fwc [m]';
caxis([-4 4]);




[TSG_INTERP,TSG_RMS,xxvec,yyvec,ngrid] = ffgridrms(Longitude,Latitude,fwc_2011,0.25,0.25,min(smos.lon),min(smos.lat),max(smos.lon),max(smos.lat));
[xx2,yy2] = meshgrid(smos.lon,smos.lat);

SITU_INTERP = flipud(TSG_INTERP);
DIFF_INTERP = fwc_year_topaz(:,:,1) - SITU_INTERP';
DIFF_INTERP_16 = fwc_year_topaz16(:,:,1) - SITU_INTERP';
DIFF_INTERP_25 = fwc_year_topaz25(:,:,1) - SITU_INTERP';
DIFF_INTERP_29 = fwc_year_topaz29(:,:,1) - SITU_INTERP';

figure 
set(gcf, 'position',[149 108 900 900], 'color', 'w');
set(gca,'fontsize',36) 
m_proj('Azimuthal Equal-area','lat',74,'long',-150,'radius',12, 'rectbox','on'); 
m_pcolor(xx2,yy2,SITU_INTERP);shading flat;hold on;colormap('jet');
m_coast('patch',[.7 .7 .7],'edgecolor','k');
m_gshhs_l('color','b');
m_grid('xtick',20,'tickdir','in','ytick',[50 60 70 80],'linest','-', 'linewidth', 1); %Para comparacion 
hold on 
c = colorbar('fontsize',32,'location','southoutside');
c.Label.String = 'fwc [m]';
caxis([10 30]);
%m_plot(bx,by,'--k','Linewidth',2)

figure 
set(gcf, 'position',[149 108 900 900], 'color', 'w');
set(gca,'fontsize',36) 
m_proj('Azimuthal Equal-area','lat',74,'long',-150,'radius',12, 'rectbox','on'); 
m_pcolor(xx2,yy2,DIFF_INTERP');shading flat;hold on;colormap(cmocean('balance'));
m_coast('patch',[.7 .7 .7],'edgecolor','k');
m_gshhs_l('color','b');
m_grid('xtick',20,'tickdir','in','ytick',[50 60 70 80],'linest','-', 'linewidth', 1); %Para comparacion 
hold on 
c = colorbar('fontsize',32,'location','southoutside');
c.Label.String = 'fwc [m]';
caxis([-4 4]);

%% 3. Mean temporal evolution over the region
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

figure
set(gcf,'position',[249 608 1200 600],'color','w');
errorbar(years,mean_fwc,mean_err,'Linewidth',2,'Marker','diamond')
set(gca,'fontsize',26)
grid on; 
ylabel('FWC [m]')
