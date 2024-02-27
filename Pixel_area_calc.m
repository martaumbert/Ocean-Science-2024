%%%%%% SCRIPT TO CALCULATE PIXEL AREA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% 1. Load the file
% 2. Select the Earth's Radius & excentricity
% 3. Calculate de area in km2
% 4. Save the result as a matlab files and csv file
% -------------------------------------------------------------------------
% Author    : Eva de Andres
% Email     : eva.deandres@upm.es
% Creation  : 19/12/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;       % clear Command Window
clear;     % clear Workspace
close all  % Close all figures

%addpath('/mnt/lustre/users/BECPolar/analysis/fwc/');
path = 'D:/Arctic/data/fwc/';
path_save = 'D:/Arctic/data/fwc/';
filename = ([path,'smos_topaz4b_monthly_fullArctic_fwc2011_2019_1m.nc']);

% Lat-lon data
yy = ncread(filename, 'latitude');  % common to all data sets
xx = ncread(filename, 'longitude'); % common to all data sets

% Earth's Radius & excentricity
Rearth=6356; % in km, so the area, A, will be in km2
Eearth=0.017;

% extract vectors
lons = xx(:,1);
lats = yy(:,1);

area = zeros(length(lats), length(lons));

for i=1:length(lats)
    for j= 1:length(lons)
        if (j <length(lons) && i<length(lats))
            LAT=[lats(i), lats(i+1),lats(i+1),lats(i)];
            LON=[lons(j), lons(j), lons(j+1),lons(j+1)];
            area(i,j)=areaint(LAT,LON,[Rearth, Eearth]);

        else
            LAT=[lats(end), lats(end),lats(end),lats(end)];
            LON=[lons(end), lons(end), lons(end),lons(end)];
        end
    end
end

area(:,j)=area(:,1); % Last coumnn filled qith the values of the first one (note that area only depends on latitude)

%cd('/mnt/lustre/users/BECPolar/analysis/fwc/')
save((fullfile(path_save,'area.mat')),'area')
save(([path_save,'area.mat']),area)
%plot_BG(xx,yy,area,'area [km2]','haline',bx,by);

% Save as csv ----------------
writematrix(area,([path_save,'area2_v2.csv']))
