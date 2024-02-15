%%%%%% SCRIPT TO COMPUTE FWC from TOPAZ4B OUTPUTS (MLDm) %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% 1. List TOPAZ files
% 2. Load common variables once (lat, lon and depth)
% 3. Read SMOS SSS regrided in lustre
% 3. Calculate FWC (WARNING! dz = dz(z))
% 4. Create netcdf file
% -------------------------------------------------------------------------
% Author    : Eva De Andr?s & Marta Umbert & Maria Sanchez
% Email     : eva.deandres@upm.es / mumbert@icm.csic.es 
% Creation  : 21/12/2022
% Last modified : 21/04/2023 Marta Umbert & Maria S?nchez 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;       % clear Command Window
clear;     % clear Workspace
close all  % Close all figures

%--------------------------------------------------------------------------
% Specify the space where you are running this script :
gaia = true; % true = using gaia / false = using your personal computer

% 1. List TOPAZ files -----------------------------------------------------

% Path
if gaia==true
    addpath('/mnt/lustre/users/BECPolar/data/regrid_cdo/TOPAZ4/');
    pathfile = '/mnt/lustre/users/BECPolar/data/regrid_cdo/TOPAZ4/';
    path_smos = '/mnt/lustre/users/BECPolar/data/regrid_cdo/artic_sss_smos_bec/';
    path_save='/mnt/lustre/users/BECPolar/analysis/fwc/smos_topaz/data/';
else
    pathfile = 'D:/Arctic/data/topaz4b/mothly_files/';
    path_smos = 'D:/Arctic/data/sss_smos/regrid_025/';
    path_save='D:/Arctic/data/fwc/';
end

fileID = fopen([pathfile,'list_files.txt']);


% Each name file will be saved as a cell from the cell array 
C = textscan(fileID,'%s'); % %s'=read as a cell array of character vectors 
fclose(fileID);

% 2. Load common variables once (lat, lon and depth) ----------------------
% from the first file on list

filename = cell2mat(C{1}(1)); % convert cellarray to ordinary array
lat      = ncread(filename,'lat');
lon      = ncread(filename,'lon');
dep      = ncread(filename,'depth');

% 3. Load SMOS SSS (monthly data) ----------------------------------------- 

smosfile = fullfile(path_smos,'sss_smos_bec_2011_2019_regrided_025_monthly_data_40N.nc');

slat     = ncread(smosfile,'lat');
slon     = ncread(smosfile,'lon');
stime    = ncread(smosfile,'time');
ssss     = ncread(smosfile,'sss');

% Modify length of C (topaz files) - same as sss time
C{1} = C{1}(1:length(stime));

% 4. Calculate FWC (WARNING! dz = dz(z)) ----------------------------------
% establish reference Salinity
S0 = 34.8;

% Compute dz vector as a function of depth (ie. dz(z)). The distanct 
% between one level and the next one is not always the same, so here
% we compute the depth difference in each step:

ll = 1;
dz = zeros(1,length(dep)-1);

for l = 1:(length(dep)-1)
    dz(l) = dep(ll+1) - dep(l);
    ll = ll+1;
end

% Add 1 at the beginning of the matrix to have 40 levels and transpone dz
dz = [1; dz'];

% Estimate FWC for all files on list:

time = zeros(1,numel(C{1}));
fwc  = zeros(length(lon),length(lat),numel(C{1}));
sss  = zeros(length(lon),length(lat),numel(C{1}));
sst  = zeros(length(lon),length(lat),numel(C{1}));
sss_smos_topaz = zeros(length(lon),length(lat),numel(C{1}));

for it = 1:(numel(C{1})) % Total number of elements in cellarray  

    % 1/ Read salinity and temperature from topaz data
    filename = cell2mat(C{1}(it));
    filedate = datenum(filename(1:8),'yyyymmdd');
    time(it) = filedate;
    sal      = ncread(filename,'so');
    T        = ncread(filename,'thetao');

    % Depth 0 
    sal_l0 = sal(:,:,1);
    sal_l0(~isnan(ssss(:,:,it))) = NaN;
    tmp = cat(3,sal_l0,ssss(:,:,it));
    sal_l0_sum = sum(tmp,3,'omitnan');
    sal_l0_sum(sal_l0_sum==0) = NaN;
    topaz_smos = sal;
    topaz_smos(:,:,1) = sal_l0_sum;
    
    % Depth 2
    sal_l2 = topaz_smos(:,:,2);
    sal_l2(~isnan(ssss(:,:,it))) = NaN;
    tmp = cat(3,sal_l2,ssss(:,:,it));
    sal_l2_sum = sum(tmp,3,'omitnan');
    sal_l2_sum(sal_l2_sum==0) = NaN;
    topaz_smos(:,:,2) = sal_l2_sum;

    % Depth 3
    sal_l3 = topaz_smos(:,:,3);
    sal_l3(~isnan(ssss(:,:,it))) = NaN;
    tmp = cat(3,sal_l3,ssss(:,:,it));
    sal_l3_sum = sum(tmp,3,'omitnan');
    sal_l3_sum(sal_l3_sum==0) = NaN;
    topaz_smos(:,:,3) = sal_l3_sum;

    % Depth 4
    sal_l4 = topaz_smos(:,:,4);
    sal_l4(~isnan(ssss(:,:,it))) = NaN;
    tmp = cat(3,sal_l4,ssss(:,:,it));
    sal_l4_sum = sum(tmp,3,'omitnan');
    sal_l4_sum(sal_l4_sum==0) = NaN;
    topaz_smos(:,:,4) = sal_l4_sum;

    % Depth 5
    sal_l5 = topaz_smos(:,:,5);
    sal_l5(~isnan(ssss(:,:,it))) = NaN;
    tmp = cat(3,sal_l5,ssss(:,:,it));
    sal_l5_sum = sum(tmp,3,'omitnan');
    sal_l5_sum(sal_l5_sum==0) = NaN;
    topaz_smos(:,:,5) = sal_l5_sum;

    % Depth 6 
    sal_l6 = topaz_smos(:,:,6);
    sal_l6(~isnan(ssss(:,:,it))) = NaN;
    tmp = cat(3,sal_l6,ssss(:,:,it));
    sal_l6_sum = sum(tmp,3,'omitnan');
    sal_l6_sum(sal_l6_sum==0) = NaN;
    topaz_smos(:,:,6) = sal_l6_sum;

   % Depth 8
    sal_l8 = topaz_smos(:,:,7);
    sal_l8(~isnan(ssss(:,:,it))) = NaN;
    tmp = cat(3,sal_l8,ssss(:,:,it));
    sal_l8_sum = sum(tmp,3,'omitnan');
    sal_l8_sum(sal_l8_sum==0) = NaN;
    topaz_smos(:,:,7) = sal_l8_sum;

    % Delth 10
    sal_l10 = topaz_smos(:,:,8);
    sal_l10(~isnan(ssss(:,:,it))) = NaN;
    tmp = cat(3,sal_l10,ssss(:,:,it));
    sal_l10_sum = sum(tmp,3,'omitnan');
    sal_l10_sum(sal_l10_sum==0) = NaN;
    topaz_smos(:,:,8) = sal_l10_sum;

    % Depth 11
    sal_l11 = topaz_smos(:,:,9);
    sal_l11(~isnan(ssss(:,:,it))) = NaN;
    tmp = cat(3,sal_l11,ssss(:,:,it));
    sal_l11_sum = sum(tmp,3,'omitnan');
    sal_l11_sum(sal_l11_sum==0) = NaN;
    topaz_smos(:,:,9) = sal_l11_sum;

    % Depth 13 
    sal_l13 = topaz_smos(:,:,10);
    sal_l13(~isnan(ssss(:,:,it))) = NaN;
    tmp = cat(3,sal_l13,ssss(:,:,it));
    sal_l13_sum = sum(tmp,3,'omitnan');
    sal_l13_sum(sal_l13_sum==0) = NaN;
    topaz_smos(:,:,10) = sal_l13_sum;

    % Depth 16
    sal_l16 = topaz_smos(:,:,11);
    sal_l16(~isnan(ssss(:,:,it))) = NaN;
    tmp = cat(3,sal_l16,ssss(:,:,it));
    sal_l16_sum = sum(tmp,3,'omitnan');
    sal_l16_sum(sal_l16_sum==0) = NaN;
    topaz_smos(:,:,11) = sal_l16_sum;
    
    % Depth 18
    sal_l18 = topaz_smos(:,:,12);
    sal_l18(~isnan(ssss(:,:,it))) = NaN;
    tmp = cat(3,sal_l18,ssss(:,:,it));
    sal_l18_sum = sum(tmp,3,'omitnan');
    sal_l18_sum(sal_l18_sum==0) = NaN;
    topaz_smos(:,:,12) = sal_l18_sum;
    
    % Depth 22
    sal_l22 = topaz_smos(:,:,13);
    sal_l22(~isnan(ssss(:,:,it))) = NaN;
    tmp = cat(3,sal_l22,ssss(:,:,it));
    sal_l22_sum = sum(tmp,3,'omitnan');
    sal_l22_sum(sal_l22_sum==0) = NaN;
    topaz_smos(:,:,13) = sal_l22_sum;
    
    % Depth 25
    sal_l25 = topaz_smos(:,:,14);
    sal_l25(~isnan(ssss(:,:,it))) = NaN;
    tmp = cat(3,sal_l25,ssss(:,:,it));
    sal_l25_sum = sum(tmp,3,'omitnan');
    sal_l25_sum(sal_l25_sum==0) = NaN;
    topaz_smos(:,:,14) = sal_l25_sum;
    
    % Depth 29
    sal_l29 = topaz_smos(:,:,15);
    sal_l29(~isnan(ssss(:,:,it))) = NaN;
    tmp = cat(3,sal_l29,ssss(:,:,it));
    sal_l29_sum = sum(tmp,3,'omitnan');
    sal_l29_sum(sal_l29_sum==0) = NaN;
    topaz_smos(:,:,15) = sal_l29_sum;
    
    % Depth 34
    sal_l34 = topaz_smos(:,:,16);
    sal_l34(~isnan(ssss(:,:,it))) = NaN;
    tmp = cat(3,sal_l34,ssss(:,:,it));
    sal_l34_sum = sum(tmp,3,'omitnan');
    sal_l34_sum(sal_l34_sum==0) = NaN;
    topaz_smos(:,:,16) = sal_l34_sum;
    
    % Depth 40
    sal_l40 = topaz_smos(:,:,17);
    sal_l40(~isnan(ssss(:,:,it))) = NaN;
    tmp = cat(3,sal_l40,ssss(:,:,it));
    sal_l40_sum = sum(tmp,3,'omitnan');
    sal_l40_sum(sal_l40_sum==0) = NaN;
    topaz_smos(:,:,17) = sal_l40_sum;

   % truncamos la matriz de salinidad a los valores inferiores a S0
    sal_trunc = topaz_smos;
    %sal_trunc = sal;
    sal_trunc(sal_trunc > S0) = nan;

    % Para el c?lculo de fwc en la columna, debemos ponderar: Sz*dz
    for k = 1:length(dz)
        Sz_weighted(:,:,k) = (1-sal_trunc(:,:,k)/S0).*dz(k);
    end
    % integrate vertically (within the water column)
    FW = sum(Sz_weighted,3,'omitnan');
    FW(FW == 0) = nan; % avoiding zeros from 'omitnan' in the previous sum.

    % saving variables of interest
    fwc(:,:,it)=FW;
    sss(:,:,it)=sal(:,:,1);
    sst(:,:,it)=T(:,:,1);
    sss_smos_topaz(:,:,it)=topaz_smos(:,:,1);

    % cleaning
    clear sal T sal_trunc FW Sz_weighted
end

% Manage time variable -----------

if gaia==true
    tt = datetime(time,'ConvertFrom','datenum')';
    epoch_start = datenum('1970-01-01 00:00:00');
    tt_epoch = etime(datevec(tt),datevec(epoch_start));
else
    tt = datetime(time,'ConvertFrom','datenum')';
    tt_epoch = convertTo(tt,'epochtime','Epoch','1970-01-01');
end

% Saving a .mat file with selected variables
%cd ('/EVA/Documents/ICM/beaufort_FWC/outputs/topaz_original_monthly_fullArctic2011-2020/')
%save ('smos_topaz4b_monthly_fullArctic_fwc2011_2019_16m.mat',"sst","sss","fwc","time","lon","lat")

% 5. Create a netCDF file -------------------------------------------------

% create netcdf file
ncid = netcdf.create([path_save,'smos_topaz4b_monthly_fullArctic_fwc2011_2019_40m.nc'],'NC_CLOBBER');

% define dimensions
latdimID = netcdf.defDim(ncid,'lat',length(lat));
londimID = netcdf.defDim(ncid,'lon',length(lon));
ttdimID  = netcdf.defDim(ncid,'time',length(time));

% define variables in the new file
latid = netcdf.defVar(ncid,'latitude','double',latdimID);
lonid = netcdf.defVar(ncid,'longitude','double',londimID);
ttid  = netcdf.defVar(ncid,'time','double',ttdimID);

ccid  = netcdf.defVar(ncid,'sss','double',[londimID latdimID ttdimID]);
ddid  = netcdf.defVar(ncid,'sst','double',[londimID latdimID ttdimID]);
fwcid = netcdf.defVar(ncid,'fwc','double',[londimID latdimID ttdimID]);
ssid  = netcdf.defVar(ncid,'sss_smos_topaz','double',[londimID latdimID ttdimID]);

% create variable attributes
% units

netcdf.putAtt(ncid,latid,'valid_min','-90');
netcdf.putAtt(ncid,latid,'valid_max','90');
netcdf.putAtt(ncid,latid,'axis','Y');
netcdf.putAtt(ncid,latid,'units','degrees_north');
netcdf.putAtt(ncid,latid,'long_name','latitude');
netcdf.putAtt(ncid,latid,'standard_name','latitude');

netcdf.putAtt(ncid,lonid,'valid_min','-180');
netcdf.putAtt(ncid,lonid,'valid_max','180');
netcdf.putAtt(ncid,lonid,'axis','X');
netcdf.putAtt(ncid,lonid,'units','degrees_east');
netcdf.putAtt(ncid,lonid,'long_name','longitude');
netcdf.putAtt(ncid,lonid,'standard_name','longitude');

netcdf.putAtt(ncid,ttid,'units','seconds since 1970-01-01');
netcdf.putAtt(ncid,ttid,'long_name','time');
netcdf.putAtt(ncid,ttid,'standard_name','time');
netcdf.putAtt(ncid,ttid,'axis','T');
netcdf.putAtt(ncid,ttid,'calendar','gregorian');

netcdf.putAtt(ncid,ccid,'units','psu');
netcdf.putAtt(ncid,ddid,'units','Celsius degrees');
netcdf.putAtt(ncid,fwcid,'units','meters');
netcdf.putAtt(ncid,ssid,'units','psu');

% description
netcdf.putAtt(ncid,fwcid,'description','meters of freshwater within the water column (up to 34.8 psu or sea bottom');

% leave define mode and enter data mode to write data
netcdf.endDef(ncid)

% write data to variable.
netcdf.putVar(ncid,latid,double(lat));
netcdf.putVar(ncid,lonid,double(lon));
netcdf.putVar(ncid,ttid,double(tt_epoch));
netcdf.putVar(ncid,ccid,sss);
netcdf.putVar(ncid,ddid,sst);
netcdf.putVar(ncid,fwcid,fwc);
netcdf.putVar(ncid,ssid,sss_smos_topaz);

% close file to save changes
netcdf.close(ncid)


