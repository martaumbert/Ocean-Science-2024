%%%%%% SCRIPT TO COMPUTE FWC from TOPAZ4B OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% 1. List TOPAZ files
% 2. Load common variables once (lat, lon and depth)
% 3. Calculate FWC (WARNING! dz = dz(z))
% 4. Create netcdf file
% -------------------------------------------------------------------------
% Author    : Maria Sanchez from Eva De Andrés script
% Email     : eva.deandres@upm.es
% Creation  : January 2023
% Modified  : Marta Umbert April 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;       % clear Command Window
clear;     % clear Workspace
close all  % Close all figures

% 1. List TOPAZ files -----------------------------------------------------

% Path
addpath('/mnt/lustre/users/BECPolar/data/regrid_cdo/TOPAZ4/');
pathfile = '/mnt/lustre/users/BECPolar/data/regrid_cdo/TOPAZ4/';

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

% 3. Calculate FWC (WARNING! dz = dz(z)) ----------------------------------
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

% Estimate FWC for all files on list
for it=1:numel(C{1})
    filename = cell2mat(C{1}(it));
    filedate = datenum(filename(1:8),'yyyymmdd');
    time(it) = filedate;
    sal      = ncread(filename,'so');
    T        = ncread(filename,'thetao');

    % truncamos la matriz de salinidad a los valores inferiores a S0
    sal_trunc = sal;
    sal_trunc(sal_trunc > S0) =nan;
   
    % Para el cálculo de fwc en la columna, debemos ponderar: Sz*dz
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
    
    % cleaning
    clear sal T sal_trunc FW Sz_weighted
end

% convert time in 'seconds from 10970-01-01'
tt = datetime(time,'ConvertFrom','datenum')';
epoch_start = datenum('1970-01-01 00:00:00');
tt_epoch = etime(datevec(tt),datevec(epoch_start));

%tt_epoch = convertTo(tt,'epochtime','Epoch','1970-01-01'); % no funciona en gaia

% Saving a .mat file with selected variables
%save ('/mnt/lustre/users/BECPolar/analysis/fwc/topaz4b_monthly_fullArctic_fwc2011_2019_latlon.mat',"sst","sss","fwc","time","lon","lat")

% 5. Create a netCDF file -------------------------------------------------

% path settings
pathfile='/mnt/lustre/users/BECPolar/analysis/fwc/';

% create netcdf file
ncid = netcdf.create([pathfile,'topaz4b_monthly_fullArctic_fwc2011_2019_latlon.nc'],'NC_CLOBBER');

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

% create variable attributes
% units
netcdf.putAtt(ncid,latid,'units','degree_north');
netcdf.putAtt(ncid,latid,'axis','Y');
netcdf.putAtt(ncid,lonid,'units','degree_east');
netcdf.putAtt(ncid,lonid,'axis','X');
netcdf.putAtt(ncid,ttid,'units','seconds since 1970-01-01');
netcdf.putAtt(ncid,ttid,'standard_name','time');
netcdf.putAtt(ncid,ttid,'axis','T');
netcdf.putAtt(ncid,ttid,'calendar','gregorian');

netcdf.putAtt(ncid,ccid,'units','psu');
netcdf.putAtt(ncid,ddid,'units','Celsius degrees');
netcdf.putAtt(ncid,fwcid,'units','meters');

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

% close file to save changes
netcdf.close(ncid)
