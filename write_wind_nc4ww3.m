clear
close all
run('/scale_wlg_persistent/filesets/home/santanarc/scripts/niwa/matlab/startup.m')

% Daily time
time_lima=datenum(1990,1,4,0,0,0):1:datenum(2014,12,31,0,0,0); % Graham Harrington ECAN request on 18/09/2023

% first merge u and v files
% cdo merge ua_historical_GFDL-ESM4_CCAM_6hrly_Global_raw.nc va_historical_GFDL-ESM4_CCAM_6hrly_Global_raw.nc wind_historical_GFDL-ESM4_CCAM_6hrly_Global_raw.nc

check_file=1;
load_atm  =0;
set_xlim  =0;
proc_obs  =0;
plot_fig  =1;
save_fig  =0;
save_nc   =0;
save_csv  =1;

path_fig=['/scale_wlg_nobackup/filesets/nobackup/niwa03150/santanarc/figures/']; % GLOBALWAVE/'];%2018/01/05/00'];

path_efs='/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/projections/historical_proc/6hrly_global/wind/';

prefix  ={'wind_historical_',''};
midfix  ={'GFDL-ESM4_CCAM_' ,''};
sufix   ={'6hrly_Global_raw',''};

scrsz=[1 1 1366 768];
scrsz=[2 42 958 953];
%scrsz=get(0,'screensize');

k=30;

fname=[path_efs,prefix{1},midfix{1},sufix{1},'.nc'];
display(['Loading: ',fname]);

% Loop in model time
for t=time_lima

  ptime=datestr(t,'YYYY/');
  ftime=datestr(t,'YYYYmmDDHH');

  if t==time_lima(1)
    display(['Reading lon, lat, time']); tic;
    FillValue=ncreadatt(fname,'ua','_FillValue'); 
    lon=(ncread(fname,'lon'));
    lat=(ncread(fname,'lat'));
    timew=squeeze(double(ncread(fname,'time'))); 
    time_mat=timew./(24*60)+datenum(1959,1,1); 
    timed=datevec(time_mat); 
    toc
  end

  % find indexes of selected day on nc file
  td=datevec(t);
  it=find(timed(:,1)==td(:,1) & timed(:,2)==td(:,2) & timed(:,3)==td(:,3));


  display(['Reading winds']); tic;
  ua=squeeze((ncread(fname,'ua', [1 1 1 it(1)],[Inf Inf Inf it(end)-it(1)+1]))); 
  va=squeeze((ncread(fname,'ua', [1 1 1 it(1)],[Inf Inf Inf it(end)-it(1)+1]))); 


  save_nc=1;
  if save_nc==1

    path_out=[path_efs,ptime];
    system(['mkdir -p ',path_out]); 
    fda=[path_out,prefix{1},midfix{1},sufix{1},'_',ftime,'.nc'];
    display(['Saving: ',fda]);

    if exist(fda)==2
       system(['rm -rf ',fda]);
    end

    % creating nc and variables
    ncid = netcdf.create(fda,'NETCDF4');
    % time
    dtime =  netcdf.defDim(ncid,'time',length(it));
    varid = netcdf.defVar(ncid,'time','double',dtime);
    netcdf.putVar(ncid,varid,timew(it))
    netcdf.putAtt(ncid,varid,'standard_name',"time");
    netcdf.putAtt(ncid,varid,'calendar',"gregorian");
    netcdf.putAtt(ncid,varid,'units',"minutes since 1959-01-01 00:00:00");
    % lon
    dlon =  netcdf.defDim(ncid,'lon',length(lon));
    varid = netcdf.defVar(ncid,'lon','double',dlon);
    netcdf.putVar(ncid,varid,lon)
    % lat
    dlat =  netcdf.defDim(ncid,'lat',length(lat));
    varid = netcdf.defVar(ncid,'lat','double',dlat);
    netcdf.putVar(ncid,varid,lat)
    % ua 
    uvar = netcdf.defVar(ncid,'ua','float', [dlon dlat dtime]);
    netcdf.defVarFill(ncid,uvar,false,FillValue);
    netcdf.putVar(ncid,uvar,ua)
    netcdf.putAtt(ncid,uvar,'standard_name',"eastward_wind");
    netcdf.putAtt(ncid,uvar,'units',"m s-1");
    % va 
    vvar = netcdf.defVar(ncid,'va','float', [dlon dlat dtime]);
    netcdf.defVarFill(ncid,vvar,false,FillValue);
    netcdf.putVar(ncid,vvar,va)
    netcdf.putAtt(ncid,vvar,'standard_name',"northward_wind");
    netcdf.putAtt(ncid,vvar,'units',"m s-1");

    netcdf.close(ncid);

    %ncdisp(fda)
%return
%
%
%    nccreate(fda,'lon','Dimensions' , {'lon',  length(lon)});
%    nccreate(fda,'lat','Dimensions' , {'lat',  length(lat)});
%    nccreate(fda,'time','Dimensions', {'time', length(it)});
%    nccreate(fda,'ua','Dimensions'  , {'lon',  length(lon), 'lat', length(lat), 'time',  length(it)}); 
%    nccreate(fda,'va','Dimensions'  , {'lon',  length(lon), 'lat', length(lat), 'time',  length(it)}); 
%
%
%    ncwriteatt(fda,'time','standard_name',"time");
%    ncwriteatt(fda,'time','calendar',"gregorian");
%    ncwriteatt(fda,'time','units',"minutes since 1959-01-01 00:00:00");
%    ncwriteatt(fda,'time','axis',"T");
%
%    ncwriteatt(fda,'ua','standard_name',"eastward_wind");
%    ncwriteatt(fda,'ua','units',"m s-1");
%    %ncwriteatt(fda,'ua','_FillValue', "9.96921e+36f");
%    %ncwriteatt(fda,'ua','missing_value'," 9.96921e+36f");
%
%    ncwriteatt(fda,'va','standard_name',"northward_wind");
%    ncwriteatt(fda,'va','units',"m s-1");
%    %ncwriteatt(fda,'va','_FillValue',"9.96921e+36f");
%    %ncwriteatt(fda,'va','missing_value',"9.96921e+36f");
%
%    ncwrite(fda,'time',timew(it));
%    ncwrite(fda,'ua',ua)
%    ncwrite(fda,'va',ua)
%
%    netcdf.defVarFill(fda,false,FillValue);
%    netcdf.defVarFill(fda,false,FillValue);

  end

end

