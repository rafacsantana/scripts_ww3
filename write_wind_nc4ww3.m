clear
close all
run('/scale_wlg_persistent/filesets/home/santanarc/scripts/niwa/matlab/startup.m')

% Daily time
time_lima=datenum(1985,1,1,0,0,0):1:datenum(1990,1,1,0,0,0); % Graham Harrington ECAN request on 18/09/2023
%time_lima=datenum(1990,1,1,0,0,0):1:datenum(2014,12,31,0,0,0); % Graham Harrington ECAN request on 18/09/2023
%time_lima=datenum(1992,2,29,0,0,0):1:datenum(1992,2,29,0,0,0); % Graham Harrington ECAN request on 18/09/2023
time_lima=datenum(2014,12,16,0,0,0):1:datenum(2014,12,31,0,0,0); % Graham Harrington ECAN request on 18/09/2023
time_lima=datenum(2020,12,1,0,0,0):1:datenum(2022,01,31,0,0,0); % Graham Harrington ECAN request on 18/09/2023

%/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/projections/historical_proc/6hrly_global/wind/2014/wind_historical_GFDL-ESM4_CCAM_6hrly_Global_raw_2014121700.nc
%Index exceeds the number of array elements (0).
%
%Error in write_wind_nc4ww3 (line 60)
%  ua=squeeze((ncread(fname,'ua', [1 1 1 it(1)],[Inf Inf Inf it(end)-it(1)+1])));

% first merge u and v files
% cdo merge ua_historical_GFDL-ESM4_CCAM_6hrly_Global_raw.nc va_historical_GFDL-ESM4_CCAM_6hrly_Global_raw.nc wind_historical_GFDL-ESM4_CCAM_6hrly_Global_raw.nc

path_in='/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/projections/ssp370/6hrly_global/';
path_ou='/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/projections/historical_proc/6hrly_global/';

vars=[2];

vnames   ={'ua/','sic/'}; % (ua for wind) always with /
prefixs  ={'ua_ssp370_','sic_ssp370_'};
midfixs  ={'GFDL-ESM4_CCAM_' ,'GFDL-ESM4_CCAM_'};
sufixs   ={'6hrly_Global_raw','6hrly_Global_raw'};

scrsz=[1 1 1366 768];
scrsz=[2 42 958 953];
%scrsz=get(0,'screensize');

for vr=vars

  vname ={vnames{vr}};
  prefix={prefixs{vr}};
  midfix={midfixs{vr}};
  sufix ={sufixs{vr}};
  
  fname=[path_in,vname{1},prefix{1},midfix{1},sufix{1},'.nc'];
  display(['Loading: ',fname]);
  
  
  % Loop in model time
  for t=time_lima
    vn=vname{1};
  
    ptime=datestr(t,'YYYY/');
    ftime=datestr(t,'YYYYmmDDHH');
  
    if t==time_lima(1)
      display(['Reading lon, lat, time']); tic;
      FillValue=ncreadatt(fname,vn(1:end-1),'_FillValue'); 
      lon=(ncread(fname,'lon'));
      lat=(ncread(fname,'lat'));
      timew=squeeze(double(ncread(fname,'time'))); 
      if strcmp(prefix{1},'ua_historical_') || strcmp(prefix{1},'sic_historical_')
        time_ini=datenum(1959,1,1); 
      else
        time_ini=datenum(2015,1,1); 
      end
      time_mat=timew./(24*60)+time_ini; 
      timed=datevec(time_mat); 
      toc
      % climate models dont have leap years. removing 29/2 and adding a day to get to 31/21 of 2015/2100
      % ww3 will read 29/02 as 28/02
      timeww=timew;
      time_matt=timeww./(24*60)+time_ini; 
      timedd=datevec(time_matt); 
      for i=1:length(timeww)
        tt=time_matt(i);
        td=datevec(tt); 
        if td(:,2)==2 & td(:,3)==29
          it=find(timedd(:,1)==td(:,1) & timedd(:,2)==td(:,2) & timedd(:,3)==td(:,3));
          timeww(it:end)=timeww(it:end)+24*60; % plus a day (minutes)
          time_matt=timeww./(24*60)+time_ini; 
          timedd=datevec(time_matt); 
        end
      end 
      time_matt=timeww./(24*60)+time_ini; 
      timedd=datevec(time_matt); 
      timed=timedd;
      time_mat=time_matt;
      timew=timeww;
    end
  
    % find indexes of selected day on nc file
    td=datevec(t);
    %it=find(timed(:,1)==td(:,1) & timed(:,2)==td(:,2) & timed(:,3)==td(:,3));
    [dif it]=nanmin(abs(time_mat-t));
    if strcmp(sufix{1},'6hrly_Global_raw')
      it=it:it+4;
    else % hourly
      it=it:it+24;
    end
  
    % time that will be written
    timewit=(t:1/(length(it)-1):t+1)-time_ini;
    timewit=timewit(1:end)*(24*60);

    %display(['Processing times: ']) % ,num2str(timew(it)./(24*60)+time_ini)])
    %datevec(timewit./(24*60)+time_ini)
  
 
    % need to create data for 2015/01/01 00:00
    %if td(:,1)==2014 & td(:,2)==12 & td(:,3)==31
    %end
  
    if strcmp(prefix{1},'ua_historical_') || strcmp(prefix{1},'ua_ssp370_')
  
      display(['Reading winds']); tic;
      fname=[path_in,vname{1},prefix{1},midfix{1},sufix{1},'.nc'];
      %fname=[path_in,'ua/','ua_historical_',midfix{1},sufix{1},'.nc'];
      ua=squeeze((ncread(fname,'ua', [1 1 1 it(1)],[Inf Inf Inf it(end)-it(1)+1]))); 

      if strcmp(prefix{1},'ua_historical_') 
        fname=[path_in,'va/','va_historical_',midfix{1},sufix{1},'.nc'];
      elseif strcmp(prefix{1},'ua_ssp370_')
        fname=[path_in,'va/','va_ssp370_',midfix{1},sufix{1},'.nc'];
      end

      va=squeeze((ncread(fname,'va', [1 1 1 it(1)],[Inf Inf Inf it(end)-it(1)+1]))); 
      
      path_out=[path_ou,'wind/',ptime];
      system(['mkdir -p ',path_out]); 
      fda=[path_out,'wind_historical_',midfix{1},sufix{1},'_',ftime,'.nc'];
      display(['Saving: ',fda]);
  
      if exist(fda)==2
         system(['rm -rf ',fda]);
      end
  
      % creating nc and variables
      ncid = netcdf.create(fda,'NETCDF4');
      % time
      dtime =  netcdf.defDim(ncid,'time',length(it));
      varid = netcdf.defVar(ncid,'time','double',dtime);
      netcdf.putVar(ncid,varid,timewit)
      netcdf.putAtt(ncid,varid,'standard_name',"time");
      netcdf.putAtt(ncid,varid,'calendar',"gregorian");
      if strcmp(prefix{1},'ua_historical_') 
        netcdf.putAtt(ncid,varid,'units',"minutes since 1959-01-01 00:00:00");
      elseif strcmp(prefix{1},'ua_ssp370_')
        netcdf.putAtt(ncid,varid,'units',"minutes since 2015-01-01 00:00:00");
      end
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
  
    elseif strcmp(prefix{1},'sic_historical_') || strcmp(prefix{1},'sic_ssp370_')
  
      display(['Reading sic']); tic;
      sic=squeeze((ncread(fname,'sic', [1 1 it(1)],[Inf Inf it(end)-it(1)+1]))); 
      sic=sic./100;
  
      path_out=[path_ou,vname{1},ptime];
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
      netcdf.putVar(ncid,varid,timewit)
      netcdf.putAtt(ncid,varid,'standard_name',"time");
      netcdf.putAtt(ncid,varid,'calendar',"gregorian");
      if strcmp(prefix{1},'sic_historical_') 
        netcdf.putAtt(ncid,varid,'units',"minutes since 1959-01-01 00:00:00");
      elseif strcmp(prefix{1},'sic_ssp370_')
        netcdf.putAtt(ncid,varid,'units',"minutes since 2015-01-01 00:00:00");
      end
      % lon
      dlon =  netcdf.defDim(ncid,'lon',length(lon));
      varid = netcdf.defVar(ncid,'lon','double',dlon);
      netcdf.putVar(ncid,varid,lon)
      % lat
      dlat =  netcdf.defDim(ncid,'lat',length(lat));
      varid = netcdf.defVar(ncid,'lat','double',dlat);
      netcdf.putVar(ncid,varid,lat)
      % sic 
      sicvar = netcdf.defVar(ncid,'sic','float', [dlon dlat dtime]);
      netcdf.defVarFill(ncid,sicvar,false,FillValue);
      netcdf.putVar(ncid,sicvar,sic)
      netcdf.putAtt(ncid,sicvar,'standard_name',"sea_ice_area_fraction");
      netcdf.putAtt(ncid,sicvar,'units',"1");
  
      netcdf.close(ncid);
      
    end % if strcmp(prefix
  
  end

end
