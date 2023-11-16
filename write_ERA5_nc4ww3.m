clear
close all
%run('/scale_wlg_persistent/filesets/home/santanarc/scripts/niwa/matlab/startup.m')

% Daily time
time_lima=datenum(1985,4,1,0,0,0):1:datenum(1987,12,31,0,0,0); % Gr
%time_lima=datenum(1990,1,1,0,0,0):1:datenum(2014,12,31,0,0,0); %
%time_lima=datenum(1992,2,29,0,0,0):1:datenum(1992,2,29,0,0,0); %
%time_lima=datenum(2014,01,01,0,0,0):1:datenum(2014,12,31,0,0,0);
%time_lima=datenum(2015,01,01,0,0,0):1:datenum(2015,12,31,0,0,0);
%time_lima=datenum(2020,4,1,0,0,0):1:datenum(2020,12,1,0,0,0); 
%time_lima=datenum(2020,12,1,0,0,0):1:datenum(2020,12,31,0,0,0); 
%time_lima=datenum(2021,1,1,0,0,0):1:datenum(2021,12,31,0,0,0); 

path_ou='/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/ERA5/';
path_in=path_ou; 

vars=[1,2];

% ERA5_siconc_y1985.nc
% ERA5_u10_y1985.nc
% ERA5_v10_y1985.nc

prefixs  ={'ERA5_'}; % prior to 2015/01/01
vnames   ={'u10_y','siconc_y'}; % (ua for wind) always with /

scrsz=[1 1 1366 768];
scrsz=[2 42 958 953];
%scrsz=get(0,'screensize');

fold=''; % temporary old filename 
for vr=vars

  vname ={vnames{vr}};
  prefix={prefixs{1}};
  
  % Loop in model time
  for t=time_lima
    vn=vname{1};

    ptime=datestr(t,'YYYY');
    ftime=datestr(t,'YYYYmmDDHH');

    fname=[path_in,prefix{1},vname{1},ptime,'.nc'];
    display(['Loading: ',fname]);
 
    if strcmp(fname,fold)~=1
      display(['Reading lon, lat, time']); tic;
      FillValue=ncreadatt(fname,vn(1:end-2),'_FillValue'); 
      lon=(ncread(fname,'longitude'));
      lat=(ncread(fname,'latitude'));
      % flipping latitude so it starts at -90
      lat=flipud(lat);
      timew=squeeze(double(ncread(fname,'time'))); 

      time_ini=datenum(1900,1,1); 
      time_mat=(timew./(24))+time_ini; 
      timed=datevec(time_mat); 
      toc
    end
  
    % find indexes of the selected day on nc file
    td=datevec(t);
    it=find(timed(:,1)==td(:,1) & timed(:,2)==td(:,2) & timed(:,3)==td(:,3));

    if it(end)~=length(timew)
      it(end+1)=it(end)+1;
    end

    % time variable that will be written in the nc file
    timewit=timew(it);

    if strcmp(vname{1},'u10_y') 
 
      display(['Reading winds']); tic;
      fname=[path_in,prefix{1},vname{1},ptime,'.nc'];
      ua=squeeze((ncread(fname,'u10', [1 1 it(1)],[Inf Inf it(end)-it(1)+1]))); 

      fname=[path_in,prefix{1},'v10_y',ptime,'.nc'];
      va=squeeze((ncread(fname,'v10', [1 1 it(1)],[Inf Inf it(end)-it(1)+1]))); 
     
      % flipping u and v according to latitude 
      ua=fliplr(ua);
      va=fliplr(va);
 
      path_out=[path_ou,ptime,'/'];
      system(['mkdir -p ',path_out]); 
      fda=[path_out,'winds_10m_era5_',ftime,'.nc'];
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
      netcdf.putAtt(ncid,varid,'units',"hours since 1900-01-01 00:00:00.0");

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
  
    elseif strcmp(vname{1},'siconc_y') 
  
      display(['Reading sic']); tic;
      sic=squeeze((ncread(fname,'siconc', [1 1 it(1)],[Inf Inf it(end)-it(1)+1]))); 

      % flipping u and v according to latitude 
      sic=fliplr(sic);
  
      path_out=[path_ou,ptime,'/'];
      system(['mkdir -p ',path_out]); 
      fda=[path_out,'sic_era5_',ftime,'.nc'];
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
      netcdf.putAtt(ncid,varid,'units',"hours since 1900-01-01 00:00:00.0");
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
