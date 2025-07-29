clear
close_all=1;
if close_all==1; close all; end
warning off
tic
  
% Daily time
time_lima=datenum(1984,1,1,0,0,0):1:datenum(1993,12,31,0,0,0);    % 1-year 
time_lima=datenum(1994,1,1,0,0,0):1:datenum(2003,12,31,0,0,0);    % 1-year 
time_lima=datenum(2004,1,1,0,0,0):1:datenum(2013,12,31,0,0,0);    % 1-year 
time_lima=datenum(2014,1,1,0,0,0):1:datenum(2023,12,31,0,0,0);    % 1-year 
time_lima=datenum(1984,1,1,0,0,0):1:datenum(2023,12,31,0,0,0);    % 40 years 


time_lima=datenum(2020,12,2,0,0,0):1:datenum(2020,12,26,0,0,0);    % 1-year
time_lima=datenum(2021,1,1,0,0,0):1:datenum(2021,12,31,0,0,0);    % 1-year
%time_lima=datenum(2018,1,31,0,0,0):1:datenum(2018,2,1,0,0,0);    % 1-year 
%
%time_lima=datenum(2018,1,6,0,0,0):1:datenum(2024,8,20,0,0,0);    % NZWAVE-HR 
%time_lima=datenum(2019,1,1,0,0,0):1:datenum(2019,12,31,0,0,0);    % Graham Harrington 
%time_lima=datenum(2020,1,1,0,0,0):1:datenum(2023,12,31,0,0,0);    % Graham Harrington 
%time_lima=datenum(2018,7,1,0,0,0):1:datenum(2018,7,22,0,0,0);    % Graham Harrington 
%time_lima=datenum(2018,1,26,0,0,0):1:datenum(2018,2,5,0,0,0);    % Graham Harrington 
%time_lima=datenum(2018,7,23,0,0,0):1:datenum(2018,12,31,0,0,0);    % Graham Harrington 
%time_lima=datenum(2018,7,22,0,0,0):1:datenum(2018,7,24,0,0,0);    % Graham Harrington 

%time_lima=datenum(2018,1,6,0,0,0):1:datenum(2024,8,20,0,0,0); % Graham Harrington ECAN request on 18/09/2023
%time_lima=datenum(1992,1,1,0,0,0):1:datenum(1996,12,31,0,0,0); % Graham Harrington ECAN request on 18/09/2023
%time_lima=datenum(1997,1,1,0,0,0):1:datenum(1999,12,31,0,0,0); % Graham Harrington ECAN request on 18/09/2023
%time_lima=datenum(2000,1,1,0,0,0):1:datenum(2005,12,31,0,0,0);    % 1-year 
%time_lima=datenum(1992,1,1,0,0,0):1:datenum(2010,12,31,0,0,0);    % 1-year 

%
%time_lima=datenum(2014,7,9,0,0,0):1:datenum(2014,9,3,0,0,0);     % initial date NZLAM hindcast

%time_lima=datenum(2018,1,6,0,0,0):1:datenum(2024,8,20,0,0,0); % NZWAVE-HR
time_lima=datenum(2005,1,1,0,0,0):1:datenum(2014,12,30,0,0,0);   % NZWAVE-GFDL-CCAM

%time_lima=datenum(2030,1,1,0,0,0):1:datenum(2039,12,31,0,0,0); % climate projections
%time_lima=datenum(2050,1,1,0,0,0):1:datenum(2069,12,31,0,0,0); % climate projections
%time_lima=datenum(2080,1,1,0,0,0):1:datenum(2099,12,31,0,0,0); % climate projections

files={'Banks_Peninsula','Baring_Head'            ,'wairewa_lake_forsyth','Pelorus_Sound'  ,'Taumutu',...
       'Mangawhai'      ,'BoP_leatherback_turtles','Hjort_Trench'        ,'Puysegur_Trench','Port_Waikato',...
       'Blenheim'}; 

% line below cannot be commented (variable used down there even if not intended to be used)
scase='graveyard_hills'; %dixon-anderson'; % cant be commented (variable used down there even if not intended to be used)
%[tstation,lat_obss,lon_obss]=read_wave_stations(scase); % when no obs are available. it's used in data request for a group of stations
%files=tstation;
%stations=1:length(files);

%stations=2:length(files);
stations=[1]; % 0 = coastal stations (points)

runs=[2,15,16,14];
runs=[9];
%runs=[2:5];
expt_names={'GLOBALWAVE'          ,'NZWAVE'             ,'NZWAVE-ST6'         ,'NZWAVE-HR-NOTIDES','NZWAVE-HR',... 5
            'GLOBALWAVE-GFDL-CCAM','NZWAVE-ST4-GLOBALUM','NZWAVE-ST6-GLOBALUM','NZWAVE-GFDL-CCAM' ,'NZWAVE-ERA5',...10
            'GLOBALWAVE-ERA5'     ,'NZWAVE-ERA5-2021'   ,'NZWAVE-ST6-NZLAM'   ,'NZWAVE-ST6-BOM'   ,'NZWAVE-ST6-DEFAULT',... 15
            'NZWAVE-ST6-BOM-WCOR' ,}; % 20            
% switches
plot_series=1; % stations=0, plot_series=1, plot_coastm=1, plot_obs==0 to save the coastal points
plot_rose  =0;
plot_disp  =0;
plot_eva   =0;

plot_coastm=0; % statistical maps for coastal points
plot_map   =0; % daily station maps
plot_matm  =0; % maps with atmospheric forcing
plot_downs =0; % daily downscaled maps
plot_mmap  =0; % mean model maps
plot_mserie=0; % mean model series
plot_amap  =0; % altimeter maps
plot_aserie=0; % altimeter series

proc_obs   =0;
ck_mod_mat =0;
display_txt=1;
plot_obs   =1; % stations=0, plot_series=1, plot_coastm=1, plot_obs==0 to save the coastal points
plot_mod   =1; % model data
plot_harm  =0; % harmonic analysis
plot_atm   =1; % plot atmospheric forcing
save_csv   =0;
save_mat   =1;
save_fig   =0;
save_video =0;
save_gif   =0;

set_xlim   =1;
obs_first  =0;
add_stats  =1; % to legend
portrait   =1;

interp_half_hour=1;
correct_model=2; % -1=scatter_kde, 0=binscatter, 1=nothing, >1 correct it
select_cpoint=0; % select one coastal point based on lon_c lat_c variables given below
time_limab=time_lima;

% series, rose, dispersion, eva

% variables to be plotted in the time series
%NEVER CHANGE THIS ORDER
% 1 {'PeakPerid'}   2 {'PeakDirection'}   3 {'Directional spread'}   4 {'Tm01'}   5 {'Tm02'}   6 {'Hs'}   7 {'Qp'}  8 Tm
% 9 wlv 10 ucurr 11 vcurr
dcol=[1,4,5]; % select variables from model data
%dcol=[6,5,8,9]; %[6,1,2]; % Tp, Tm01, Tm02
dcol=[6,2,1]; % Hs, wlv
dcol=[6];

% maps of snapshots
pcname ='hs'; % m_pcolor: wlv, hs, curr
vecname='hs'; % m_vec: hs, curr
conname='psl'; % m_contour: psl

% mmaps = model mean maps
mcol=[1,3]; % selected variables from mvalm
mvalm={'hs_m_mean','tp_m_mean','we_m_mean'}; 
mvals={'hs_mean','tp_mean','we_mean'}; % 1=hs, 2=tp, 3=we

% altimeter stats
acol=[2]; % selected variables from avalm, avals, ovalm
%acol=[7:12];
%series
crop_l4=0;
avals={'hs_mean_rmse','hs_maxi_rmse','hs_mean_bias','hs_maxi_bias','hs_mean_corr','hs_maxi_corr',...
       'hs_mean_frmse','hs_maxi_frmse','hs_mean_fbias','hs_maxi_fbias','hs_mean_fcorr','hs_maxi_fcorr'};
%maps
avalm={'hs_mean_mrmse','hs_maxi_mrmse','hs_mean_mbias','hs_maxi_mbias','hs_mean_mcorr','hs_maxi_mcorr',... 
       'hs_mean_fmrmse','hs_maxi_fmrmse','hs_mean_fmbias','hs_maxi_fmbias','hs_mean_fmcorr','hs_maxi_fmcorr'}; % 
% obs maps
ovalm={'hs_mean_ostd','hs_maxi_ostd','hs_mean_omean','hs_maxi_omean','hs_mean_omax','hs_maxi_omax',... 
       'hs_mean_ostd','hs_maxi_ostd','hs_mean_omean','hs_maxi_omean','hs_mean_omax','hs_maxi_omax'}; % 


if length(runs)==1
  colors={'b'};
else
	colors={'c','r','m','b','k','c'};
end

for i=1:length(runs)
  expts{i}=expt_names{runs(i)};
end

%sufixes in netcdf files
gnames={'globalwave+globalum','nzwave+nzlam','nzwave_hr+nzcsm','tongawave+globalum','nzwave_era5'};

path_source=['/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/']; % GLOBALWAVE/'];%2018/01/05/00'];
%path_source=['~/EFS/']; 
path_santanarc='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/';
path_obs=[path_santanarc,'data/obs/'];
path_fig='/scale_wlg_nobackup/filesets/nobackup/niwa03150/santanarc/figures/';

% stations or grids - stations are bad in Baring_Head
ww3pre={'ww3p_interp','ww3g'};
ww3pre=ww3pre{2};

%fname=[pname,ww3pre,'_',ftime,'-utc_',gnames{ig},'.nc'];
fname='/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/NZWAVE/2021/01/01/00/ww3g_2021010100-utc_nzwave+nzlam.nc';
depth_mod=double(ncread(fname,'depth',[1 1 1],[Inf Inf 1]));
lon_mod=double(ncread(fname,'lon'));
lat_mod=double(ncread(fname,'lat'));

fname='/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/NZWAVE-HR/2021/01/01/00/ww3g_2021010100-utc_nzwave_hr+nzcsm.nc';
depth_hr=double(ncread(fname,'depth',[1 1 1],[Inf Inf 1]));
lon_hr=double(ncread(fname,'lon'));
lat_hr=double(ncread(fname,'lat'));

% Loading satellite obs
if plot_aserie || plot_amap
  atype='cmems_nrt';
  % cmems_l4  
  if strcmp(atype,'cmems_l4')
    path_cop='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/data/copernicus-l4-wave-data/';
    file_cop=[path_cop,'cmems_obs-wave_glo_phy-swh_my_multi-l4-2deg_P1D_multi-vars_143.00E-185.00E_53.00S-21.00S_2002-01-15-2020-12-31.nc'];
    time_cop=ncread(file_cop,'time')+datenum('1950-01-01');
    lon_cop=ncread(file_cop,'longitude');
    lat_cop=ncread(file_cop,'latitude');
    hs_cop_mean=ncread(file_cop,'VAVH_DAILY_MEAN');
    hs_cop_maxi=ncread(file_cop,'VAVH_DAILY_MAX');

  elseif strcmp(atype,'cmems_nrt')
    path_cop='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/data/copernicus-2021-wave-data/';
    file_cop=[path_cop,'cmems_obs-wave_glo_phy-swh_nrt_multi-l4-2deg_P1D_multi-vars_143.00E-185.00E_53.00S-21.00S_2021-01-01-2021-12-31.nc'];
    %ncdisp(file_cop)
    %return
    time_cop=ncread(file_cop,'time')+datenum('1950-01-01');
    lon_cop=ncread(file_cop,'longitude');
    lat_cop=ncread(file_cop,'latitude');
    hs_cop_mean=ncread(file_cop,'VAVH_DAILY_MEAN');
    hs_cop_maxi=ncread(file_cop,'VAVH_DAILY_MAX');
    % cropping l4 data for NZWAVE-HR area
  end

  [dif ilonc]=nanmin(abs(lon_cop-min(lon_hr))); [dif flonc]=nanmin(abs(lon_cop-max(lon_hr)));
  [dif ilatc]=nanmin(abs(lat_cop-min(lat_hr))); [dif flatc]=nanmin(abs(lat_cop-max(lat_hr)));
  if crop_l4
    lon_cop=lon_cop(ilonc:flonc);
    lat_cop=lat_cop(ilatc:flatc);
    hs_cop_mean=hs_cop_mean(ilonc:flonc,ilatc:flatc,:);
	  hs_cop_maxi=hs_cop_maxi(ilonc:flonc,ilatc:flatc,:);
  end

end

if plot_downs %|| plot_coastm
%if plot_map || plot_mmap || plot_downs || plot_coastm

  % gebco bathy
  gebco_nc='/scale_wlg_persistent/filesets/project/niwa00020/data/topography/gebco/1-min/gridone.nc';
  spac=double(ncread(gebco_nc,'spacing'));
  lon_gebco=double(ncread(gebco_nc,'x_range')); lon_gebco=lon_gebco(1):spac(1):lon_gebco(end);
  lat_gebco=double(ncread(gebco_nc,'y_range')); lat_gebco=lat_gebco(end):-spac(2):lat_gebco(1);
  bat_gebco=double(ncread(gebco_nc,'z'));
  bat_gebco=reshape(bat_gebco,length(lon_gebco),length(lat_gebco));
  lon_gebco=lon_gebco(1:10:end);
  lat_gebco=lat_gebco(1:10:end);
  bat_gebco=bat_gebco(1:10:end,1:10:end);
  lon_gebco=wrapTo360(lon_gebco);

end

if plot_matm

  path_ssp='/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/projections/historical_proc/6hrly_global/';
  file_ssp='psl_ssp370_GFDL-ESM4_CCAM_6hrly_Global_raw_2030033000.nc';
  file_ssp=[path_ssp,'psl/2030/',file_ssp];
  lon_ssp=ncread(file_ssp,'lon');
  lat_ssp=ncread(file_ssp,'lat');
  [lon_sspm,lat_sspm]=meshgrid(lon_ssp,lat_ssp);

  % reading track file
  path_tr='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/data/ccam_cyclone_tracks/maxes/';
  if time_lima(1)<datenum(2015,1,1)
    file_tr=[path_tr,'GFDL-ESM4_historical_nz_maxes_1985_2014'];
  else
    file_tr=[path_tr,'GFDL-ESM4_ssp370_nz_maxes_2070_2099'];
  end
  tr=importdata([file_tr,'.csv']);
  tr.textdata(1,:)
  for i=2:length(tr.textdata)
    time_tr(i-1)=datenum(tr.textdata{i,4},'YYYY-mm-ddTHH:MM:SS.000000000');
    id_tr(i-1)=str2num(tr.textdata{i,2});
  end
  
end

if plot_coastm || stations(1)==0

  % finding the closest points to the coast
  np=1; % number of the closest points 
  xl=[390,605]; yl=[180,520];
  xl=[390,700]; yl=[180,520]; % including the Chatham Islands
  dep=depth_mod(xl(1):xl(2),yl(1):yl(2));

  % finding y points in the x direction
  k=0;
  for i=1:size(dep,1)
    [ngap,igap,fgap,sizegap] = findgap(dep(i,:)');
    if ngap>0
      for j=1:ngap; k=k+1;
        points=[igap(j)-np:igap(j)-1,fgap(j)+1:fgap(j)+np];
        yp(k).p=points;
        xp(k).p(1:length(points))=i;
        end; end; end
  % finding x points in the y direction
  k=0;
  for i=1:size(dep,2)
    [ngap,igap,fgap,sizegap] = findgap(dep(:,i));
    if ngap>0
      for j=1:ngap; k=k+1;
        points=[igap(j)-np:igap(j)-1,fgap(j)+1:fgap(j)+np];
        xpp(k).p=points;
        ypp(k).p(1:length(points))=i;
        end; end; end
  % concatenating x and y points
  xpoints=[]; ypoints=[];
  for i=1:length(xp)
    xpoints=[xpoints,xp(i).p]; 
		ypoints=[ypoints,yp(i).p];
  end
  for i=1:length(xpp)
    xpoints=[xpoints,xpp(i).p]; 
		ypoints=[ypoints,ypp(i).p];
  end

  % cleaning repeated pair of points
	points=[xpoints',ypoints'];
  points=unique(points,'rows');
  points(:,1)=points(:,1)+xl(1)-1;
  points(:,2)=points(:,2)+yl(1)-1;

  % selecting from spectral points
  spec_points=0;
  if spec_points
    fname='/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/NZWAVE-ERA5/2021/01/01/00/ww3p_2021010100-utc_nzwave_era5.nc';
    lon_spe=double(ncread(fname,'lon'));
    lat_spe=double(ncread(fname,'lat'));
    depth_spe=double(ncread(fname,'depth'));

    a=find(lon_spe(:,1)<=178.7 & lon_spe(:,1)>=166.2 & lat_spe(:,1)<=-34 & lat_spe(:,1)>=-48);
    b=find(depth_spe(a,1)<=1000);
    c=a(b);
    a=find(lon_spe(c,1)<=170.6 & lon_spe(c,1)>=168 & lat_spe(c,1)<=-37.5 & lat_spe(c,1)>=-42.5);
    c(a)=nan; return
    c(isnan(c))=[];
    figure; hold on; plot(lon_spe(c,1),lat_spe(c,1),'.'); xlim([165 180]); ylim([-50 -30])
	end


end

lon_obss=[]; lat_obss=[];

% selected coastal point 
lon_c=174.697967; % port waikato
lat_c=-37.385506;
if select_cpoint
  [dis]=sqrt((lon_mod(points(:,1))-lon_c).^2 + (lat_mod(points(:,2))-lat_c).^2);
  [dif iloc]=nanmin(abs(dis));
end


if plot_coastm==1 & plot_series==0
	stations=1:length(points);
  if select_cpoint
   	stations=iloc;
  end
end


write_coastal_points_nc=0;
if write_coastal_points_nc==1;

  path_expt=[path_source,'NZWAVE-ERA5/']; 
  path_dm=[path_expt,'matlab/'];
  system(['mkdir -p ',path_dm]);
  filename=[path_dm,'nzwave_coastal_points_np_',num2str(np),'.nc'];
  display(['Saving: ',filename])
  if exist(filename)==2
   system(['rm -rf ',filename]);
  end
  
  lon_points=lon_mod(points(:,1)); lat_points=lat_mod(points(:,2));
  coastal_points=points;
  
  nccreate(filename,'lon_points','Dimensions',{'lon_points',length(lon_points)},'Datatype','double','Format','classic');
  ncwrite(filename,'lon_points',lon_points);
  nccreate(filename,'lat_points','Dimensions',{'lat_points',length(lat_points)},'Datatype','double','Format','classic');
  ncwrite(filename,'lat_points',lat_points);
  nccreate(filename,'coastal_points','Dimensions',{'lon_points',length(lon_points),'pair_ij',2},'Datatype','double','Format','classic');
  ncwrite(filename,'coastal_points',coastal_points);

end

% biggest loop on stations
for fobs=stations
    
  % OBS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if plot_coastm==1 & plot_series==0
    file=['i',num2str(points(fobs,1)),'_j',num2str(points(fobs,2)),'_grid_point_'];
  elseif plot_coastm==1 & plot_series==1
    file=['bla'];
    %file=['i',num2str(points(fobs,1)),'_j',num2str(points(fobs,2)),'_grid_point_'];
  else
    file=files{fobs};
    file=replace(file,' ','_');
    file_obs=[path_obs,file];
    [obs,time_obs,lon_obs,lat_obs]=proc_wave_obs(file,proc_obs,time_lima,scase); %,plot_obs);   
  end
  
  if plot_series
   
    if portrait
      scrsz=[2    42   958   953];
    else
      scrsz=[1 41 1920 962];
      %scrsz=[1 1 1366 768];
      %scrsz=get(0,'screensize');
    end
 
    if close_all==1; figure('position',scrsz,'color',[1 1 1],'visible','on'); end
    hold on
    set(gca,'fontsize',12,'fontweight','bold')
    
    if obs_first
      legh={'Obs'}; legt={'Obs'};
    else
      legh=[]; legt=[];
    end

    ke=0;
    for expt=expts
      ke=ke+1;
    
      expt=expt{1};
      [ig,ltime]=grab_nc_sufix(expt);
    
      path_expt=[path_source,expt,'/']; % GLOBALWAVE/'];%2018/01/05/00'];
      path_dm=[path_expt,'matlab/'];
      system(['mkdir -p ',path_dm]);
      %filename=[path_dm,file,'_',datestr(time_lima(1),'YYYYmmDDHH'),'_',datestr(time_lima(end),'YYYYmmDDHH'),'.mat'];
      hs_mod=[];tp_mod=[]; pd_mod=[];time_mod=[];
      tm01_mod=[];tm_mod=[]; md_mod=[]; we_mod=[];
      tm02_mod=[];ds_mod=[]; wlv_mod=[];time_mod=[];
      ucur_mod=[];vcur_mod=[]; 

      % checking if the whole time series is saved as matlab file
      filename=[path_dm,file,'_',datestr(time_lima(1),'YYYYmmDDHH'),'_',datestr(time_lima(end),'YYYYmmDDHH'),'.mat'];
      if ck_mod_mat==1 && exist(filename)==2

        disp(['Loading: ',filename]); % loading long timeseries rather than daily files
        load(filename); % loading long timeseries rather than daily files

      else

        checkin=0;  
        % Loop in model time
        for t=time_lima
          %tic
          if plot_coastm==1 & plot_series==0 % plot_series=0 needed to save the coastal points
            file=['i',num2str(points(fobs,1)),'_j',num2str(points(fobs,2)),'_grid_point_'];
          elseif plot_coastm==1 & plot_series==1
            file=['bla'];
            %file=['i',num2str(points(fobs,1)),'_j',num2str(points(fobs,2)),'_grid_point_'];
          elseif ke==1 & t==time_lima(1)
            %file=files{fobs};
            %file=replace(file,' ','_');
            %file_obs=[path_obs,file];
            %[obs,time_obs,lon_obs,lat_obs]=proc_wave_obs(file,proc_obs); %,plot_obs);   
          end

          filename=[path_dm,file,'_',datestr(t,'YYYYmmDDHH'),'.mat'];
    
          if ck_mod_mat==1 && exist(filename)==2
            if display_txt==1; display(['Loading: ',filename]); end
            load(filename)%,'time_mod','model')
            if exist('ucur','var')~=1
              wlv=nan(size(hs)); ucur=nan(size(hs)); vcur=nan(size(hs));
              we=nan(size(hs)); ucur=nan(size(hs)); vcur=nan(size(hs));
            end 
 
          elseif fobs==0 %plot_coastm==1, working on the coastal points

            ptime=datestr(t,'YYYY/mm/DD/HH/');
            ftime=datestr(t,'YYYYmmDDHH');
            pname=[path_expt,ptime];
            fname=[pname,ww3pre,'_',ftime,'-utc_',gnames{ig},'.nc'];
            display(['Loading: ',fname]);

            if checkin==0; % t==time_lima(1)
              checkin=1; 
              lon_mod=double(ncread(fname,'lon'));
              lat_mod=double(ncread(fname,'lat'));
            end

           	stations=1:length(points);
            for fobss=stations
              %tic; 

              file=['i',num2str(points(fobss,1)),'_j',num2str(points(fobss,2)),'_grid_point_'];
              filename=[path_dm,file,'_',datestr(t,'YYYYmmDDHH'),'.mat'];
              %display(['Checking: ',filename])

              if ck_mod_mat==0 || exist(filename)~=2

                lon_obs=lon_mod(points(fobss,1));
			  			  lat_obs=lat_mod(points(fobss,2));
                [dif ilon]=nanmin(abs(lon_mod-lon_obs));
                [dif ilat]=nanmin(abs(lat_mod-lat_obs));
                dis=sqrt((lon_mod(ilon)-lon_obs).^2 + (lat_mod(ilat)-lat_obs).^2)*111;
                %display(['Distance between obs ',file,' and grid point is: ',num2str(dis),' km'])
            
                time=squeeze(double(ncread(fname,'time')))./24+t; time=time(1:ltime-1);
                if strncmp(ww3pre,'ww3p',4)
                elseif strncmp(ww3pre,'ww3g',4)
                  hs=squeeze(double(ncread(fname,'hsig',          [ilon ilat 1],[1 1 ltime]))); hs=hs(1:end-1);
                  tp=squeeze(double(ncread(fname,'tpeak',         [ilon ilat 1],[1 1 ltime]))); tp=tp(1:end-1);
                  %tm01=squeeze(double(ncread(fname,'tmean01',     [ilon ilat 1],[1 1 ltime]))); tm01=tm01(1:end-1);
                  %tm02=squeeze(double(ncread(fname,'tmean02',     [ilon ilat 1],[1 1 ltime]))); tm02=tm02(1:end-1);
                  %tm=squeeze(double(ncread(fname,'tmean',         [ilon ilat 1],[1 1 ltime]))); tm=tm(1:end-1);
                  pd=squeeze(double(ncread(fname,'peak_direction',[ilon ilat 1],[1 1 ltime]))); pd=pd(1:end-1);
                  %ds=squeeze(double(ncread(fname,'directional_spread',[ilon ilat 1],[1 1 ltime]))); ds=ds(1:end-1);
                  %if strcmp(expt,'NZWAVE-HR')
                  %  wlv=squeeze(double(ncread(fname,'wlv',          [ilon ilat 1],[1 1 ltime-1]))); 
                  %  model(ke).wlv=wlv;
                  %  ucur=squeeze(double(ncread(fname,'ucur',        [ilon ilat 1],[1 1 ltime-1]))); 
                  %  model(ke).ucur=ucur;
                  %  vcur=squeeze(double(ncread(fname,'vcur',        [ilon ilat 1],[1 1 ltime-1]))); 
                  %  model(ke).vcur=vcur;
                  %else
                  %  wlv=nan(size(hs));
                  %  ucur=nan(size(hs));
                  %  vcur=nan(size(hs));
                  %end
                end
                %wave energy (we) = (1/16)*rho*g*(hs.^2); % https://www.coastalwiki.org/wiki/Shallow-water_wave_theory#Introduction
                we=(1/16)*1025*9.8*(hs.^2);    

                display(['Saving:   ',filename])
                save(filename,'time','tp','hs','pd')%,'ds','tm01','tm02','pd','tm','wlv','ucur','vcur','we')
              end %if ck_mod_mat==0 || exist(filename)~=2
              %toc
            end % for fobss=stations
        
          else

            ptime=datestr(t,'YYYY/mm/DD/HH/');
            ftime=datestr(t,'YYYYmmDDHH');
          
            pname=[path_expt,ptime];
          
            fname=[pname,ww3pre,'_',ftime,'-utc_',gnames{ig},'.nc'];
            %fname=[pname,'ww3g_',ftime,'-utc_',gnames{1},'.nc'];
            display(['Loading: ',fname]);
    
            if checkin==0; % t==time_lima(1)
              checkin=1; 
              lon_mod=double(ncread(fname,'lon'));
              lat_mod=double(ncread(fname,'lat'));
              % coastal points
              if plot_coastm==1
                lon_obs=lon_mod(points(fobs,1));
			  				lat_obs=lat_mod(points(fobs,2));
              end
              [dif ilon]=nanmin(abs(lon_mod-lon_obs));
              [dif ilat]=nanmin(abs(lat_mod-lat_obs));
    
              if strcmp(file,'SteepHead_SP_FullRecord_QC')
                %ilon=ilon+3;
              elseif strcmp(file,'wairewa_lake_forsyth') & ig==2 % NZWAVE
                ilat=ilat-2;
              end
              ilon, ilat
              dis=sqrt((lon_mod(ilon)-lon_obs).^2 + (lat_mod(ilat)-lat_obs).^2)*111;
              display(['Distance between obs ',file,' and grid point is: ',num2str(dis),' km'])
            end
            
            time=squeeze(double(ncread(fname,'time')))./24+t; time=time(1:ltime-1);
            if strncmp(ww3pre,'ww3p',4)
              hs=squeeze(double(ncread(fname,'hsig',          [ilon 1],[1 ltime]))); hs=hs(1:end-1)';
              tp=squeeze(double(ncread(fname,'tpeak',         [ilon 1],[1 ltime]))); tp=tp(1:end-1)';
              tm01=squeeze(double(ncread(fname,'tmean01',     [ilon 1],[1 ltime]))); tm01=tm01(1:end-1)';
              tm=squeeze(double(ncread(fname,'tmean',         [ilon 1],[1 ltime]))); tm=tm(1:end-1)';
              pd=squeeze(double(ncread(fname,'peak_direction',[ilon 1],[1 ltime]))); pd=pd(1:end-1)';
            elseif strncmp(ww3pre,'ww3g',4)
              hs=squeeze(double(ncread(fname,'hsig',          [ilon ilat 1],[1 1 ltime]))); hs=hs(1:end-1);
              tp=squeeze(double(ncread(fname,'tpeak',         [ilon ilat 1],[1 1 ltime]))); tp=tp(1:end-1);
              tm01=squeeze(double(ncread(fname,'tmean01',     [ilon ilat 1],[1 1 ltime]))); tm01=tm01(1:end-1);
              tm02=squeeze(double(ncread(fname,'tmean02',     [ilon ilat 1],[1 1 ltime]))); tm02=tm02(1:end-1);
              tm=squeeze(double(ncread(fname,'tmean',         [ilon ilat 1],[1 1 ltime]))); tm=tm(1:end-1);
              pd=squeeze(double(ncread(fname,'peak_direction',[ilon ilat 1],[1 1 ltime]))); pd=pd(1:end-1);
              md=squeeze(double(ncread(fname,'mean_direction',[ilon ilat 1],[1 1 ltime]))); md=md(1:end-1);
              ds=squeeze(double(ncread(fname,'directional_spread',[ilon ilat 1],[1 1 ltime]))); ds=ds(1:end-1);
              try
                vinfo = ncinfo(fname, 'wlv'); wlv_exist=1;
              catch; wlv_exist=0; disp('wlv not found');
              end 
              if wlv_exist==1; %strcmp(expt,'NZWAVE-HR')
                wlv=squeeze(double(ncread(fname,'wlv',          [ilon ilat 1],[1 1 ltime-1]))); 
                model(ke).wlv=wlv;
                ucur=squeeze(double(ncread(fname,'ucur',        [ilon ilat 1],[1 1 ltime-1]))); 
                model(ke).ucur=ucur;
                vcur=squeeze(double(ncread(fname,'vcur',        [ilon ilat 1],[1 1 ltime-1]))); 
                model(ke).vcur=vcur;
              else
                wlv=nan(size(hs));
                ucur=nan(size(hs));
                vcur=nan(size(hs));
              end
            end
            %wave energy (we) = (1/16)*rho*g*(hs.^2); % https://www.coastalwiki.org/wiki/Shallow-water_wave_theory#Introduction
            we=(1/16)*1025*9.8*(hs.^2);    

            display(['Saving: ',filename])
            save(filename,'time','tp','pd','md','ds','tm01','tm02','hs','tm','wlv','ucur','vcur','we')
  
          end % if ck_mod_mat==1 && exist(filename)==2
          
          if plot_coastm==0 
            time_mod=[time_mod;time];
            hs_mod=[hs_mod;hs];
            tp_mod=[tp_mod;tp];
            tm01_mod=[tm01_mod;tm01];
            tm02_mod=[tm02_mod;tm02];
            tm_mod=[tm_mod;tm];
            pd_mod=[pd_mod;pd];
            md_mod=[md_mod;md];
            ds_mod=[ds_mod;ds];
            wlv_mod=[wlv_mod;wlv];
            if exist('ucur','var')==1
              ucur_mod=[ucur_mod;ucur];
	  	        vcur_mod=[vcur_mod;vcur];
		          we_mod=[we_mod;we];
            end  

            clear model
            model.data(:,1)=tp_mod;
            model.data(:,2)=pd_mod;
            model.data(:,6)=hs_mod;
            model.data(:,4)=tm01_mod;
            model.data(:,5)=tm02_mod;
            model.data(:,3)=ds_mod;
		        model.data(:,7)=md_mod;
            model.data(:,8)=tm_mod;
            model.data(:,9)=wlv_mod; % we_mod
            if exist('ucur','var')==1
		          model.data(:,10)=ucur_mod;
		          model.data(:,11)=vcur_mod;
		          model.data(:,12)=we_mod;
	          end
  
          end
 
        end % for t=time_lima

        if plot_coastm==1
          display(['Finished processing coastal stations from: ',datestr(time_lima(1)),' to ',datestr(time_lima(end))]) 
          toc
          return
        end

        % interpolating to half hourly intervals
        interp_half_hour=0;
        if interp_half_hour==1;
          time_modi=time_mod(1):1/24/2:time_mod(end);
          for i=1:size(model.data,2)
        	  datai(:,i)=interp1(time_mod,model.data(:,i),time_modi);
        	end
          time_mod=time_modi;
        	model.data=datai;
        end

        if save_mat==1
          fda=[filename(1:end-14),datestr(time_lima(1),'YYYYmmDDHH'),'_',datestr(time_lima(end),'YYYYmmDDHH'),'.mat'];
          display(['Saving: ',fda]);
          save(fda,'time_mod','tp_mod','pd_mod','md_mod','ds_mod','tm01_mod','tm02_mod','hs_mod','tm_mod','wlv_mod','ucur_mod','vcur_mod')% ,'we',model)
        end

      end % if the whole time series file is saved just load it

      clear model
      model.data(:,1)=tp_mod;
      model.data(:,2)=pd_mod;
      model.data(:,6)=hs_mod;
      model.data(:,4)=tm01_mod;
      model.data(:,5)=tm02_mod;
      model.data(:,3)=ds_mod;
		  model.data(:,7)=md_mod;
      model.data(:,8)=tm_mod;
      if exist('ucur','var')==1
        model.data(:,9)=wlv_mod; % we_mod
		    model.data(:,10)=ucur_mod;
		    model.data(:,11)=vcur_mod;
		    %model.data(:,12)=we_mod;
	    end
 
      if save_csv==1

    
        fda=[filename(1:end-14),datestr(time_lima(1),'YYYYmmDDHH'),'_',datestr(time_lima(end),'YYYYmmDDHH'),'.csv'];
        display(['Saving: ',fda]);
    
        if exist(fda)==2
        system(['rm -rf ',fda]);
        end
    
        csvdata=model.data;
        csvdata=num2cell(csvdata);
        for i=1:length(time_mod); csvtime{i,1}=datestr(time_mod(i),'yyyy-mm-ddTHH:MM:SS'); end
        if exist('ucur','var')==1
          csvhead={'date','Peak Period','Peak Direction','Directional spread','Mean period Tm01','Mean period Tm02','Hs','Meand Direction','Mean period Tm','Water level','U current','V current','Wave Energy (J/m2)'};
        else
          %csvhead={'date','Peak Period','Peak Direction','Hs','Mean period Tm01','Mean period Tm02','Directional spread','Meand Direction','Mean period Tm'};
          csvhead={'date','Peak Period','Peak Direction','Hs'};
        end

        % concatenating str matrices
        csvall=[csvtime,csvdata];
        %csvall=[csvhead;csvall];
        %csvwrite(fda,csvall);
        T = cell2table(csvall,'VariableNames',csvhead);
        writetable(T,fda)
        fclose('all');

      end

      % Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      tnames={'Peak Period',    'Peak Direction',    'Directional spread',   'Tm01','Tm02','Significant Wave Height (Hs)','Qp','Tm','Water Level','U current','V current','Wave Energy (J/m2)'};
      %tnames={'Peak Period','Peak Direction','Hs','Mean period Tm01','Mean period Tm02','Directional spread','Meand Direction','Mean period Tm','Water level','U current','V current','Wave Energy (J/m2)'};
      for i=1:length(dcol)
        
        subplot(length(dcol),1,i);
        set(gca,'fontsize',12,'fontweight','bold')
        hold on

        if ke==1 & plot_obs==1
  	      if obs_first==1;
				    if dcol(i)<=size(obs,2) 
              plot(time_obs,obs(:,dcol(i))','.','markersize',12,'color',[0 .7 0],'linewidth',2)
	   				end
          end
          if plot_mod==1
            for eix=1:length(expts); plot(nan,nan,'.','markersize',12,'color',colors{eix},'linewidth',2); 
            end
          end
       	  plot(nan,nan,'.','markersize',12,'color',[0 .7 0],'linewidth',2)
      	end

        plot(time_mod,model.data(:,dcol(i))','.','markersize',12,'color',colors{ke},'linewidth',2)

        if obs_first==0; 
        	if ke==length(expts) & plot_obs==1
	    			if dcol(i)<=size(obs,2) 
            	plot(time_obs,obs(:,dcol(i))','.','markersize',12,'color',[0 .7 0],'linewidth',2)
	   			  end
        	end
        end

        % computing harmonic analysis of wave variables with t_tide and plotting its timeseries
        if plot_harm==1
          const_name={'M2 ','S2 '};%,'S2','N2','K2','K1','O1','P1','Q1','MF','MM','M4','MS4','MN4'};
          if plot_obs==1 & ke==1 
            % cropping observations according to time_lima
            [dif it]=nanmin(abs(time_obs-time_lima(1))); [dif ft]=nanmin(abs(time_obs-time_lima(end)));
            time_harm=time_obs(it):1/24/2:time_obs(ft);
            data_obs=interp1(time_obs,obs(1:end,dcol(i)),time_harm); % [dif it]=nanmin(abs(time_obs-time_lima(1))); [dif ft]=nanmin(abs(time_obs-time_lima(end)));
				    [tidestruc,data_oharm]=t_tide(data_obs,'interval',(time_harm(2)-time_harm(1))*24,'start time',time_harm(1),'latitude',lat_obs,'synthesis',1); %,'output','none'););
				    %[tidestruc,data_harm]=t_predic(time_obs,tidestruc,'synthesis',1);
				    plot(time_harm,data_oharm,'.','markersize',12,'color',[1 0 0],'linewidth',2)
            % saving amplitude of selected constituents (const_name) 
					  [in,jn]=size(tidestruc.name);
            for ii=1:length(const_name)
              for iii=1:in
                if strcmp([tidestruc.name(iii),tidestruc.name(iii+in),tidestruc.name(iii+(in*2))],const_name{ii})==1; c_ind=iii;
                end
              end
              tide_amp(ii,ke)=tidestruc.tidecon(c_ind,1); % mod_amp(lat,lon,constante)
            end
          end
          if plot_mod==1 
            % cropping observations according to time_lima
            data_mod=model.data(:,dcol(i)); 
				    [tidestruc,data_mharm]=t_tide(data_mod,'interval',(time_mod(2)-time_mod(1))*24,'start time',time_mod(1),'latitude',lat_obs,'synthesis',1); %,'output','none'););
				    %[tidestruc,data_harm]=t_predic(time_obs,tidestruc,'synthesis',1);
				    plot(time_mod,data_mharm,'.','markersize',12,'color',[0 0 0],'linewidth',2)
            % saving amplitude of selected constituents (const_name) for each experiment (ke)
					  [in,jn]=size(tidestruc.name);
            for ii=1:length(const_name)
              for iii=1:in
                if strcmp([tidestruc.name(iii),tidestruc.name(iii+in),tidestruc.name(iii+(in*2))],const_name{ii})==1; c_ind=iii;
                end
              end
              tide_amp(ii,ke+1)=tidestruc.tidecon(c_ind,1); % mod_amp(lat,lon,constante)
            end
          end
        end

        if set_xlim==1;
          xlim([time_lima(1) time_lima(end)])
        end
      
        %if i==3
          datetick('x','dd/mmm/yyyy','keeplimits')
        %end
        title([replace(file,'_',' '),' ',tnames{dcol(i)}])
        grid('on')
        if plot_obs==1 
	        modeli=interp1(time_mod,model.data(:,dcol(i)),time_obs);
          if (dcol(i)~=2 || dcol(i)~=1) & i==1
            mod_rmse=nanrmse(obs(:,dcol(i)),modeli');
            mod_corr=nancorr(obs(:,dcol(i)),modeli');
            mod_bias=nanmean(modeli')-nanmean(obs(:,dcol(i)));
            if add_stats
              legh=[legh,{[expt,' rmse=',num2str(mod_rmse,'%.2f'),' corr=',num2str(mod_corr,'%.2f'),', bias=',num2str(mod_bias,'%.2f')]}];
            else
              legh=[legh,{[expt]}];%,' rmse=',num2str(mod_rmse,'%.2f'),' corr=',num2str(mod_corr,'%.2f'),', bias=',num2str(mod_bias,'%.2f')]}];
            end
				  elseif dcol(i)==2 
				    ylim([0 360])
          end
        else
          legh=[legh,{[expt]}];%,' rmse=',num2str(mod_rmse,'%.2f'),' corr=',num2str(mod_corr,'%.2f'),', bias=',num2str(mod_bias,'%.2f')]}];
        end
        
      end % for i=1:length(dcol)
    
    end % for expt=expts

    if obs_first==0; 
      legh=[legh,{['Obs.']}]; %,' rmse=',num2str(mod_rmse,'%.2f'),' corr=',num2str(mod_corr,'%.2f'),', bias=',num2str(mod_bias,'%.2f')]}];
    end

    %if ~isempty(find(dcol==6))
    %  subplot(length(dcol),1,find(dcol==6));
      subplot(length(dcol),1,1);
      legend([legh],'location','best')
    %end


    path_dm=strcat(path_fig,expt,'/');
    system(['mkdir -p ',path_dm]);
    figname=[path_dm,file,'_wave_series_',replace(tnames{dcol(i)},' ','_'),'_',datestr(time_lima(1),'YYYY_mm_DD_HH_MM'),'_',datestr(time_lima(end),'YYYY_mm_DD_HH_MM')];
    display(['Plotting: ',figname]);
    if save_fig==1
      display(['Saving: ',figname]);
      export_fig(gcf,[figname,'.png'],'-png','-r150');
      saveas(gcf,[figname,'.fig'],'fig')
      %close
      %clf('reset')
    end


    if plot_harm==1
       disp(' ')
       disp([replace(file,'_',' '),' ',tnames{dcol(i)}])
      for ii=1:length(const_name)
        if plot_obs
        	display(['Amplitude of ',const_name{ii},' in obs: ',num2str(tide_amp(ii,1))])
        end
        for ke=1:length(expts)
      	  display(['Amplitude of ',const_name{ii},' in ',expts{ke},': ',num2str(tide_amp(ii,ke+1))])
      	end
      end
    end

  end % if plot_series

  if plot_rose
      
    legh={'Obs'}; legt={'Obs'};
    
    ke=0;
    for expt=expts
      ke=ke+1;
    
      expt=expt{1};
      [ig,ltime]=grab_nc_sufix(expt);
    
      path_expt=[path_source,expt,'/']; % GLOBALWAVE/'];%2018/01/05/00'];
      
      path_dm=[path_expt,'matlab/'];
      system(['mkdir -p ',path_dm]);
      %filename=[path_dm,file,'_',datestr(time_lima(1),'YYYYmmDDHH'),'_',datestr(time_lima(end),'YYYYmmDDHH'),'.mat'];
      hs_mod=[];tp_mod=[]; pd_mod=[];time_mod=[];
      tm01_mod=[];tm_mod=[]; pd_mod=[];time_mod=[];
      tm02_mod=[];ds_mod=[]; wlv_mod=[];time_mod=[];
    
      checkin=0;  
      % Loop in model time
      for t=time_lima
        filename=[path_dm,file,'_',datestr(t,'YYYYmmDDHH'),'.mat'];
    
        if ck_mod_mat==1 && exist(filename)==2
      
          display(['Loading: ',filename])
          load(filename)%,'time_mod','model')
      
        else

          display(['File doenst exist: ',filename])
          display(['Run with plot_series=1 first'])
          return
  
        end % if ck_mod_mat==1 && exist(filename)==2
    
        time_mod=[time_mod;time];
        hs_mod=[hs_mod;hs];
        tp_mod=[tp_mod;tp];
        tm01_mod=[tm01_mod;tm01];
        tm02_mod=[tm02_mod;tm02];
        tm_mod=[tm_mod;tm];
        pd_mod=[pd_mod;pd];
        ds_mod=[ds_mod;ds];
        wlv_mod=[wlv_mod;wlv];
        
      end % for t=time_lima

      clear model
     
      % 1 {'PeakPerid'}   2 {'PeakDirection'}   3 {'Directional spread'}   4 {'Tm01'}   5 {'Tm02'}   6 {'Hs'}   7 {'Qp'}  8 Tm
      model.data(:,9)=wlv_mod;
      model.data(:,8)=tm_mod;
      model.data(:,1)=tp_mod;
      model.data(:,2)=pd_mod;
      model.data(:,3)=ds_mod;
      model.data(:,4)=tm01_mod;
      model.data(:,5)=tm02_mod;
      model.data(:,6)=hs_mod;

      % wave rose options        
      Option = {'vWinds', [0:6],'anglenorth', 0, 'angleeast', 90, 'labels', {'N (0°)', 'E (90°)', 'S (180°)', 'W (270°)'}, 'freqlabelangle', 45, 'legendtype', 0};    
      %Option={'legendtype', 0, 'vWinds', [0:6]};

      tnames={'Peak Period',    'Peak Direction',    'Directional spread',   'Tm01','Tm02','Significant Wave Height','Qp','Tm','Water Level'};
      vnames={'Tp',             'Dp',                'Spr',                  'Tm01','Tm02','Hs',                         'Qp' ,'Tm','Wlv'};
      for i=1:length(dcol)

        if portrait
          scrsz=[2    42   958   953];
        else
          scrsz=[1 41 1920 962];
        end

        h=figure(i);
        set(h,'position',scrsz,'color',[1 1 1],'visible','on')

        if length(runs)==1
          lsub=2; csub=1;
        elseif length(runs)==2
          lsub=3; csub=1;
        elseif mod(length(runs),2)==1
          lsub=round(length(runs)/2); csub=round(length(runs)/2)-1;
        elseif mod(length(runs),2)==0
          lsub=round(length(runs)/2); csub=round(length(runs)/2);
        end
         
        if portrait
          ax=subplot(lsub,csub,ke);
        else
          ax=subplot(csub,lsub,ke);
        end
        Options = [Option, {'axes', ax, 'cmap', 'invbone','TitleString', {['(a) ',replace(file,'_',' '), ' obs. ',tnames{dcol(i)}];''},...
        'LegendType', 0, 'LegendPosition', 'north','LegendOrientation', 'vertical', 'LabLegend',vnames{dcol(i)} , 'LegendVariable', vnames{dcol(i)}, 'scalefactor', 1.0}];%, 'vWinds', [0:.5 1 1.2 2 5]}];
        if ke==1
          [dif is]=nanmin(abs(time_obs-time_lima(1)));
          [dif ie]=nanmin(abs(time_obs-time_lima(end)));
          dp_obs=obs(is:ie,2)';
          data_obs=obs(is:ie,dcol(i))';
          [figure_handle, count, speeds, directions, Table] = WindRose(dp_obs,data_obs,Options);
          %title('Obs')
        end
        if portrait
          ax=subplot(lsub,csub,ke+1);
        else
          ax=subplot(csub,lsub,ke+1);
        end
        lett={'(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)','(m)','(n)','(o)','(p)','(q)','(r)','(s)','(t)','(u)','(v)','(w)','(x)','(y)','(z)'};
        ltype=0; if length(runs)==ke; ltype=2; end
        Options = [Option, {'axes', ax, 'cmap', 'invbone','TitleString', {[lett{ke},' ',expt];''},...
        'LegendType', ltype, 'LegendPosition', 'west','LegendOrientation','vertical', 'LabLegend', vnames{dcol(i)} , 'LegendVariable', vnames{dcol(i)}, 'scalefactor', 1.0}];
        [figure_handle, count, speeds, directions, Table] = WindRose(model.data(:,2)',model.data(:,dcol(i))',Options);
         
      end
    
    end
    
    if save_fig==1
      %export_fig(gcf,'steephead_period','-png','-r150');
    end

    path_dm=strcat(path_fig,'paper_1/');
    system(['mkdir -p ',path_dm]);
    figname=[path_dm,file,'_wave_rose_',replace(tnames{dcol(i)},' ','_'),'_',datestr(time_lima(1),'YYYY_mm_DD_HH_MM'),'_',datestr(time_lima(end),'YYYY_mm_DD_HH_MM'),'.png'];
    display(['Plotting: ',figname]);
    if save_fig==1
      display(['Saving: ',figname]);
      export_fig(gcf,figname,'-png','-r150');
      %close
      %clf('reset')
    end

  end % if plot_rose

  if plot_disp
  
    scrsz=[1 1 1366 768];
    %scrsz=get(0,'screensize');
    if portrait==1
      scrsz=[2    42   958   953];
    else
      scrsz=[1 1 1910 990];
    end

    figure('position',scrsz,'color',[1 1 1],'visible','on')
    hold on
    set(gca,'fontsize',12,'fontweight','bold')
      
    legh={'Obs'}; legt={'Obs'};
     
    ke=0;
    for expt=expts
      ke=ke+1;
    
      expt=expt{1};
      [ig,ltime]=grab_nc_sufix(expt);
    
      path_expt=[path_source,expt,'/']; % GLOBALWAVE/'];%2018/01/05/00'];
      
      path_dm=[path_expt,'matlab/'];
      system(['mkdir -p ',path_dm]);
      %filename=[path_dm,file,'_',datestr(time_lima(1),'YYYYmmDDHH'),'_',datestr(time_lima(end),'YYYYmmDDHH'),'.mat'];
      hs_mod=[];tp_mod=[]; pd_mod=[];time_mod=[];
      tm01_mod=[];tm_mod=[]; pd_mod=[]; we_mod=[];
      tm02_mod=[];ds_mod=[]; wlv_mod=[];time_mod=[];

      % checking if the whole time series is saved as matlab file
      filename=[path_dm,file,'_',datestr(time_lima(1),'YYYYmmDDHH'),'_',datestr(time_lima(end),'YYYYmmDDHH'),'.mat'];
      if exist(filename)==2

        display(['Loading: ',filename])
        load(filename); % loading long timeseries rather than daily files

      else
    
        checkin=0;  
        % Loop in model time
        for t=time_lima

          filename=[path_dm,file,'_',datestr(t,'YYYYmmDDHH'),'.mat'];
    
          if exist(filename)==2
        
            display(['Loading: ',filename])
            load(filename)%,'time_mod','model')
        
          else
            display(['NOT FOUND: ',filename])
            display(['Run with plot_series=1 first'])
            return
  
          end % if ck_mod_mat==1 && exist(filename)==2
    
          time_mod=[time_mod;time];
          hs_mod=[hs_mod;hs];
          tp_mod=[tp_mod;tp];
          tm01_mod=[tm01_mod;tm01];
          tm02_mod=[tm02_mod;tm02];
          tm_mod=[tm_mod;tm];
          pd_mod=[pd_mod;pd];
          ds_mod=[ds_mod;ds];
          wlv_mod=[wlv_mod;wlv];
          %we_mod=[we_mod;we];
          
          
          %path_dm=[path_expt,'matlab/'];
          %system(['mkdir -p ',path_dm]);
          %display(['Saving: ',filename])
          %save(filename,'time_mod','model')
        
        end % for t=time_lima

      end

      clear model
     
      % 1 {'PeakPerid'}   2 {'PeakDirection'}   3 {'Directional spread'}   4 {'Tm01'}   5 {'Tm02'}   6 {'Hs'}   7 {'Qp'}  8 Tm
      model.data(:,9)=wlv_mod; % we_mod
      model.data(:,8)=tm_mod;
      model.data(:,1)=tp_mod;
      model.data(:,2)=pd_mod;
      model.data(:,3)=ds_mod;
      model.data(:,4)=tm01_mod;
      model.data(:,5)=tm02_mod;
      model.data(:,6)=hs_mod;
 
      tnames={'Peak Period (s)',    'Peak Direction (^o)',    'Directional spread (^o)',   'Tm01 (m)','Tm02 (m)','Hs (m)','Qp (^o)','Tm (s)','Water Level (m)'};
      for i=1:length(dcol)
        
        modeli=interp1(time_mod,model.data(:,dcol(i)),time_obs)';
        
        if length(runs)==1
          %subplot(length(dcol),1,i);
          lsub=1; csub=1;
        elseif length(runs)==2
          %subplot(length(dcol),1,i);
          lsub=2; csub=1;
        elseif mod(length(runs),2)==1
          lsub=round(length(runs)/2); csub=round(length(runs)/2)-1;
        elseif mod(length(runs),2)==0
          lsub=round(length(runs)/2); csub=round(length(runs)/2);
        end

        if portrait==1
          ax=subplot(lsub,csub,ke);
        else
          ax=subplot(csub,lsub,ke);
        end
        yplot=modeli; xplot=obs(:,dcol(i));

        if correct_model==-1;

          display(['Plotting scatter plot']); tic
          scatter_kde(xplot, yplot, 'filled', 'MarkerSize', 10); toc %,'MarkerEdgeColor','k'); toc
          colormap(flipud(cmocean('matter')))
          set(gca,'fontsize',12,'fontweight','bold')
          hold on

          % Add Color bar
          %cb = colorbar();
          %cb.Label.String = 'Probability density estimate';

          hold on
          plot([0 9],[0 9],'--b','linewidth',2)

          inan=isnan(modeli); obs(inan,dcol(i))=nan;  
          inan=isnan(obs(:,dcol(i))); modeli(inan)=nan;
          modi=modeli; inan=isnan(modi); modi(inan)=[]; 
          obsi=obs(:,dcol(i)); inan=isnan(obsi); obsi(inan)=[];  
          
          p=polyfit(obsi,modi,1);
          %plot([0 9],[0 9*p(1)+p(2)],'--b','linewidth',2)

          if ke==length(expts)
	        	h=colorbar;
          	set(h,'position',[0.915 0.24 0.011 0.44]) % [left bottom width height]
            set(get(h,'ylabel'),'string','Probability density estimate','fontsize',10,'fontweight','bold');
          else
            colorbar('off')
          end
          %colorbar
          caxis([0 1])
          %if set_xlim==1;
          %  xlim([time_lima(1)-2 time_lima(end)+2])
          %end
          %if i==3
          %end
          title([expt])% ,': ',replace(file,'_',' '),' ',tnames{dcol(i)}])
          grid('on')
          xlim([0 9]) 
          ylim([0 9])
          %axis('equal')

        elseif correct_model==0;

          hold on
          plot(xplot,yplot,'.','color',[1 1 1],'linewidth',2)
          plot(xplot,yplot,'.','color',[1 1 1],'linewidth',2)
          plot([0 9],[0 9],'r','linewidth',2)

          inan=isnan(modeli); obs(inan,dcol(i))=nan;  
          inan=isnan(obs(:,dcol(i))); modeli(inan)=nan;
          modi=modeli; inan=isnan(modi); modi(inan)=[]; 
          obsi=obs(:,dcol(i)); inan=isnan(obsi); obsi(inan)=[];  
          
          p=polyfit(obsi,modi,1);
          plot([0 9],[0 9*p(1)+p(2)],'color',[0 0 1],'linewidth',2)

          bh=binscatter(xplot,yplot,90);
          set(gca,'fontsize',12,'fontweight','bold')
          hold on
          colormap(flipud(cmocean('matter')))

          if ke==length(expts)
	          h=colorbar;
          	set(h,'position',[0.915 0.24 0.011 0.44])
            set(get(h,'ylabel'),'string','Bin counts','fontsize',14,'fontweight','bold');
          else
            colorbar('off')
          end
          %colorbar
          caxis([0 250])
          %if set_xlim==1;
          %  xlim([time_lima(1)-2 time_lima(end)+2])
          %end
          %if i==3
          %end
          title([expt])% ,': ',replace(file,'_',' '),' ',tnames{dcol(i)}])
          grid('on')
          %xlim([0 9]) 
          %ylim([0 9])
          axis('equal')

        elseif correct_model>=1;
          set(gca,'fontsize',12,'fontweight','bold')
          hold on

          plot(xplot,yplot,'.','color',[1 1 1],'linewidth',2)
          plot(xplot,yplot,'.','color',[1 1 1],'linewidth',2)
          plot([0 9],[0 9],'r','linewidth',2)
          plot(xplot,yplot,'.','color',[0 0 0],'linewidth',2)
           
          %if set_xlim==1;
          %  xlim([time_lima(1)-2 time_lima(end)+2])
          %end
          %if i==3
          %end
          title([expt])% ,': ',replace(file,'_',' '),' ',tnames{dcol(i)}])
          grid('on')
          xlim([0 9]) 
          ylim([0 9])
 
          inan=isnan(modeli); obs(inan,dcol(i))=nan;  
          inan=isnan(obs(:,dcol(i))); modeli(inan)=nan;
          modi=modeli; inan=isnan(modi); modi(inan)=[]; 
          obsi=obs(:,dcol(i)); inan=isnan(obsi); obsi(inan)=[];  
          
          p=polyfit(obsi,modi,1);
          plot([0 9],[0 9*p(1)+p(2)],'k','linewidth',2)

          % correcting model based on the linear relationship
          %cfac=(obsi-p(2))/p(1);
          if correct_model>1
            modc=modi.*(1+(1-p(1)));% + p(2);
            p=polyfit(obsi,modc,1)
            plot(obsi,modc,'.','color',[0 0 1],'linewidth',2)
            plot([0 9],[0 9*p(1)+p(2)],'b','linewidth',2)
            xplot=modc; yplot=obsi;
          end

        end

        ylabel(['Forecasted ',tnames{dcol(i)}])      
        xlabel(['Observed ',tnames{dcol(i)}])      

        %axis('equal')  
        %if dcol(i)==6
        mod_rmse=nanrmse(xplot,yplot);
        mod_corr=nancorr(xplot,yplot);
        mod_bias=nanmean(yplot)-nanmean(xplot);
        if correct_model<=0;
          yd=8.5; dsz=1;
          text(-.0,yd,['y=',num2str(p(1),'%.2f'),'*x',num2str(p(2),'%+.2f')],'fontsize',10,'fontweight','bold');
          text(-.0,yd-dsz,['rmse=',num2str(mod_rmse,'%.2f')]                ,'fontsize',10,'fontweight','bold');
          text(-.0,yd-dsz*2,['bias=',num2str(mod_bias,'%.2f')]              ,'fontsize',10,'fontweight','bold');
          text(-.0,yd-dsz*3,['r=',num2str(mod_corr,'%.2f')]                 ,'fontsize',10,'fontweight','bold');
        elseif correct_model==1;
          %  legh=[legh,{[expt,' rmse=',num2str(mod_rmse,'%.2f'),' corr=',num2str(mod_corr,'%.2f'),', bias=',num2str(mod_bias,'%.2f')]}];
          legh=[{['y=',num2str(p(1),'%.2f'),'*x',num2str(p(2),'%+.2f')];['rmse=',num2str(mod_rmse,'%.2f'),' corr=',num2str(mod_corr,'%.2f'),', bias=',num2str(mod_bias,'%.2f')]}];
          legend([legh],'location','south')
        end
        box('on') 
      end
    end
    
    if ~isempty(find(dcol==6))
      %subplot(length(dcol),1,find(dcol==6));
      %legend([legh],'location','best')
    end
    
    path_dm=strcat(path_fig,expt,'/');
    system(['mkdir -p ',path_dm]);
    figname=[path_dm,file,'_',replace(tnames{dcol(i)},' ','_'),'_',datestr(time_lima(1),'YYYY_mm_DD_HH_MM'),'_',datestr(time_lima(end),'YYYY_mm_DD_HH_MM'),'.png'];
    display(['Plotting: ',figname]);
    if save_fig==1
      display(['Saving: ',figname]);
      export_fig(gcf,figname,'-png','-r150');
      %close
      %clf('reset')
    end
    
  end % if plot_disp
    
    
  if plot_eva
  
    scrsz=[1 1 1366 768];
    %scrsz=get(0,'screensize');
    %scrsz=[2    42   958   953];
    scrsz=[1 1 1910 990];
    
    figure('position',scrsz,'color',[1 1 1],'visible','on')
    hold on
    set(gca,'fontsize',12,'fontweight','bold')
      
    legh={'Obs'}; legt={'Obs'};
     
    ke=0;
    for expt=expts
      ke=ke+1;
    
      expt=expt{1};
      [ig,ltime]=grab_nc_sufix(expt);
    
      path_expt=[path_source,expt,'/']; % GLOBALWAVE/'];%2018/01/05/00'];
      
      path_dm=[path_expt,'matlab/'];
      system(['mkdir -p ',path_dm]);
      %filename=[path_dm,file,'_',datestr(time_lima(1),'YYYYmmDDHH'),'_',datestr(time_lima(end),'YYYYmmDDHH'),'.mat'];
      hs_mod=[];tp_mod=[]; pd_mod=[];time_mod=[];
      tm01_mod=[];tm_mod=[]; pd_mod=[]; we_mod=[];
      tm02_mod=[];ds_mod=[]; wlv_mod=[];time_mod=[];
    
      checkin=0;  
      % Loop in model time
      for t=time_lima

        filename=[path_dm,file,'_',datestr(t,'YYYYmmDDHH'),'.mat'];
    
        if ck_mod_mat==1 && exist(filename)==2
      
          display(['Loading: ',filename])
          load(filename)%,'time_mod','model')
      
        else
					
          display(['NOT FOUND: ',filename])
					return

        end % if ck_mod_mat==1 && exist(filename)==2
    
        time_mod=[time_mod;time];
        hs_mod=[hs_mod;hs];
        tp_mod=[tp_mod;tp];
        tm01_mod=[tm01_mod;tm01];
        tm02_mod=[tm02_mod;tm02];
        tm_mod=[tm_mod;tm];
        pd_mod=[pd_mod;pd];
        ds_mod=[ds_mod;ds];
        wlv_mod=[wlv_mod;wlv];
        %we_mod=[we_mod;we];
        
        clear model
     
        % 1 {'PeakPerid'}   2 {'PeakDirection'}   3 {'Directional spread'}   4 {'Tm01'}   5 {'Tm02'}   6 {'Hs'}   7 {'Qp'}  8 Tm
        model.data(:,9)=wlv_mod; % we_mod
        model.data(:,8)=tm_mod;
        model.data(:,1)=tp_mod;
        model.data(:,2)=pd_mod;
        model.data(:,3)=ds_mod;
        model.data(:,4)=tm01_mod;
        model.data(:,5)=tm02_mod;
        model.data(:,6)=hs_mod;
        
        %path_dm=[path_expt,'matlab/'];
        %system(['mkdir -p ',path_dm]);
        %display(['Saving: ',filename])
        %save(filename,'time_mod','model')
      
      end % for t=time_lima
    
      tnames={'Peak Period',    'Peak Direction',    'Directional spread',   'Tm01','Tm02','Significant Wave Height (Hs)','Qp','Tm','Water Level'};
      for i=1:length(dcol)
        
        %cdfplot(hs_obs) % PF cdf = cumulative density function. find suitable thresholds to constrain GPthres.m. This is just a rough guide. looks like something between 2 and 6 metres would work (between 0.5 and 1)
        tolerance =0.0001; % sometimes gpfit doesn't like having zero values, so add a tiny (inconsequential) amount
        tthr = 3; % peaks over threshold (POT) separated by at least 3 days. Most weather systems in NZ are separated by 4-7 days
        ARI = [1 2 5 10 20 50 100]; % This specifies the Average Recurrence Intervals (ARI) for output. Change these to what you like.
         
        disp(['GPD fitted to peaks over threshold']);
        disp(['DO NOT click on figures while plots are being made']);
        %figure; clf
        %semilogx([0 0],[0 0],'r',[0 0],[0 0],'b'); hold on;
        if ke==1

          ax=figure; % subplot(lsub,csub,ke);
          set(gca,'fontsize',12,'fontweight','bold')
          hold on
          fign=99;
 
          data_obs=obs(:,dcol(i)); 
          data = data_obs(:); time = time_obs(:);
          ythr = GPthres(time,data+tolerance,2,6,0.1,tthr,1); % Automate threshold selection
          disp(['GPD threshold selected: H_s = ' num2str(ythr) ' m'])
          [GPD,gpp,gpci,Sy,ARI_POT,POT_sorted] = GPDfit(data+tolerance,time,tthr,ythr,ARI); % PF
          close
          figure(fign)
          semilogx([0 0],[0 0],'color',[0 .7 0]); hold on;
          semilogx([0 0],[0 0],colors{ke});
          if correct_model>=1
            semilogx([0 0],[0 0],colors{ke},'linewidth',2);
            hold on;
          end

          semilogx(ARI_POT,POT_sorted,'o','color',[0 .7 0]),
          hold on
          semilogx(ARI,GPD(1,:),'-', 'color',[0 .7 0])
          semilogx(ARI,GPD(2,:),'--','color',[0 .7 0])
          semilogx(ARI,GPD(3,:),'--','color',[0 .7 0])
          title(['GPD fitted to peaks over threshold. ',datestr(time_lima(1),'YYYY'),'-',datestr(time_lima(end),'YYYY')])
          xlabel('Average recurrence interval (years)'), ylabel('H_s')

        end

        data = model.data(:,dcol(i)); time = time_mod(:);
        data=data(:);
        ythr = GPthres(time,data+tolerance,2,6,0.1,tthr,1); % Automate threshold selection
        disp(['GPD threshold selected: H_s = ' num2str(ythr) ' m']);
        figure;
        [GPD,gpp,gpci,Sy,ARI_POT,POT_sorted] = GPDfit(data+tolerance,time,tthr,ythr,ARI); % PF
        close
        figure(fign)
        semilogx(ARI_POT,POT_sorted,'o','color',colors{ke}),
        hold on
        semilogx(ARI,GPD(1,:),'-', 'color',colors{ke})
        semilogx(ARI,GPD(2,:),'--','color',colors{ke})
        semilogx(ARI,GPD(3,:),'--','color',colors{ke})

        %return
 
        % Compare to annual maxima method. Don't need this but sometimes useful to
        %AMy = AMcalc(time_obs,hs_obs,1); AM = AMy(:,3);
        %[gev,gevci] = gevfit(AM);
        %GEV = gevinv(1-ARI2AEP(ARI),gev(1),gev(2),gev(3));
        %GEV(2,:) = gevinv(1-ARI2AEP(ARI),gevci(1,1),gevci(1,2),gevci(1,3));
        %GEV(3,:) = gevinv(1-ARI2AEP(ARI),gevci(2,1),gevci(2,2),gevci(2,3));
        %
        %figure;[AMari,AMsorted] = gringorten(AM,1,1);  hold on
        %semilogx(ARI,GEV(1,:),'k','LineWidth',2)
        %semilogx(ARI,GEV(2,:),'--k',ARI,GEV(3,:),'--k')
        %hold off
        %title('GEV fitted to annual maxima'), xlabel('Average recurrence interval  (years)'), ylabel('H_s')
        %
        %figure(2)
        %hold on
        %semilogx(ARI,GEV(1,:),'b'), semilogx(ARI,GEV(2,:),'--b',ARI,GEV(3,:),'--b'),
        %semilogx(AMari,AMsorted,'+b')
        %hold off
        if correct_model==0;
          legend('Obs.',expt);
          legend(replace(file,'_',' '),expt);

        elseif correct_model==1;
          modeli=interp1(time_mod,model.data(:,dcol(i)),time_obs)';
          xplot=modeli; yplot=obs(:,dcol(i));
 
          inan=isnan(modeli); obs(inan,dcol(i))=nan;  
          inan=isnan(obs(:,dcol(i))); modeli(inan)=nan;
          modi=modeli; inan=isnan(modi); modi(inan)=[]; 
          obsi=obs(:,dcol(i)); inan=isnan(obsi); obsi(inan)=[];  
          
          p=polyfit(obsi,modi,1);
          %plot([0 9],[0 9*p(1)+p(2)],'b','linewidth',2)

          % correcting model based on the linear relationship
          %cfac=(obsi-p(2))/p(1);
          %modc=modeli.*(1+(1-p(1)));% + p(2);
          modc=model.data(:,dcol(i)).*(1+(1-p(1)));% + p(2);
          %pc=polyfit(obsi,modc,1)
          %plot(obsi,modc,'.','color',[0 0 0],'linewidth',2)
          %plot([0 9],[0 9*pc(1)+pc(2)],'k','linewidth',2)

        elseif correct_model==2;

          modeli=interp1(time_mod,model.data(:,dcol(i)),time_obs)';
          xplot=modeli; yplot=obs(:,dcol(i));
 
          inan=isnan(modeli); obs(inan,dcol(i))=nan;  
          inan=isnan(obs(:,dcol(i))); modeli(inan)=nan;
          modi=modeli; inan=isnan(modi); modi(inan)=[]; 
          obsi=obs(:,dcol(i)); inan=isnan(obsi); obsi(inan)=[];  

          model_wave_height=modi; satellite_wave_height=obsi;
         
          % Assuming model_wave_height and satellite_wave_height are your data vectors
          mdl = fitlm(model_wave_height, satellite_wave_height);
          disp(mdl);
          % Plot the initial fit
          figure; hold on;
          plot(mdl);
          xlabel('Model Wave Height');
          ylabel('Observed Wave Height');
          title('Initial Linear Regression Fit');          
          plot([0 9],[0 9],'k','linewidth',2)
          % Corrected wave height using the linear model
          corrected_wave_height = predict(mdl, model_wave_height);
          scatter(corrected_wave_height, satellite_wave_height,'.k');
          new_mdl = fitlm(corrected_wave_height, satellite_wave_height);
          disp(new_mdl);

          %residuals = mdl.Residuals.Raw;
          %figure; hold on
          %plot(model_wave_height, residuals, 'o');
          %xlabel('Model Wave Height');
          %ylabel('Residuals');
          %title('Residuals vs. Model Wave Height');

          mdl_poly = fitlm(model_wave_height, satellite_wave_height, 'quadratic');
          disp(mdl_poly);
          % Plot the polynomial fit
          figure; hold on;
          plot(mdl_poly);
          xlabel('Model Wave Height');
          ylabel('Observed Wave Height');
          title('Quadratic Polynomial Regression Fit');
          plot([0 9],[0 9],'k','linewidth',2)
          % Corrected wave height using the linear model
          qcorrected_wave_height = predict(mdl_poly, model_wave_height);
          scatter(qcorrected_wave_height, satellite_wave_height,'.k');
          new_mdl = fitlm(qcorrected_wave_height, satellite_wave_height);
          disp(new_mdl);
          
          % Alternatively, if using a log-transformed model:
          %corrected_wave_height = exp(predict(mdl_log, log_model_wave_height));
          
          %figure;
          %scatter(corrected_wave_height, satellite_wave_height);
          %xlabel('Corrected Wave Height');
          %ylabel('Satellite Wave Height');
          %title('Corrected vs. Satellite Wave Height');
          %% Calculate the new R-squared value
          %new_mdl = fitlm(corrected_wave_height, satellite_wave_height);
          %disp(new_mdl.Rsquared.Adjusted);

          modc=corrected_wave_height;

return

        end


        if correct_model>=1;

          data=modc; time=time_mod;

          ythr = GPthres(time,data+tolerance,2,6,0.1,tthr,1); % Automate threshold selection
          disp(['GPD threshold selected: H_s = ' num2str(ythr) ' m']);
          figure;
          [GPD,gpp,gpci,Sy,ARI_POT,POT_sorted] = GPDfit(data+tolerance,time,tthr,ythr,ARI); % PF
          close
          figure(fign)
          semilogx(ARI_POT,POT_sorted,'o','color',colors{ke+1}),
          hold on
          semilogx(ARI,GPD(1,:),'-', 'color',colors{ke+1})
          semilogx(ARI,GPD(2,:),'--','color',colors{ke+1})
          semilogx(ARI,GPD(3,:),'--','color',colors{ke+1})
          legend('Obs.',expt,['corrected ',expt]);

        end

        disp(['GEV fitted to annual maxima']);
        disp(['DO NOT click on figures while plots are being made']);
        %figure; clf
        %semilogx([0 0],[0 0],'r',[0 0],[0 0],'b'); hold on;
        if ke==1

          ax=figure; % subplot(lsub,csub,ke);
          set(gca,'fontsize',12,'fontweight','bold')
          hold on
          fign=98;
          % Compare to annual maxima method. Don't need this but sometimes useful to
          data_obs=obs(:,dcol(i)); 
          data = data_obs(:); time = time_obs(:);
          AMy = AMcalc(time,data,1); AM = AMy(:,3);
          [gev,gevci] = gevfit(AM);
          GEV = gevinv(1-ARI2AEP(ARI),gev(1),gev(2),gev(3));
          GEV(2,:) = gevinv(1-ARI2AEP(ARI),gevci(1,1),gevci(1,2),gevci(1,3));
          GEV(3,:) = gevinv(1-ARI2AEP(ARI),gevci(2,1),gevci(2,2),gevci(2,3));

          figure;[AMari,AMsorted] = gringorten(AM,1,1);  hold on; 
					close
          semilogx(AMari,AMsorted,'xb')
          semilogx(ARI,GEV(1,:),'k','LineWidth',2)
          semilogx(ARI,GEV(2,:),'--k',ARI,GEV(3,:),'--k')
          hold off
          title('GEV fitted to annual maxima'), xlabel('Average recurrence interval  (years)'), ylabel('H_s')

          %
          %figure(2)
          %hold on
          %semilogx(ARI,GEV(1,:),'b'), semilogx(ARI,GEV(2,:),'--b',ARI,GEV(3,:),'--b'),
          %semilogx(AMari,AMsorted,'+b')
          %hold off
return
 

        end

        data = model.data(:,dcol(i)); time = time_mod(:);
        data=data(:);
        ythr = GPthres(time,data+tolerance,2,6,0.1,tthr,1); % Automate threshold selection
        disp(['GPD threshold selected: H_s = ' num2str(ythr) ' m']);
        figure;
        [GPD,gpp,gpci,Sy,ARI_POT,POT_sorted] = GPDfit(data+tolerance,time,tthr,ythr,ARI); % PF
        close
        figure(fign)
        semilogx(ARI_POT,POT_sorted,'o','color',colors{ke}),
        hold on
        semilogx(ARI,GPD(1,:),'-', 'color',colors{ke})
        semilogx(ARI,GPD(2,:),'--','color',colors{ke})
        semilogx(ARI,GPD(3,:),'--','color',colors{ke})

        %return
 
        % Compare to annual maxima method. Don't need this but sometimes useful to
        AMy = AMcalc(time_obs,hs_obs,1); AM = AMy(:,3);
        %[gev,gevci] = gevfit(AM);
        %GEV = gevinv(1-ARI2AEP(ARI),gev(1),gev(2),gev(3));
        %GEV(2,:) = gevinv(1-ARI2AEP(ARI),gevci(1,1),gevci(1,2),gevci(1,3));
        %GEV(3,:) = gevinv(1-ARI2AEP(ARI),gevci(2,1),gevci(2,2),gevci(2,3));
        %
        %figure;[AMari,AMsorted] = gringorten(AM,1,1);  hold on
        %semilogx(ARI,GEV(1,:),'k','LineWidth',2)
        %semilogx(ARI,GEV(2,:),'--k',ARI,GEV(3,:),'--k')
        %hold off
        %title('GEV fitted to annual maxima'), xlabel('Average recurrence interval  (years)'), ylabel('H_s')
        %
        %figure(2)
        %hold on
        %semilogx(ARI,GEV(1,:),'b'), semilogx(ARI,GEV(2,:),'--b',ARI,GEV(3,:),'--b'),
        %semilogx(AMari,AMsorted,'+b')
        %hold off
        if correct_model==0;
          legend('Obs.',expt);
        end

        if correct_model==1;
          modeli=interp1(time_mod,model.data(:,dcol(i)),time_obs)';
          xplot=modeli; yplot=obs(:,dcol(i));
 
          inan=isnan(modeli); obs(inan,dcol(i))=nan;  
          inan=isnan(obs(:,dcol(i))); modeli(inan)=nan;
          modi=modeli; inan=isnan(modi); modi(inan)=[]; 
          obsi=obs(:,dcol(i)); inan=isnan(obsi); obsi(inan)=[];  
          
          p=polyfit(obsi,modi,1);
          %plot([0 9],[0 9*p(1)+p(2)],'b','linewidth',2)

          % correcting model based on the linear relationship
          %cfac=(obsi-p(2))/p(1);
          %modc=modeli.*(1+(1-p(1)));% + p(2);
          modc=model.data(:,dcol(i)).*(1+(1-p(1)));% + p(2);
          %pc=polyfit(obsi,modc,1)
          %plot(obsi,modc,'.','color',[0 0 0],'linewidth',2)
          %plot([0 9],[0 9*pc(1)+pc(2)],'k','linewidth',2)

          data=modc; time=time_mod;
          ythr = GPthres(time,data+tolerance,2,6,0.1,tthr,1); % Automate threshold selection
          disp(['GPD threshold selected: H_s = ' num2str(ythr) ' m']);
          figure;
          [GPD,gpp,gpci,Sy,ARI_POT,POT_sorted] = GPDfit(data+tolerance,time,tthr,ythr,ARI); % PF
          close
          figure(fign)
          semilogx(ARI_POT,POT_sorted,'o','color',colors{ke+1}),
          hold on
          semilogx(ARI,GPD(1,:),'-', 'color',colors{ke+1})
          semilogx(ARI,GPD(2,:),'--','color',colors{ke+1})
          semilogx(ARI,GPD(3,:),'--','color',colors{ke+1})
          legend('Obs.',expt,['corrected ',expt]);

        end
        
      end
    
    end
    
    if ~isempty(find(dcol==6))
      %subplot(length(dcol),1,find(dcol==6));
      %legend([legh],'location','best')
    end
    
    if save_fig==1
      %export_fig(gcf,'steephead_period','-png','-r150');
    end

  end % if plot_eva


  if plot_map

    scrsz=[1 1 1366 768];
    scrsz=[1 1 1920 1080];
    %scrsz=get(0,'screensize');
    
    figure('position',scrsz,'color',[1 1 1],'visible','on')
    hold on
    set(gca,'fontsize',12,'fontweight','bold')
      
    legh={'Obs'}; legt={'Obs'};

    checkin=0;  
    % Loop in desired time
    for t=time_lima
    
      ptime=datestr(t,'YYYY/mm/DD/HH/');
      ftime=datestr(t,'YYYYmmDDHH');
    
      ke=0;
      for expt=expts
        ke=ke+1;
      
        expt=expt{1};

        [ig,ltime]=grab_nc_sufix(expt);
 
      
        path_expt=[path_source,expt,'/']; % GLOBALWAVE/'];%2018/01/05/00'];
    
        pname=[path_expt,ptime];
        
        fname=[pname,ww3pre,'_',ftime,'-utc_',gnames{ig},'.nc'];
        %fname=[pname,'ww3g_',ftime,'-utc_',gnames{1},'.nc'];
        display(['Loading: ',fname]);
    
        if checkin==0; % t==time_lima(1)
          checkin=1; 
          if strcmp(expt,'NZWAVE-HR')
            depth_mod=double(ncread(fname,'depth',[1 1 1],[Inf Inf 1]));
            %wlv_mod=double(ncread(fname,'wlv',[1 1 1],[Inf Inf 1]));
            %depth_mod=depth_mod-wlv_mod;
          end

          lon_mod=double(ncread(fname,'lon'));
          lat_mod=double(ncread(fname,'lat'));
          [dif ilon]=nanmin(abs(lon_mod-lon_obs));
          [dif ilat]=nanmin(abs(lat_mod-lat_obs));

          [lon_modm,lat_modm]=meshgrid(lon_mod,lat_mod);
 
          if strcmp(file,'SteepHead_SP_FullRecord_QC')
            %ilon=ilon+3;
          elseif strcmp(file,'wairewa_lake_forsyth') & ig==2 % NZWAVE
            ilat=ilat-2;
          end
    
          dis=sqrt((lon_mod(ilon)-lon_obs).^2 + (lat_mod(ilat)-lat_obs).^2)*110;
          display(['Distance between obs and grid point is: ',num2str(dis),' km'])
        end
        
        display(['Reading nc variables']); tic
        time=squeeze(double(ncread(fname,'time')))./24+t; time=time(1:ltime-1);
        hs=squeeze(double(ncread(fname,'hsig',          [1 1 1],[Inf Inf ltime-1]))); 
        tp=squeeze(double(ncread(fname,'tpeak',         [1 1 1],[Inf Inf ltime-1]))); 
        %tm01=squeeze(double(ncread(fname,'tmean01',     [1 1 1],[Inf Inf ltime-1]))); 
        %tm=squeeze(double(ncread(fname,'tmean',         [1 1 1],[Inf Inf ltime-1]))); 
        pd=squeeze(double(ncread(fname,'peak_direction',[1 1 1],[Inf Inf ltime-1]))); 
        if strcmp(expt,'NZWAVE-HR-bla')
          wlv=squeeze(double(ncread(fname,'wlv',          [1 1 1],[Inf Inf ltime-1]))); 
          model(ke).wlv=wlv;
          ucur=squeeze(double(ncread(fname,'ucur',        [1 1 1],[Inf Inf ltime-1]))); 
          model(ke).ucur=ucur;
          vcur=squeeze(double(ncread(fname,'vcur',        [1 1 1],[Inf Inf ltime-1]))); 
          model(ke).vcur=vcur;
        else
          model(ke).wlv=nan(size(hs));
          model(ke).ucur=nan(size(hs));
          model(ke).vcur=nan(size(hs));
        end
        display(['Finished reading nc variables']); toc

        model(ke).time=time;
        model(ke).hs=hs;
        model(ke).tp=tp;
        model(ke).pd=pd;
        %model(ke).=;

      end % for expt=expts

      time=model(1).time;
      tnames={'Peak Period',    'Peak Direction',    'Directional spread',   'Tm01','Tm02 or Tm','Significant Wave Height (Hs)','Qp'};

      if plot_obs
        % converting to u and v before interpolating
        pd_obs=obs(:,2); hs_obs=obs(:,6);
        pd=mod(-90-pd_obs,360); u_obs=cosd(pd).*hs_obs; v_obs=sind(pd).*hs_obs;
        u_obsl=interp1(time_obs,u_obs,time); v_obsl=interp1(time_obs,v_obs,time);
        u_obsc=interp1(time_obs,u_obs,time,'cubic'); v_obsc=interp1(time_obs,v_obs,time,'cubic');
        dp_obs=270-atan2(v_obsc,u_obsc)*180/pi;
        dp_obs(dp_obs>360)=dp_obs(dp_obs>360)-360;

        pd_obs=dp_obs;
        hs_obs=sqrt(u_obsc.^2+v_obsc.^2);
        %hs_obs=interp1(time_obs,obs(:,6),time);
        tp_obs=interp1(time_obs,obs(:,1),time);
      end

      for i=1:length(time)
 
        ke=0;
        for expt=expts
          ke=ke+1;

          display(['Ploting map']); tic
          ax(ke)=subplot(1,length(expts),ke);
          hold on
          m_proj('lambert','long', [lon_obs-.3 lon_obs+.25],'lat',[lat_obs-.25 lat_obs+.25]);
          %m_proj('lambert','long', [nanmin(lon_mod(:)) nanmax(lon_mod(:))],'lat',[nanmin(lat_mod(:)) nanmax(lat_mod(:))]);

          if strcmp(pcname,'hs')
            data=model(ke).hs(:,:,i);
            cmap=cmocean('thermal');
            cmap=avhrr; % cmocean('thermal');
            caxis([-2 2])
          elseif strcmp(pcname,'wlv')
            data=model(ke).wlv(:,:,i);
            pd_mod=model(ke).pd(:,:,i);
            pd=mod(-90-pd_mod,360);
            u_mod=cosd(pd).*data; v_mod=sind(pd).*data;
            cmap=cmocean('balance');
            caxis([-2 2])
          elseif strcmp(pcname,'curr')
            u_mod=model(ke).ucur(:,:,i);
            v_mod=model(ke).vcur(:,:,i);
            data=sqrt(u_mod.^2+v_mod.^2);
            cmap=cmocean('tempo');
          end
          %colormap(ax(1),gwbreal)
          %m_pcolor(lon_ccmp,lat_ccmp,wsc_ccmp');
          %set(get(cs,'ylabel'),'string','N/m^3','fontsize',20,'fontweight','bold');
          %caxis([-2E-06 2E-06])
          colormap(ax(ke),cmap) % parula, avhrr
          m_pcolor(lon_mod,lat_mod,data');
          %caxis([0 .13])
          cb=colorbar;%('southoutside');
          set(get(cb,'ylabel'),'string','','fontsize',12,'fontweight','bold');
          set(cb,'fontsize',12,'fontweight','bold');

          if strcmp(pcname,'hs')
            caxis([0 4])
          elseif strcmp(pcname,'wlv')
            caxis([-2 2])
          elseif strcmp(pcname,'curr')
            caxis([0 2])
          end

          %%set(cs,'position','south');
          %shading interp
          %%caxis([nanmin(t_avhrr(:)) nanmax(t_avhrr(:))]); % 12 25

          %[cs,h]=m_contour(lon_bath,lat_bath,bath',[-1000 -1000],'color',[1 1 1],'linewidth',2);
          %clabel(cs,h,'color',[1 1 1],'fontsize',14,'fontweight','bold','LabelSpacing',2000)
          if strncmp(expt,'NZWAVE',6)
            [cs,h]=m_contour(lon_mod,lat_mod,-depth_mod',[-2 -20 -50 -100 -200],'color',[1 1 1],'linewidth',2);
            clabel(cs,h,'color',[1 1 1],'fontsize',10,'fontweight','bold','LabelSpacing',2000)
          else
            [cs,h]=m_contour(lon_geb,lat_geb,-depth_geb',[-200 -1000],'color',[1 1 1],'linewidth',2);
            clabel(cs,h,'color',[0 0 0],'fontsize',10,'fontweight','bold','LabelSpacing',2000)
          end

          %[cs,h]=m_contour(lon_aviso,lat_aviso,z_aviso_avg,'color','k');
          %[cs,h]=m_contour(lon_aviso,lat_aviso,z_aviso_avg,[.76:.04:.88],'color','k');
          %clabel(cs,h,'color',[0 0 0],'fontsize',14,'fontweight','bold','LabelSpacing',2000) % each contour line represents 0.02 increase/decrease

          m_gshhs_i('patch',[.7 .7 .7],'edgecolor','k','linewidth',.5);
          %%m_grid('fontname','helvetica','fontsize',22,'fontweight','bold');
          %%m_grid('xtick',[round(loni):2:round(lonf)],'xticklabels',[round(loni):2:round(lonf)],'fontname','helvetica','fontsize',20,'fontweight','bold');
          m_grid('fontname','helvetica','fontsize',14,'fontweight','bold');
          title([replace(file,'_',' '),' Hs and ',pcname,' at ',datestr(time(i),'HH:MM dd/mmm/yyyy')],'fontsize',14,'fontweight','bold')

          %set(gca,'fontsize',14,'fontname','helvetica','fontweight','bold')
          %%m_contour(plon,plat,bathy',[200 2000],'w','linewidth',1);
          %%caxis([0 .5])

          % getting u and v
          if strcmp(vecname,'wlv')
          elseif strcmp(vecname,'hs')
            dat=model(ke).hs(:,:,i);
            pd_mod=model(ke).pd(:,:,i);
            pd=mod(-90-pd_mod,360);
            u_mod=cosd(pd).*dat; v_mod=sind(pd).*dat;
          elseif strcmp(vecname,'curr')
            u_mod=model(ke).ucur(:,:,i);
            v_mod=model(ke).vcur(:,:,i);
          end

          v_spa=1; vsize=2;
          if strcmp(vecname,'hs')
            [hp ht]=m_vec(vsize,lon_modm(1:v_spa:end,1:v_spa:end),lat_modm(1:v_spa:end,1:v_spa:end),u_mod(1:v_spa:end,1:v_spa:end)',v_mod(1:v_spa:end,1:v_spa:end)','k','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5);
            if plot_obs
              [hp ht]=m_vec(vsize,lon_obs+.1,lat_obs+.14,2,0,'k','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5,'key',[num2str(2.0),' m']);
            end
            set(ht,'fontsize',14,'fontweight','bold')
          elseif strcmp(vecname,'curr')
            [hp ht]=m_vec(5,lon_modm(1:v_spa:end,1:v_spa:end),lat_modm(1:v_spa:end,1:v_spa:end),u_mod(1:v_spa:end,1:v_spa:end)',v_mod(1:v_spa:end,1:v_spa:end)','k','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5)
            if plot_obs
              [hp ht]=m_vec(5,lon_obs+.1,lat_obs+.14,2,0,'k','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5,'key',[num2str(2.0),' m/s']);
            end
            set(ht,'fontsize',14,'fontweight','bold')
          end

          if plot_obs
            pd_mod=model(ke).pd(:,:,i);
            m_text(lon_obs+.1,lat_obs+.2,['Mod dp: ',num2str(pd_mod(ilon,ilat),'%.1f'),'^o'],'fontsize',14,'fontweight','bold') 
            m_text(lon_obs+.1,lat_obs+.18,['Obs dp: ',num2str(pd_obs(i),'%.2f'),'^o'],'fontsize',14,'fontweight','bold') 
          end

          % model single vector 
          data=model(ke).hs(:,:,i);
          pd_mod=model(ke).pd(:,:,i);
          pd=mod(-90-pd_mod,360);
          u_mod=cosd(pd).*data; v_mod=sind(pd).*data;
      	  vsize=2;
          m_vec(vsize,lon_mod(ilon),lat_mod(ilat),u_mod(ilon,ilat),v_mod(ilon,ilat),'m','shaftwidth',2,'headwidth',6,'headangle',60,'headlength',8)
          [hp ht]=m_vec(vsize,lon_obs+.1,lat_obs+.12,2,0                           ,'m','shaftwidth',2,'headwidth',6,'headangle',60,'headlength',8,'key',[num2str(2.0),' m']);
          set(ht,'fontsize',14,'fontweight','bold')
          
          % obs single vector 
          if plot_obs
            pd=mod(-90-pd_obs,360);
            u_obs=cosd(pd).*hs_obs; v_obs=sind(pd).*hs_obs;
            m_vec(vsize,lon_obs,lat_obs,u_obsc(i),v_obsc(i),'m'                          ,'shaftwidth',2,'headwidth',6,'headangle',60,'headlength',8)
            %m_vec(10,lon_obs,lat_obs,u_obsl(i),v_obsl(i),'c','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5)
					end	

          display(['Finished ploting map']); toc

          %if dcol(i)==6
          %  mod_rmse=nanrmse(obs(:,dcol(i)),modeli');
          %  mod_bias=nanmean(modeli')-nanmean(obs(:,dcol(i)));
          %  legh=[legh,{[expt,' rmse=',num2str(mod_rmse,'%.2f'),', bias=',num2str(mod_bias,'%.2f')]}];
          %end

        end % for expt=expts
				pause(0.5)

        path_dm=strcat(path_fig,expt{1},'/');
        system(['mkdir -p ',path_dm]);
        figname=[path_dm,file,'_',pcname,'_',vecname,'_',datestr(time(i),'YYYY_mm_DD_HH_MM'),'.png'];
        display(['Plotting: ',figname]);
        if save_fig==1
          display(['Saving: ',figname]);
          export_fig(gcf,figname,'-png','-r150');
          %close
          %clf('reset')	
				end

        if save_video==1 

          if i==1 & t==time_lima(1)
            videoName = [path_dm,'/',file,'_',pcname,'_',vecname,'_',datestr(time_lima(1),'yyyymmdd_HHMMSS'),'_',datestr(time_lima(end),'yyyymmdd_HHMMSS'),''];
            display(['Saving video name: ',videoName])
            vidObj = VideoWriter(videoName,'Uncompressed AVI');
            %vidObj = VideoWriter(videoName,'Motion JPEG 2000');
            set(vidObj,'FrameRate',5)
            open(vidObj);
          end
          display(['Printing frame time: ',datestr(time(i),'yyyymmdd_HHMMSS')]);
          pause(0.5)
          M = getframe(gcf);
          writeVideo(vidObj,M);

        end

        if save_gif==1 
          gifName = [path_dm,'/',file,'_',pcname,'_',vecname,'_',datestr(time_lima(1),'yyyymmdd_HHMMSS'),'_',datestr(time_lima(end),'yyyymmdd_HHMMSS'),'.gif'];
          display(['Saving GIF: ',gifName]);
          drawnow
          frame = getframe(1);
          im = frame2im(frame);
          [imind,cm] = rgb2ind(im,256);
          if i == 1 & t==time_lima(1);
            imwrite(imind,cm,gifName,'gif', 'Loopcount',inf,'DelayTime',.2);
          else
            imwrite(imind,cm,gifName,'gif','WriteMode','append','DelayTime',.2);
          end
				end

        if save_fig==1 || save_video==1 || save_gif==1
          clf('reset')
          set(gcf,'color',[1 1 1])
				else
				  pause(0.5)
				end

      end % for i=1:length(time)

    end % for t=time_lima
    
    if save_video==1 
      close(vidObj);
    end

    %close 

  end % if plot_map

  if plot_matm

    scrsz=[1 1 1366 768];
    scrsz=[1 1 1920 1080];
    %scrsz=get(0,'screensize');
    
    figure('position',scrsz,'color',[1 1 1],'visible','on')
    hold on
    set(gca,'fontsize',12,'fontweight','bold')
      
    legh={'Obs'}; legt={'Obs'};

    checkin=0;  
    % Loop in desired time
    for t=time_lima
    
      ptime=datestr(t,'YYYY/mm/DD/HH/');
      ftime=datestr(t,'YYYYmmDDHH');
    
      ke=0;
      for expt=expts
        ke=ke+1;
      
        expt=expt{1};

        [ig,ltime]=grab_nc_sufix(expt);
      
        path_expt=[path_source,expt,'/']; % GLOBALWAVE/'];%2018/01/05/00'];
    
        pname=[path_expt,ptime];
        
        fname=[pname,ww3pre,'_',ftime,'-utc_',gnames{ig},'.nc'];
        %fname=[pname,'ww3g_',ftime,'-utc_',gnames{1},'.nc'];
        display(['Loading: ',fname]);
    
        if checkin==0; % t==time_lima(1)
          checkin=1; 
          if strcmp(expt,'NZWAVE-HR')
            depth_mod=double(ncread(fname,'depth',[1 1 1],[Inf Inf 1]));
            wlv_mod=double(ncread(fname,'wlv',[1 1 1],[Inf Inf 1]));
            depth_mod=depth_mod-wlv_mod;
          end

          lon_mod=double(ncread(fname,'lon'));
          lat_mod=double(ncread(fname,'lat'));
          [dif ilon]=nanmin(abs(lon_mod-lon_obs));
          [dif ilat]=nanmin(abs(lat_mod-lat_obs));

          [lon_modm,lat_modm]=meshgrid(lon_mod,lat_mod);
 
          if strcmp(file,'SteepHead_SP_FullRecord_QC')
            %ilon=ilon+3;
          elseif strcmp(file,'wairewa_lake_forsyth') & ig==2 % NZWAVE
            ilat=ilat-2;
          end
    
          dis=sqrt((lon_mod(ilon)-lon_obs).^2 + (lat_mod(ilat)-lat_obs).^2)*110;
          display(['Distance between obs and grid point is: ',num2str(dis),' km'])
        end
        
        display(['Reading nc variables']); tic
        time=squeeze(double(ncread(fname,'time')))./24+t; time=time(1:ltime-1);
        hs=squeeze(double(ncread(fname,'hsig',          [1 1 1],[Inf Inf ltime-1]))); 
        tp=squeeze(double(ncread(fname,'tpeak',         [1 1 1],[Inf Inf ltime-1]))); 
        %tm01=squeeze(double(ncread(fname,'tmean01',     [1 1 1],[Inf Inf ltime-1]))); 
        %tm=squeeze(double(ncread(fname,'tmean',         [1 1 1],[Inf Inf ltime-1]))); 
        pd=squeeze(double(ncread(fname,'peak_direction',[1 1 1],[Inf Inf ltime-1]))); 
        if strcmp(expt,'NZWAVE-HR')
          wlv=squeeze(double(ncread(fname,'wlv',          [1 1 1],[Inf Inf ltime-1]))); 
          model(ke).wlv=wlv;
          ucur=squeeze(double(ncread(fname,'ucur',        [1 1 1],[Inf Inf ltime-1]))); 
          model(ke).ucur=ucur;
          vcur=squeeze(double(ncread(fname,'vcur',        [1 1 1],[Inf Inf ltime-1]))); 
          model(ke).vcur=vcur;
        else
          model(ke).wlv=nan(size(hs));
          model(ke).ucur=nan(size(hs));
          model(ke).vcur=nan(size(hs));
        end
        display(['Finished reading nc variables']); toc

        model(ke).time=time;
        model(ke).hs=hs;
        model(ke).tp=tp;
        model(ke).pd=pd;
        %model(ke).=;

      end % for expt=expts

      time=model(1).time;
      tnames={'Peak Period',    'Peak Direction',    'Directional spread',   'Tm01','Tm02 or Tm','Significant Wave Height (Hs)','Qp'};

      if plot_obs
        % converting to u and v before interpolating
        pd_obs=obs(:,2); hs_obs=obs(:,6);
        pd=mod(-90-pd_obs,360); u_obs=cosd(pd).*hs_obs; v_obs=sind(pd).*hs_obs;
        u_obsl=interp1(time_obs,u_obs,time); v_obsl=interp1(time_obs,v_obs,time);
        u_obsc=interp1(time_obs,u_obs,time,'cubic'); v_obsc=interp1(time_obs,v_obs,time,'cubic');
        dp_obs=270-atan2(v_obsc,u_obsc)*180/pi;
        dp_obs(dp_obs>360)=dp_obs(dp_obs>360)-360;

        pd_obs=dp_obs;
        hs_obs=sqrt(u_obsc.^2+v_obsc.^2);
        %hs_obs=interp1(time_obs,obs(:,6),time);
        tp_obs=interp1(time_obs,obs(:,1),time);
      end

      for i=1:6:length(time)

        ptime=datestr(time(i),'YYYY/');
        ftime=datestr(time(i),'YYYYmmDD');

        path_sp=[path_ssp,'psl/',ptime];
        %file_sp=[path_sp,prefix{1},midfix{1},sufix{1},'_',ftime,'00.nc'];
        if time_lima(1)<datenum('2015-01-01 00:00:00')
      	  file_sp=[path_sp,'psl_historical_GFDL-ESM4_CCAM_6hrly_Global_raw_',ftime,'00.nc'];
          time_sp=ncread(file_sp,'time')./(24*60) + datenum('1959-01-01 00:00:00');
      	else
      	  file_sp=[path_sp,'psl_ssp370_GFDL-ESM4_CCAM_6hrly_Global_raw_',ftime,'00.nc'];
          time_sp=ncread(file_sp,'time')./(24*60) + datenum('2015-01-01 00:00:00');
      	end
        psl=ncread(file_sp,'psl')./100;
        
        [dif i_sp]=nanmin(abs(time_sp-time(i))); 
        psl=psl(:,:,i_sp);
 
        ke=0;
        for expt=expts
          ke=ke+1;

          display(['Ploting map']); tic
          ax(ke)=subplot(1,length(expts),ke);
          hold on
          %m_proj('lambert','long', [lon_obs-.3 lon_obs+.25],'lat',[lat_obs-.25 lat_obs+.25]);
          m_proj('lambert','long', [nanmin(lon_mod(:)) nanmax(lon_mod(:))],'lat',[nanmin(lat_mod(:)) nanmax(lat_mod(:))]);

          if strcmp(pcname,'hs')
            data=model(ke).hs(:,:,i);
            cmap=cmocean('thermal');
            cmin=0; cint=.5; cmax=7;
            cmap=mycmap(cmap,cmin,cint,cmax);

          elseif strcmp(pcname,'wlv')
            data=model(ke).wlv(:,:,i);
            pd_mod=model(ke).pd(:,:,i);
            pd=mod(-90-pd_mod,360);
            u_mod=cosd(pd).*data; v_mod=sind(pd).*data;
            cmap=cmocean('balance');
            caxis([-2 2])
          elseif strcmp(pcname,'curr')
            u_mod=model(ke).ucur(:,:,i);
            v_mod=model(ke).vcur(:,:,i);
            data=sqrt(u_mod.^2+v_mod.^2);
            cmap=cmocean('tempo');
          end
          %colormap(ax(1),gwbreal)
          %m_pcolor(lon_ccmp,lat_ccmp,wsc_ccmp');
          %set(get(cs,'ylabel'),'string','N/m^3','fontsize',20,'fontweight','bold');
          %caxis([-2E-06 2E-06])
          colormap(ax(ke),cmap) % parula, avhrr
          m_pcolor(lon_mod,lat_mod,data');
          %caxis([0 .13])
          cb=colorbar;%('southoutside');
          set(get(cb,'ylabel'),'string','','fontsize',12,'fontweight','bold');
          set(cb,'fontsize',12,'fontweight','bold');

          if strcmp(pcname,'hs')
          %  caxis([0 4])
          elseif strcmp(pcname,'wlv')
            caxis([-2 2])
          elseif strcmp(pcname,'curr')
            caxis([0 2])
          end
          colormap(cmap)
          caxis([cmin cmax+cint]);

          if strncmp(conname,'psl',3)
            [cs,h]=m_contour(lon_sspm,lat_sspm,psl','color',[1 1 1],'linewidth',2);
            clabel(cs,h,'color',[1 1 1],'fontsize',10,'fontweight','bold','LabelSpacing',2000)
          end

          plot_track=1;
          if plot_track
            [dif i_tr]=nanmin(abs(time_tr-time(i))); 
            if time_tr(i_tr)==time(i); 
              loc=find(id_tr(:)==id_tr(i_tr));
              m_plot(tr.data(loc(1):i_tr,1),tr.data(loc(1):i_tr,2),'m','linewidth',2)
              m_plot(tr.data(i_tr,1),tr.data(i_tr,2),'.m','markersize',20)
            end
          end


          %%set(cs,'position','south');
          %shading interp
          %%caxis([nanmin(t_avhrr(:)) nanmax(t_avhrr(:))]); % 12 25

          %[cs,h]=m_contour(lon_bath,lat_bath,bath',[-1000 -1000],'color',[1 1 1],'linewidth',2);
          %clabel(cs,h,'color',[1 1 1],'fontsize',14,'fontweight','bold','LabelSpacing',2000)
          if strncmp(expt,'NZWAVE',6)
            %[cs,h]=m_contour(lon_mod,lat_mod,-depth_mod',[-200 -200],'color',[0 .7 0],'linewidth',2);
            %clabel(cs,h,'color',[0 .7 0],'fontsize',10,'fontweight','bold','LabelSpacing',2000)
          else
            [cs,h]=m_contour(lon_geb,lat_geb,-depth_geb',[-200 -1000],'color',[0 0 0],'linewidth',2);
            clabel(cs,h,'color',[0 0 0],'fontsize',10,'fontweight','bold','LabelSpacing',2000)
          end

          %[cs,h]=m_contour(lon_aviso,lat_aviso,z_aviso_avg,'color','k');
          %[cs,h]=m_contour(lon_aviso,lat_aviso,z_aviso_avg,[.76:.04:.88],'color','k');
          %clabel(cs,h,'color',[0 0 0],'fontsize',14,'fontweight','bold','LabelSpacing',2000) % each contour line represents 0.02 increase/decrease

          m_gshhs_i('patch',[.7 .7 .7],'edgecolor',[0 0 0],'linewidth',.5);
          %%m_grid('fontname','helvetica','fontsize',22,'fontweight','bold');
          %%m_grid('xtick',[round(loni):2:round(lonf)],'xticklabels',[round(loni):2:round(lonf)],'fontname','helvetica','fontsize',20,'fontweight','bold');
          m_grid('fontname','helvetica','fontsize',14,'fontweight','bold');
          title([replace(file,'_',' '),' Hs and ',pcname,' at ',datestr(time(i),'HH:MM dd/mmm/yyyy')],'fontsize',14,'fontweight','bold')
          title(['Hs and ',conname,' at ',datestr(time(i),'HH:MM dd/mmm/yyyy')],'fontsize',14,'fontweight','bold')

          %set(gca,'fontsize',14,'fontname','helvetica','fontweight','bold')
          %%m_contour(plon,plat,bathy',[200 2000],'w','linewidth',1);
          %%caxis([0 .5])

          % getting u and v
          if strcmp(vecname,'wlv')
          elseif strcmp(vecname,'hs')
            dat=model(ke).hs(:,:,i);
            pd_mod=model(ke).pd(:,:,i);
            pd=mod(-90-pd_mod,360);
            u_mod=cosd(pd).*dat; v_mod=sind(pd).*dat;
          elseif strcmp(vecname,'curr')
            u_mod=model(ke).ucur(:,:,i);
            v_mod=model(ke).vcur(:,:,i);
          end

          v_spa=40;
          if strcmp(vecname,'hs')
            [hp ht]=m_vec(10,lon_modm(1:v_spa:end,1:v_spa:end),lat_modm(1:v_spa:end,1:v_spa:end),u_mod(1:v_spa:end,1:v_spa:end)',v_mod(1:v_spa:end,1:v_spa:end)','k','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5);
            if plot_obs
              [hp ht]=m_vec(10,lon_obs+.1,lat_obs+.14,2,0,'k','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5,'key',[num2str(2.0),' m']);
            end
            set(ht,'fontsize',14,'fontweight','bold')
          elseif strcmp(vecname,'curr')
            [hp ht]=m_vec(5,lon_modm(1:v_spa:end,1:v_spa:end),lat_modm(1:v_spa:end,1:v_spa:end),u_mod(1:v_spa:end,1:v_spa:end)',v_mod(1:v_spa:end,1:v_spa:end)','k','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5)
            if plot_obs
              [hp ht]=m_vec(5,lon_obs+.1,lat_obs+.14,2,0,'k','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5,'key',[num2str(2.0),' m/s']);
            end
            set(ht,'fontsize',14,'fontweight','bold')
          end

          if plot_obs
            pd_mod=model(ke).pd(:,:,i);
            m_text(lon_obs+.1,lat_obs+.2,['Mod dp: ',num2str(pd_mod(ilon,ilat),'%.1f'),'^o'],'fontsize',14,'fontweight','bold') 
            m_text(lon_obs+.1,lat_obs+.18,['Obs dp: ',num2str(pd_obs(i),'%.2f'),'^o'],'fontsize',14,'fontweight','bold') 
          end

          if plot_obs
            % model single vector 
            data=model(ke).hs(:,:,i);
            pd_mod=model(ke).pd(:,:,i);
            pd=mod(-90-pd_mod,360);
            u_mod=cosd(pd).*data; v_mod=sind(pd).*data;
      	    vsize=2;
            m_vec(vsize,lon_mod(ilon),lat_mod(ilat),u_mod(ilon,ilat),v_mod(ilon,ilat),'b','shaftwidth',2,'headwidth',6,'headangle',60,'headlength',8)
            [hp ht]=m_vec(vsize,lon_obs+.1,lat_obs+.12,2,0,'b',                           'shaftwidth',2,'headwidth',6,'headangle',60,'headlength',8,'key',[num2str(2.0),' m']);
            set(ht,'fontsize',14,'fontweight','bold')
            
            % obs single vector 
            %pd=mod(-90-pd_obs,360);
            %u_obs=cosd(pd).*hs_obs; v_obs=sind(pd).*hs_obs;
            m_vec(vsize,lon_obs,lat_obs,u_obsc(i),v_obsc(i),'m'                          ,'shaftwidth',2,'headwidth',6,'headangle',60,'headlength',8)
            %m_vec(10,lon_obs,lat_obs,u_obsl(i),v_obsl(i),'c','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5)
					end

          display(['Finished ploting map']); toc

          %if dcol(i)==6
          %  mod_rmse=nanrmse(obs(:,dcol(i)),modeli');
          %  mod_bias=nanmean(modeli')-nanmean(obs(:,dcol(i)));
          %  legh=[legh,{[expt,' rmse=',num2str(mod_rmse,'%.2f'),', bias=',num2str(mod_bias,'%.2f')]}];
          %end

        end % for expt=expts
				pause(0.5)

        path_dm=strcat(path_fig,expt{1},'/');
        system(['mkdir -p ',path_dm]);
        figname=[path_dm,file,'_',pcname,'_',vecname,'_',datestr(time(i),'YYYY_mm_DD_HH_MM'),'.png'];
        display(['Plotting: ',figname]);
        if save_fig==1
          display(['Saving: ',figname]);
          export_fig(gcf,figname,'-png','-r150');
          %close
          %clf('reset')	
				end

        if save_video==1 

          if i==1 & t==time_lima(1)
            videoName = [path_dm,'/',file,'_',pcname,'_',vecname,'_',datestr(time_lima(1),'yyyymmdd_HHMMSS'),'_',datestr(time_lima(end),'yyyymmdd_HHMMSS'),''];
            display(['Saving video name: ',videoName])
            vidObj = VideoWriter(videoName,'Uncompressed AVI');
            %vidObj = VideoWriter(videoName,'Motion JPEG 2000');
            set(vidObj,'FrameRate',5)
            open(vidObj);
          end
          display(['Printing frame time: ',datestr(time(i),'yyyymmdd_HHMMSS')]);
          pause(0.5)
          M = getframe(gcf);
          writeVideo(vidObj,M);

        end

        if save_gif==1 
          gifName = [path_dm,'/',file,'_',pcname,'_',vecname,'_',datestr(time_lima(1),'yyyymmdd_HHMMSS'),'_',datestr(time_lima(end),'yyyymmdd_HHMMSS'),'.gif'];
          display(['Saving GIF: ',gifName]);
          drawnow
          frame = getframe(1);
          im = frame2im(frame);
          [imind,cm] = rgb2ind(im,256);
          if i == 1 & t==time_lima(1);
            imwrite(imind,cm,gifName,'gif', 'Loopcount',inf,'DelayTime',.5);
          else
            imwrite(imind,cm,gifName,'gif','WriteMode','append','DelayTime',.5);
          end
				end

        if save_fig==1 || save_video==1 || save_gif==1
          clf('reset')
          set(gcf,'color',[1 1 1])
				else
				  pause(0.5)
				end

      end % for i=1:length(time)

    end % for t=time_lima
    
    if save_video==1 
      close(vidObj);
    end

    %close 

  end % if plot_matm


  if plot_downs

    lon_obss=[lon_obss;lon_obs];
    lat_obss=[lat_obss;lat_obs];

    if fobs==stations(1)

      scrsz=[1 1 1366 768];
      scrsz=[1 1 1920 1080];
      %scrsz=get(0,'screensize');
      
      figure('position',scrsz,'color',[1 1 1],'visible','on')
      hold on
      set(gca,'fontsize',12,'fontweight','bold')

    elseif fobs==stations(end)
        
      legh={'Obs'}; legt={'Obs'};

      checkin=0;  
      % Loop in desired time
      for t=time_lima%(1)
       
        ptime=datestr(t,'YYYY/mm/DD/00/');
        ftime=datestr(t,'YYYYmmDD00');
      
        ax(1)=subplot(1,1,1);
        ke=0;

        for expt=expts
          ke=ke+1;
        
          expt=expt{1};

          [ig,ltime]=grab_nc_sufix(expt);
        
          path_expt=[path_source,expt,'/']; % GLOBALWAVE/'];%2018/01/05/00'];
      
          pname=[path_expt,ptime];
          
          fname=[pname,ww3pre,'_',ftime,'-utc_',gnames{ig},'.nc'];
          %fname=[pname,'ww3g_',ftime,'-utc_',gnames{1},'.nc'];
          display(['Loading: ',fname]);

          % getting time
          time=squeeze(double(ncread(fname,'time')))./24+floor(t); time=time(1:ltime-1);
          % finding the right hour
          [dif i]=nanmin(abs(time-t)); %:length(time)

          if strcmp(expt,'NZWAVE-HR')
            %depth_mod=double(ncread(fname,'depth',[1 1 1],[Inf Inf 1]));
            %wlv_mod=double(ncread(fname,'wlv',[1 1 1],[Inf Inf 1]));
            %depth_mod=depth_mod-wlv_mod;
          end

          lon_mod=double(ncread(fname,'lon'));
          lat_mod=double(ncread(fname,'lat'));
          [dif ilon]=nanmin(abs(lon_mod-lon_obs));
          [dif ilat]=nanmin(abs(lat_mod-lat_obs));

          [lon_modm,lat_modm]=meshgrid(lon_mod,lat_mod);
 
          if strcmp(file,'SteepHead_SP_FullRecord_QC')
            %ilon=ilon+3;
          elseif strcmp(file,'wairewa_lake_forsyth') & ig==2 % NZWAVE
            ilat=ilat-2;
          end
      
          dis=sqrt((lon_mod(ilon)-lon_obs).^2 + (lat_mod(ilat)-lat_obs).^2)*110;
          display(['Distance between obs and grid point is: ',num2str(dis),' km'])
          
          display(['Reading nc variables']); tic
          hs=squeeze(double(ncread(fname,'hsig',          [1 1 i],[Inf Inf 1]))); 
          tp=squeeze(double(ncread(fname,'tpeak',         [1 1 i],[Inf Inf 1]))); 
          %tm01=squeeze(double(ncread(fname,'tmean01',     [1 1 1],[Inf Inf ltime-1]))); 
          %tm=squeeze(double(ncread(fname,'tmean',         [1 1 1],[Inf Inf ltime-1]))); 
          pd=squeeze(double(ncread(fname,'peak_direction',[1 1 i],[Inf Inf 1]))); 
          if strcmp(expt,'NZWAVE-HR')
            wlv=squeeze(double(ncread(fname,'wlv',          [1 1 i],[Inf Inf 1]))); 
            model(ke).wlv=wlv;
            ucur=squeeze(double(ncread(fname,'ucur',        [1 1 i],[Inf Inf 1]))); 
            model(ke).ucur=ucur;
            vcur=squeeze(double(ncread(fname,'vcur',        [1 1 i],[Inf Inf 1]))); 
            model(ke).vcur=vcur;
          else
            model(ke).wlv=nan(size(hs));
            model(ke).ucur=nan(size(hs));
            model(ke).vcur=nan(size(hs));
          end
          display(['Finished reading nc variables']); toc

          model(ke).time=time;
          model(ke).hs=hs;
          model(ke).tp=tp;
          model(ke).pd=pd;
          %model(ke).=;

          time=model(ke).time;
          tnames={'Peak Period',    'Peak Direction',    'Directional spread',   'Tm01','Tm02 or Tm','Significant Wave Height (Hs)','Qp'};

          % converting to u and v before interpolating
          pd_obs=obs(:,2); hs_obs=obs(:,6);
          pd=mod(-90-pd_obs,360); u_obs=cosd(pd).*hs_obs; v_obs=sind(pd).*hs_obs;
          u_obsl=interp1(time_obs,u_obs,time); v_obsl=interp1(time_obs,v_obs,time);
          u_obsc=interp1(time_obs,u_obs,time,'cubic'); v_obsc=interp1(time_obs,v_obs,time,'cubic');
          dp_obs=270-atan2(v_obsc,u_obsc)*180/pi;
          dp_obs(dp_obs>360)=dp_obs(dp_obs>360)-360;

          pd_obs=dp_obs;
          hs_obs=sqrt(u_obsc.^2+v_obsc.^2);
          %hs_obs=interp1(time_obs,obs(:,6),time);
          tp_obs=interp1(time_obs,obs(:,1),time);
            
          %for i=1%:length(time)
 
            %ke=0;
            %for expt=expts
            %  ke=ke+1;

          display(['Ploting map']); tic
          if strcmp(pcname,'hs')
            data=model(ke).hs;%(:,:,i);
            cmap=avhrr; % cmocean('matter');
            %cmap=cmocean('ice');
            cmap = mycmap(cmap,1,1,14);
          elseif strcmp(pcname,'wlv')
            data=model(ke).wlv;%(:,:,i);
            pd_mod=model(ke).pd;%(:,:,i);
            pd=mod(-90-pd_mod,360);
            u_mod=cosd(pd).*data; v_mod=sind(pd).*data;
            cmap=cmocean('balance');
            caxis([-2 2])
          elseif strcmp(pcname,'curr')
            u_mod=model(ke).ucur;%(:,:,i);
            v_mod=model(ke).vcur;%(:,:,i);
            data=sqrt(u_mod.^2+v_mod.^2);
            cmap=cmocean('thermal');
          end

          hold on
          m_proj('lambert','long', [143.3203125-2 184.570495605469+2],'lat',[-54.453125-2 -20.8589458465576+2]);

          %colormap(ax(1),gwbreal)
          %m_pcolor(lon_ccmp,lat_ccmp,wsc_ccmp');
          %set(get(cs,'ylabel'),'string','N/m^3','fontsize',20,'fontweight','bold');
          %caxis([-2E-06 2E-06])
          colormap(ax(1),cmap) % parula, avhrr
          m_pcolor(lon_mod,lat_mod,data');
          %m_contourf(lon_mod,lat_mod,data',[0:.5:8]);
          %m_contour(lon_mod,lat_mod,data',[0:.5:8]);
          %caxis([0 .13])

          cb=colorbar;%('southoutside');
          set(get(cb,'ylabel'),'string','Sig. Wave Height (m)','fontsize',12,'fontweight','bold');
          set(cb,'fontsize',12,'fontweight','bold');

          if strcmp(pcname,'hs')
            caxis([0 14])
            %caxis([0 20])
          elseif strcmp(pcname,'wlv')
            caxis([-2 2])
          elseif strcmp(pcname,'curr')
            cmap=cmocean('thermal');
          end

          %%set(cs,'position','south');
          %shading interp
          %%caxis([nanmin(t_avhrr(:)) nanmax(t_avhrr(:))]); % 12 25

          %[cs,h]=m_contour(lon_bath,lat_bath,bath',[-1000 -1000],'color',[1 1 1],'linewidth',2);
          %clabel(cs,h,'color',[1 1 1],'fontsize',14,'fontweight','bold','LabelSpacing',2000)
          if ke==length(expts)
            %[cs,h]=m_contour(lon_mod,lat_mod,-depth_mod',[-2 -20 -50 -100 -200],'color',[1 1 1],'linewidth',2);
            %clabel(cs,h,'color',[0 0 0],'fontsize',10,'fontweight','bold','LabelSpacing',2000)
            [cs,h]=m_contour(lon_gebco,lat_gebco,bat_gebco',[-200 -200],'color',[146 36 40]./255,'linewidth',2);
            %clabel(cs,h,'color',[0 0 0],'fontsize',10,'fontweight','bold','LabelSpacing',2000)
          end

          %[cs,h]=m_contour(lon_aviso,lat_aviso,z_aviso_avg,'color','k');
          %[cs,h]=m_contour(lon_aviso,lat_aviso,z_aviso_avg,[.76:.04:.88],'color','k');
          %clabel(cs,h,'color',[0 0 0],'fontsize',14,'fontweight','bold','LabelSpacing',2000) % each contour line represents 0.02 increase/decrease

          if ke==1%length(expts)
            %%m_grid('fontname','helvetica','fontsize',22,'fontweight','bold');
            %%m_grid('xtick',[round(loni):2:round(lonf)],'xticklabels',[round(loni):2:round(lonf)],'fontname','helvetica','fontsize',20,'fontweight','bold');
          %elseif ke==length(expts)
            m_gshhs_i('patch',[.7 .7 .7],'edgecolor','k','linewidth',.5);
            m_grid('fontname','helvetica','fontsize',14,'fontweight','bold');
            %title([replace(file,'_',' '),' Hs and water level at ',datestr(time(i),'HH:MM dd/mmm/yyyy')],'fontsize',14,'fontweight','bold')
            title(['Sig. Wave Height at ',datestr(time(i),'HH:MM dd/mmm/yyyy'),' UTC'],'fontsize',14,'fontweight','bold')
          end
          %set(gca,'fontsize',14,'fontname','helvetica','fontweight','bold')
          %%m_contour(plon,plat,bathy',[200 2000],'w','linewidth',1);
          %%caxis([0 .5])

          % getting u and v
          if strcmp(vecname,'wlv')
          elseif strcmp(vecname,'hs')
            dat=model(ke).hs;%(:,:,i);
            pd_mod=model(ke).pd;%(:,:,i);
            pd=mod(-90-pd_mod,360);
            u_mod=cosd(pd).*dat; v_mod=sind(pd).*dat;
          elseif strcmp(vecname,'curr')
            u_mod=model(ke).ucur;%(:,:,i);
            v_mod=model(ke).vcur;%(:,:,i);
          end

          v_spa=40;
          if strcmp(vecname,'hs')
            m_vec(50,lon_modm(10:v_spa:end,10:v_spa:end),lat_modm(10:v_spa:end,10:v_spa:end),u_mod(10:v_spa:end,10:v_spa:end)',v_mod(10:v_spa:end,10:v_spa:end)','w','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5)
            if ke==1
              [hp ht]=m_vec(50,144,-30,8,0,'w','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5,'key',[num2str(8.0),' m']);
              set(ht,'fontsize',14,'fontweight','bold')
            end
          end
          if ke==1
            m_text(144,-29.5,['Sig. Wave Height'],'color','w','fontsize',14,'fontweight','bold') 
          end

          pd_mod=model(ke).pd;%(:,:,i);
          for o=1:length(stations)
            filo=files{stations(o)};
            filo=replace(filo,'_',' ');
            m_plot(lon_obss(o),lat_obss(o),'.w','markersize',30)
            m_plot(lon_obss(o),lat_obss(o),'.r','markersize',20)
            if o==1 
              %m_text(lon_obss(o)-6,lat_obss(o)+.2,filo,'color','k','fontsize',14,'fontweight','bold') 
            elseif o==2 
              %m_text(lon_obss(o)-8,lat_obss(o)+.2,filo,'color','k','fontsize',14,'fontweight','bold') 
            end
          end
          %m_text(lon_obs+.1,lat_obs+.18,['Obs dp: ',num2str(pd_obs(i),'%.2f'),'^o'],'fontsize',14,'fontweight','bold') 

          % model single vector 
          data=model(ke).hs;%(:,:,i);
          pd_mod=model(ke).pd;%(:,:,i);
          pd=mod(-90-pd_mod,360);
          u_mod=cosd(pd).*data; v_mod=sind(pd).*data;
          %m_vec(50,lon_mod(ilon),lat_mod(ilat),u_mod(ilon,ilat),v_mod(ilon,ilat),'b','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5)
          %[hp ht]=m_vec(50,lon_obs+.1,lat_obs+.12,2,0,'b','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5,'key',[num2str(2.0),' m']);
          %set(ht,'fontsize',14,'fontweight','bold')
          
          %% obs single vector 
          %%pd=mod(-90-pd_obs,360);
          %%u_obs=cosd(pd).*hs_obs; v_obs=sind(pd).*hs_obs;
          %m_vec(50,lon_obs,lat_obs,u_obsc(i),v_obsc(i),'m','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5)
          %%m_vec(10,lon_obs,lat_obs,u_obsl(i),v_obsl(i),'c','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5)


          if ke>1%==length(expts)
            %[range,ln,lt]=m_lldist([nanmax(lon_hr) nanmin(lon_hr)],[nanmin(lat_hr) nanmin(lat_hr)],140);
            %m_line(wrapTo360(ln),lt,'color','k','linewi',2)
            %[range,ln,lt]=m_lldist([nanmin(lon_hr) nanmax(lon_hr)],[nanmax(lat_hr) nanmax(lat_hr)],140);
            %m_line(wrapTo360(ln),lt,'color','k','linewi',2)
            %[range,ln,lt]=m_lldist([nanmin(lon_hr) nanmin(lon_hr)],[nanmin(lat_hr) nanmax(lat_hr)],140);
            %m_line(wrapTo360(ln),lt,'color','k','linewi',2)
            %[range,ln,lt]=m_lldist([nanmax(lon_hr) nanmax(lon_hr)],[nanmin(lat_hr) nanmax(lat_hr)],140);
            %m_line(wrapTo360(ln),lt,'color','k','linewi',2)

            m_plot([lon_modm(:,1)],[lat_modm(:,1)],'k','linewidth',2);
            m_plot([lon_modm(:,end)],[lat_modm(:,end)],'k','linewidth',2);
            m_plot([lon_modm(1,:)],[lat_modm(1,:)],'k','linewidth',2);
            m_plot([lon_modm(end,:)],[lat_modm(end,:)],'k','linewidth',2);
          end

          display(['Finished ploting map']); toc

          %if dcol(i)==6
          %  mod_rmse=nanrmse(obs(:,dcol(i)),modeli');
          %  mod_bias=nanmean(modeli')-nanmean(obs(:,dcol(i)));
          %  legh=[legh,{[expt,' rmse=',num2str(mod_rmse,'%.2f'),', bias=',num2str(mod_bias,'%.2f')]}];
          %end

            %end % for expt=expts

          %end % for i=1:length(time)

        end % for expt=expts

        path_dm=strcat(path_fig,'paper_1/');
        system(['mkdir -p ',path_dm]);
        figname=[path_dm,file,'_',pcname,'_',vecname,'_',datestr(time(i),'YYYY_mm_DD_HH_MM'),'.png'];
        display(['Plotting: ',figname]);
        if save_fig==1
          display(['Saving: ',figname]);
          export_fig(gcf,figname,'-png','-r150');
          display(['End of script']);
          return       
          %close
          %clf('reset')
        else
          %pause
          %clf('reset')
        end

        if save_video==1
          if t==time_lima(1)
            videoName = [path_dm,'/',pcname,'_',datestr(time_lima(1),'yyyymmdd_HHMMSS'),'_',datestr(time_lima(end),'yyyymmdd_HHMMSS'),''];
            display(['Saving video name: ',videoName])
            vidObj = VideoWriter(videoName,'Uncompressed AVI');
            set(vidObj,'FrameRate',5)
            open(vidObj);
          end
          display(['Printing frame time: ',figname]);
          pause(0.5)
          M = getframe(gcf);
          writeVideo(vidObj,M);
        end

        clf('reset')
        set(gcf,'color',[1 1 1])

      end % for t=time_lima

    end % fobs==stations(end)

  end % if plot_downs

  % PLOTTING COASTAL MEAN FIELDS
  if plot_coastm

    scrsz=[1 1 1366 768];
    scrsz=get(0,'screensize');
    
    if fobs==stations(1); figure('position',scrsz,'color',[1 1 1],'visible','on'); end
    %hold on
    %set(gca,'fontsize',12,'fontweight','bold')
    
    ke=0;
    for expt=expts
      ke=ke+1;
    
      expt=expt{1};
      [ig,ltime]=grab_nc_sufix(expt);
    
      path_expt=[path_source,expt,'/']; % GLOBALWAVE/'];%2018/01/05/00'];
      
      path_dm=[path_expt,'matlab/'];
      system(['mkdir -p ',path_dm]);
     
      % creating model time_lima pair
      time_lima=time_limab;
      time_modp=time_lima(1):1/24:time_lima(end)+1; time_modp=time_modp';
      hs_mod=[]; % nan(1,length(time_lima)); 
	    tp_mod=hs_mod; pd_mod=hs_mod; time_mod=hs_mod;

      % creating annual time vector 
      time_limav=datevec(time_lima);
      s_year=time_limav(1,1); f_year=time_limav(end,1);
      time_years=s_year:f_year;

      consider_annual_files=1;
      if consider_annual_files==1

        % loop over years to look for annual mat files
        for y=time_years
          filename=[path_dm,file,'_',num2str(y),'.mat'];
          if ck_mod_mat==1 && exist(filename)==2
            display(['Loading: ',filename])
            load(filename) %,'time_mod','model') 

            time_mod=[time_mod;time_m];
            hs_mod  =[hs_mod;hs_m];
            tp_mod  =[tp_mod;tp_m];
            pd_mod  =[pd_mod;pd_m];

						% cutting time_lima to account for the previous loaded annual mat files
						[dif iloc]=nanmin(abs(time_mod(end)-time_lima));
			      time_lima=time_lima(iloc+1:end);
          end

        end

        % saving a netcdf file containing hs, tp, and pd, its lon and lat for each point with the whole timeseries if s_year==1983 and f_year==2022
				if s_year==1984 && f_year==2023

	      	filename=[path_dm,'nzwave_',file,num2str(s_year),'_',num2str(f_year),'.nc'];
	      	display(['Saving: ',filename])
          if exist(filename)==2
           system(['rm -rf ',filename]);
          end

          lon_obs=lon_mod(points(fobs,1)); lat_obs=lat_mod(points(fobs,2));

	      	nccreate(filename,'time','Dimensions',{'time',length(time_mod)},'Datatype','double','Format','classic');
	      	ncwrite(filename,'time',time_mod);
	      	nccreate(filename,'lon','Dimensions',{'lon',length(lon_obs)},'Datatype','double','Format','classic');
	      	ncwrite(filename,'lon',lon_obs);
	      	nccreate(filename,'lat','Dimensions',{'lat',length(lat_obs)},'Datatype','double','Format','classic');
	      	ncwrite(filename,'lat',lat_obs);
	      	%nccreate(filename,'hs','Dimensions',{'lon',length(lon_obs),'lat',length(lat_obs),'time',length(time_mod)},'Datatype','double','Format','classic');
	      	nccreate(filename,'hs','Dimensions',{'time',length(time_mod)},'Datatype','double','Format','classic');
	      	ncwrite(filename,'hs',hs_mod);
	      	nccreate(filename,'tp','Dimensions',{'time',length(time_mod)},'Datatype','double','Format','classic');
	      	ncwrite(filename,'tp',tp_mod);
	      	nccreate(filename,'pd','Dimensions',{'time',length(time_mod)},'Datatype','double','Format','classic');
	      	ncwrite(filename,'pd',pd_mod);

				end

      else
      	hs_mod=[]; % nan(1,length(time_lima)); 
	      tp_mod=hs_mod; pd_mod=hs_mod; time_mod=hs_mod;
			end


      checkin=0;  
      % Loop in model time
      for t=time_lima

        filename=[path_dm,file,'_',datestr(t,'YYYYmmDDHH'),'.mat'];
    
        if ck_mod_mat==1 && exist(filename)==2
      
          display(['Loading: ',filename])
          load(filename)%,'time_mod','model')
          if exist('ucur','var')~=1
            wlv=nan(size(hs)); ucur=nan(size(hs)); vcur=nan(size(hs));
          end  

        else

          %display(['NOT FOUND: ',filename])
          %display(['Run with plot_series=1 first'])
          %return

          ptime=datestr(t,'YYYY/mm/DD/HH/');
          ftime=datestr(t,'YYYYmmDDHH');
        
          pname=[path_expt,ptime];
        
          fname=[pname,ww3pre,'_',ftime,'-utc_',gnames{ig},'.nc'];
          %fname=[pname,'ww3g_',ftime,'-utc_',gnames{1},'.nc'];
          display(['Loading: ',fname]);
    
          if checkin==0; % t==time_lima(1)
            checkin=1; 
            % coastal points
            if plot_coastm==1
              lon_obs=lon_mod(points(fobs,1));
							lat_obs=lat_mod(points(fobs,2));
            end
            [dif ilon]=nanmin(abs(lon_mod-lon_obs));
            [dif ilat]=nanmin(abs(lat_mod-lat_obs));
    
            dis=sqrt((lon_mod(ilon)-lon_obs).^2 + (lat_mod(ilat)-lat_obs).^2)*111;
            display(['Distance between obs ',file,' and grid point is: ',num2str(dis),' km'])
          end
          
          time=squeeze(double(ncread(fname,'time')))./24+t; time=time(1:ltime-1);
          hs=squeeze(double(ncread(fname,'hsig',          [ilon ilat 1],[1 1 ltime]))); hs=hs(1:end-1);
          tp=squeeze(double(ncread(fname,'tpeak',         [ilon ilat 1],[1 1 ltime]))); tp=tp(1:end-1);
          pd=squeeze(double(ncread(fname,'peak_direction',[ilon ilat 1],[1 1 ltime]))); pd=pd(1:end-1);

          %wave energy (we) = (1/16)*rho*g*(hs.^2); % https://www.coastalwiki.org/wiki/Shallow-water_wave_theory#Introduction
          we=(1/16)*1025*9.8*(hs.^2);    

          display(['Saving: ',filename])
          save(filename,'time','tp','hs','pd') % ,'ds','tm01','tm02','hs','tm','wlv','ucur','vcur','we')

        end % if ck_mod_mat==1 && exist(filename)==2

        pd=pd(:);
 
        time_mod=[time_mod;time];
        hs_mod=[hs_mod;hs];
        tp_mod=[tp_mod;tp];
				pd_mod=[pd_mod;pd];
        clear model

        % check if it is the 31st of Dec to save the annual mat file
        tvec=datevec(t);
				if tvec(2)==12 & tvec(3)==31 

          if length(time_mod)>=365*length(time)
						[dif imodp]=nanmin(abs(datenum(tvec(1),1,1)-time_mod));
            if dif<0.05
						  [dif fmodp]=nanmin(abs(datenum(tvec(1),12,31,23,0,0)-time_mod));
							time_m=time_mod(imodp:fmodp);
							hs_m=hs_mod(imodp:fmodp);
							tp_m=tp_mod(imodp:fmodp);
							pd_m=pd_mod(imodp:fmodp);

						  filename=[path_dm,file,'_',num2str(tvec(1)),'.mat'];
						  display(['Saving: ',filename])
						  save(filename,'time_m','hs_m','tp_m','pd_m') % ,'ds','tm01','tm02','hs','tm','wlv','ucur','vcur','we')
            end
          end
				end

      end % for t=time_lima
    

      % Statistical analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ka=0;
      ccol=[3];
      for i=ccol
        ka=ka+1;

        % selected statistical variable        
        if i==1 % max wave height
					data=nanmax(hs_mod(:));
          tname='Max wave height (m)';
				elseif i==2 % mean wave height	
					data=nanmean(hs_mod(:));
          tname='Mean wave height (m)';
				elseif i==3 % 99th percentile wave height
          data=prctile(hs_mod,99);
          tname='99th percentile wave height (m)';
				elseif i==4 % X-year return period wave height (m)
          return_p=1;
          [data,GEV]=GPDfitwaves_funcion(hs_mod,time_mod,1,return_p); % =GPDfitwaves_funcion(hs_obs,time_obs,kind,rp)
          tname=[num2str(return_p),'-year return period wave height (m)'];
        end

        if fobs==stations(1)
          ax=subplot(1,length(ccol),ka);
          %subplot(length(mcol),1,i);
          set(gca,'fontsize',12,'fontweight','bold')
          hold on
          m_proj('lambert','long', [min(lon_hr) max(lon_hr)],'lat',[min(lat_hr) max(lat_hr)]);
          %if select_cpoint
          %  m_proj('lambert','long', [min(lon_mod(points(:,1)))-.5 max(lon_mod(points(:,1)))+.5],'lat',[min(lat_mod(points(:,2)))-.5 max(lat_mod(points(:,2)))+.5]);
          %else
          %  m_proj('lambert','long', [lon_c-.5 lon_c+.5],'lat',[lat_c-.5 lat_c+.5]);
          %end
          if ke==1 & plot_obs==1
            %plot(time_obs,obs(:,dcol(i))','.','color',[0 .7 0],'linewidth',2)
          end
        end
        %m_scatter(lon_mod(points(fobs,1)),lat_mod(points(fobs,2)),data,10,'filled')

        cmin=1; cmax=9;
        cmap=avhrr; 
        dmap=linspace(cmin,cmax,length(cmap))';
        [dif ij]=nanmin(abs(dmap-data));
        h=m_plot(lon_mod(points(fobs,1)),lat_mod(points(fobs,2)),'.','markersize',20.0);
        set(h,'color',[cmap(ij,:)])

        if fobs==stations(1)
          colormap(ax,cmap)
          caxis([cmin cmax])
          cb=colorbar;%('southoutside');

          m_gshhs_i('patch',[.7 .7 .7],'edgecolor','k','linewidth',.5);
          %m_grid('xtick',[round(loni):2:round(lonf)],'xticklabels',[round(loni):2:round(lonf)],'fontname','helvetica','fontsize',20,'fontweight','bold');
          m_grid('box','fancy','fontname','helvetica','fontsize',14,'fontweight','bold');

          title([expt,' ',tname,' between ',datestr(time_limab(1),'DD mmm YYYY'),' and ',datestr(time_limab(end),'DD mmm YYYY')])%,', mean=',num2str(nanmean(data(:)),'%.2f')]);
          
        end
 
      end % for i=1:length(mcol)
    
    end % for expt=expts

    if fobs==length(points)
      path_dm=strcat(path_fig,expt,'/');
      system(['mkdir -p ',path_dm]);
      %figname=[path_dm,expt,'_wave_coastm_',datestr(time_limab(1),'YYYY_mm_DD_HH_MM'),'_',datestr(time_limab(end),'YYYY_mm_DD_HH_MM')];
      figname=[path_dm,expt,'_wave_coastm_',replace(replace(tname,' ','_'),'(m)',''),'years_',datestr(time_limab(1),'YYYY'),'_',datestr(time_limab(end),'YYYY')];
      display(['Plotting: ',figname]);
      if save_fig==1
        display(['Saving: ',figname]);
        export_fig(gcf,[figname,'.png'],'-png','-r150');
        saveas(gcf,[figname,'.fig'],'fig')
        %close
        %clf('reset')
      end
    end

  end % if plot_coastm


  % SPATIAL MEAN FUNCTION
  if plot_mserie || plot_mmap
  
    scrsz=[1 1 1366 768];
    scrsz=[1 1 1920 1080];
    scrsz=[9 1 1910 990];
    %scrsz=get(0,'screensize');
      
    legh={'Obs'}; legt={'Obs'};
    
    ke=0; ka=0;
    for expt=expts
      ke=ke+1;
    
      expt=expt{1};
      [ig,ltime]=grab_nc_sufix(expt);
      path_expt=[path_source,expt,'/']; 
      path_dm=[path_expt,'matlab/'];
      system(['mkdir -p ',path_dm]);

      hs_m_mean=[];hs_m_maxi=[]; hs_m_meanf=[];hs_m_maxif=[]; 
    
      checkin=0;  
      % Loop in model time
      kt=0;
      for t=time_lima
        kt=kt+1;
        filename=[path_dm,'daily_mean_',datestr(t,'YYYYmmDDHH'),'.mat'];
    
        if ck_mod_mat==1 && exist(filename)==2
      
          display(['Loading: ',filename])
          load(filename)%,'time_mod','model')
      
        else

          ptime=datestr(t,'YYYY/mm/DD/HH/');
          ftime=datestr(t,'YYYYmmDDHH');
          pname=[path_expt,ptime];
          fname=[pname,ww3pre,'_',ftime,'-utc_',gnames{ig},'.nc'];
          display(['Loading: ',fname]);
    
          if checkin==0; % t==time_lima(1)
            checkin=1; 
            lon_mod=double(ncread(fname,'lon'));
            lat_mod=double(ncread(fname,'lat'));
          end
          
          time=squeeze(double(ncread(fname,'time')))./24+t; time=time(1:ltime-2);
          hs=squeeze(double(ncread(fname,'hsig',          [1 1 1],[Inf Inf ltime-1]))); 
          tp=squeeze(double(ncread(fname,'tpeak',         [1 1 1],[Inf Inf ltime-1])));
          %tm01=squeeze(double(ncread(fname,'tmean01',     [1 1 1],[Inf Inf ltime-1])));
          %tm=squeeze(double(ncread(fname,'tmean',         [1 1 1],[Inf Inf ltime-1])));
          %pd=squeeze(double(ncread(fname,'peak_direction',[1 1 1],[Inf Inf ltime-1])));
          %md=squeeze(double(ncread(fname,'mean_direction',[1 1 1],[Inf Inf ltime-1])));
           
          hs_mod_mean=nanmean(hs,3);
          tp_mod_mean=nanmean(tp,3);
          
          display(['Saving: ',filename])
          save(filename,'hs_mod_mean','tp_mod_mean')
 
        end % if ck_mod_mat==1 && exist(filename)==2

        % sum of the wave fields
        if kt==1
					hs_sum=hs_mod_mean; tp_sum=tp_mod_mean;
				else
          hs_sum=hs_sum+hs_mod_mean; tp_sum=tp_sum+tp_mod_mean;
				end

      end % for t=time_lima

      % computing stats of mean fields
			hs_m_mean=hs_sum/kt; tp_m_mean=tp_sum/kt;
      we_m_mean=(1/16)*1025*9.8*(hs_m_mean.^2);


      if plot_mserie

         figure('position',scrsz,'color',[1 1 1],'visible','on');
         hold on
         set(gca,'fontsize',12,'fontweight','bold')

        for i=1:length(mcol)

          display(['data=',mvals{mcol(i)},';']);
          eval(['data=',mvals{mcol(i)},';']);

          subplot(length(mcol),1,i);
          set(gca,'fontsize',12,'fontweight','bold')
          hold on
          if ke==1 & plot_obs==1
            %plot(time_obs,obs(:,dcol(i))','.','color',[0 .7 0],'linewidth',2)
          end
          plot(time_lima,data,'-','color',colors{ke},'linewidth',2)
          
          if set_xlim==1;
            xlim([time_lima(1)-2 time_lima(end)+2])
          end
        
          %if i==3
            datetick('x','dd/mmm/yyyy','keeplimits')
          %end
          title([replace(atype,'_',' '),' ',replace(mvals{mcol(i)},'_',' ')])
          grid('on')
          
          %modeli=interp1(time_mod,model.data(:,dcol(i)),time_obs);
          
          %if dcol(i)==6
          %  mod_rmse=nanrmse(obs(:,dcol(i)),modeli');
          %  mod_corr=nancorr(obs(:,dcol(i)),modeli');
          %  mod_bias=nanmean(modeli')-nanmean(obs(:,dcol(i)));
          %  legh=[legh,{[expt,' rmse=',num2str(mod_rmse,'%.2f'),' corr=',num2str(mod_corr,'%.2f'),', bias=',num2str(mod_bias,'%.2f')]}];
          %end
          legh=[expt,' mean=',num2str(nanmean(data),'%.2f')];
          legend([legh],'location','best')
           
        end

      end % if plot_mserie

      if plot_mmap

        if ka==0
          scrsz=[2    42   958   953];
          figure('position',scrsz,'color',[1 1 1],'visible','on');
          hold on
          set(gca,'fontsize',12,'fontweight','bold')
        end

        for i=1:length(mcol)
          ka=ka+1;
          
          display(['data=',mvalm{mcol(i)},';']);
          eval(['data=',mvalm{mcol(i)},';']);

          if length(runs)==1 
            lsub=1; csub=1;
            if length(mcol)==1
              lsub=1; csub=1;
            elseif length(mcol)==2
              lsub=2; csub=1;
            elseif mod(length(mcol),2)==1
              lsub=round(length(mcol)/2); csub=round(length(mcol)/2);
            elseif mod(length(mcol),2)==0 
              lsub=round(length(mcol)/2); csub=round(length(mcol)/2)-1;
            end
          else
            if mod(length(runs),2)==1 
              lsub=round(length(runs)/2); csub=round(length(runs)/2)-1;
            elseif mod(length(runs),2)==0 
              lsub=round(length(runs)/2); csub=round(length(runs)/2);
            end
          end

          ax=subplot(lsub,csub,ka);
          %subplot(length(mcol),1,i);
          set(gca,'fontsize',12,'fontweight','bold')
          hold on
          m_proj('lambert','long', [min(lon_hr) max(lon_hr)],'lat',[min(lat_hr) max(lat_hr)]);
          if ke==1 & plot_obs==1
            %plot(time_obs,obs(:,dcol(i))','.','color',[0 .7 0],'linewidth',2)
          end
          m_pcolor(lon_mod,lat_mod,data')
          shading flat
          avalname=mvalm{mcol(i)};
          if strcmp(avalname(end-3:end),'mean')
            colormap(ax,avhrr)
            %caxis([0 1.5])
          elseif strcmp(avalname(end-3:end),'bias')
            colormap(ax,cmocean('balance'))
            caxis([-1 1])
          elseif strcmp(avalname(end-3:end),'corr')
            colormap(ax,cmocean('matter'))
            [cs,h]=contour(lon_mod,lat_mod,data','w');%[.5:.2:1],'w');
            clabel(cs,h,'color',[1 1 1],'fontsize',8,'fontweight','bold','LabelSpacing',2000)
            caxis([0 1])
          end

          cb=colorbar;%('southoutside');
          %set(get(cb,'ylabel'),'string','','fontsize',12,'fontweight','bold');
          %set(cb,'fontsize',12,'fontweight','bold');

          %[cs,h]=m_contour(lon_mod,lat_mod,-depth_mod',[-1000 -1000],'color',[1 1 1],'linewidth',2);
          %clabel(cs,h,'color',[0 0 0],'fontsize',10,'fontweight','bold','LabelSpacing',2000)

          %[cs,h]=m_contour(lon_aviso,lat_aviso,z_aviso_avg,'color','k');
          %[cs,h]=m_contour(lon_aviso,lat_aviso,z_aviso_avg,[.76:.04:.88],'color','k');
          %clabel(cs,h,'color',[0 0 0],'fontsize',14,'fontweight','bold','LabelSpacing',2000) % each contour line represents 0.02 increase/decrease

          m_gshhs_i('patch',[.7 .7 .7],'edgecolor','k','linewidth',.5);
          %%m_grid('xtick',[round(loni):2:round(lonf)],'xticklabels',[round(loni):2:round(lonf)],'fontname','helvetica','fontsize',20,'fontweight','bold');
          m_grid('fontname','helvetica','fontsize',14,'fontweight','bold');

          tmvalm=replace(mvalm{mcol(i)},'_',' '); tmvalm=replace(tmvalm,'m ',''); 
          title([expt,' ',tmvalm,', mean=',num2str(nanmean(data(:)),'%.2f')]);
          
          %modeli=interp1(time_mod,model.data(:,dcol(i)),time_obs);
          
          %if dcol(i)==6
          %  mod_rmse=nanrmse(obs(:,dcol(i)),modeli');
          %  mod_corr=nancorr(obs(:,dcol(i)),modeli');
          %  mod_bias=nanmean(modeli')-nanmean(obs(:,dcol(i)));
          %  legh=[legh,{[expt,' rmse=',num2str(mod_rmse,'%.2f'),' corr=',num2str(mod_corr,'%.2f'),', bias=',num2str(mod_bias,'%.2f')]}];
          %end
           
        end

      end % if plot_mmap
    
    end % for expt=expts
    
    if ~isempty(find(dcol==6))
      %subplot(length(dcol),1,find(dcol==6));
      %legend([legh],'location','best')
    end
    
    if save_fig==1
      %export_fig(gcf,'steephead_period','-png','-r150');
    end

  end % if plot_mserie || plot_mmap


  % ALTIMETER FUNCTIONS

  if plot_aserie || plot_amap
  
    scrsz=[1 1 1366 768];
    scrsz=[1 1 1920 1080];
    scrsz=[9 1 1910 990];

    %scrsz=get(0,'screensize');
    
      
    legh={'Obs'}; legt={'Obs'};
    
    ke=0; ka=0;
    for expt=expts
      ke=ke+1;
    
      expt=expt{1};
      [ig,ltime]=grab_nc_sufix(expt);
    
      path_expt=[path_source,expt,'/']; 
      
      path_dm=[path_expt,'matlab/'];
      system(['mkdir -p ',path_dm]);

      hs_m_mean=[];hs_m_maxi=[]; hs_m_meanf=[];hs_m_maxif=[]; 

      tvec=datevec(time_lima(1));
      filenamey=[path_dm,atype,'_',num2str(tvec(1,1)),'.mat'];
      if ck_mod_mat==1 && exist(filenamey)==2

				display(['Loading: ',filenamey])
				load(filenamey)%,'time_mod','model')

      else
    
        checkin=0;  
        % Loop in model time
        kt=0;
        for t=time_lima
          kt=kt+1;
          filename=[path_dm,atype,'_',datestr(t,'YYYYmmDDHH'),'.mat'];
    
          if ck_mod_mat==1 && exist(filename)==2
        
            display(['Loading: ',filename])
            load(filename)%,'time_mod','model')
        
          else

            ptime=datestr(t,'YYYY/mm/DD/HH/');
            ftime=datestr(t,'YYYYmmDDHH');
          
            pname=[path_expt,ptime];
          
            fname=[pname,ww3pre,'_',ftime,'-utc_',gnames{ig},'.nc'];
            display(['Loading: ',fname]);
    
            if checkin==0; % t==time_lima(1)
              checkin=1; 
              lon_mod=double(ncread(fname,'lon'));
              lat_mod=double(ncread(fname,'lat'));
            end
            
            time=squeeze(double(ncread(fname,'time')))./24+t; time=time(1:ltime-2);
            hs=squeeze(double(ncread(fname,'hsig',          [1 1 1],[Inf Inf ltime-1]))); 
            
            hs_mod_mean=nanmean(hs,3);
            hs_mod_maxi=nanmax(hs,3);
            
            % filter mod
            filt_mod_spa=1;
            if filt_mod_spa
              gauss = fspecial('average',[54 54]); % NZWAVE
              hs_mod_meanf=imfilter(hs_mod_mean,gauss,'replicate');
              hs_mod_maxif=imfilter(hs_mod_mean,gauss,'replicate');
            end
             
            %interpolating 
            hs_mod_mean=griddata(lon_mod,lat_mod',hs_mod_mean',lon_cop,lat_cop')';
            hs_mod_maxi=griddata(lon_mod,lat_mod',hs_mod_maxi',lon_cop,lat_cop')';
            hs_mod_meanf=griddata(lon_mod,lat_mod',hs_mod_meanf',lon_cop,lat_cop')';
            hs_mod_maxif=griddata(lon_mod,lat_mod',hs_mod_maxif',lon_cop,lat_cop')';
            display(['Saving: ',filename])
            save(filename,'hs_mod_mean','hs_mod_maxi','hs_mod_meanf','hs_mod_maxif')
  
          end % if ck_mod_mat==1 && exist(filename)==2

          if crop_l4
            %if strncmp(expt,'NZWAVE-HR',9)~=1
              %[dif ilon]=nanmin(abs(lon_cop-min(lon_hr))); [dif flon]=nanmin(abs(lon_cop-max(lon_hr)));
              %[dif ilat]=nanmin(abs(lat_cop-min(lat_hr))); [dif flat]=nanmin(abs(lat_cop-max(lat_hr)));
              hs_mod_mean =hs_mod_mean(ilonc:flonc,ilatc:flatc);
              hs_mod_meanf=hs_mod_meanf(ilonc:flonc,ilatc:flatc);
		          hs_mod_maxi =hs_mod_maxi(ilonc:flonc,ilatc:flatc);
		          hs_mod_maxif=hs_mod_maxif(ilonc:flonc,ilatc:flatc);
            %end
          end

          % ALTIMETER FUNCTIONS

          % concatenating 
          hs_m_mean=cat(3,hs_m_mean,hs_mod_mean); 
          hs_m_maxi=cat(3,hs_m_maxi,hs_mod_maxi); 
          hs_m_meanf=cat(3,hs_m_meanf,hs_mod_meanf); 
          hs_m_maxif=cat(3,hs_m_maxif,hs_mod_maxif); 
     
          icop=find(time_cop==t); 
 
          % time series stats
          hs_mean_rmse(kt)=nanrmse(hs_cop_mean(:,:,icop),hs_mod_mean);
          hs_maxi_rmse(kt)=nanrmse(hs_cop_maxi(:,:,icop),hs_mod_maxi);
          hs_diff=(hs_mod_mean-hs_cop_mean(:,:,icop)); 
          hs_mean_bias(kt)=nanmean(hs_diff(:));
          hs_diff=(hs_mod_maxi-hs_cop_maxi(:,:,icop));
          hs_maxi_bias(kt)=nanmean(hs_diff(:));
          hs_mean_corr(kt)=nancorr(hs_cop_mean(:,:,icop),hs_mod_mean);
          hs_maxi_corr(kt)=nancorr(hs_cop_maxi(:,:,icop),hs_mod_maxi);

          hs_mean_frmse(kt)=nanrmse(hs_cop_mean(:,:,icop),hs_mod_meanf);
          hs_maxi_frmse(kt)=nanrmse(hs_cop_maxi(:,:,icop),hs_mod_maxif);
          hs_diff=(hs_mod_meanf-hs_cop_mean(:,:,icop)); 
          hs_mean_fbias(kt)=nanmean(hs_diff(:));
          hs_diff=(hs_mod_maxif-hs_cop_maxi(:,:,icop));
          hs_maxi_fbias(kt)=nanmean(hs_diff(:));
          hs_mean_fcorr(kt)=nancorr(hs_cop_mean(:,:,icop),hs_mod_meanf);
          hs_maxi_fcorr(kt)=nancorr(hs_cop_maxi(:,:,icop),hs_mod_maxif);

        end % for t=time_lima

        % save annual file if time starts on 1st of Jan and finishes on the 31st of Dec
        tvec=datevec(time_lima);
			  if tvec(1,2)==1 & tvec(1,3)==1 & tvec(end,2)==12 & tvec(end,3)==31 
          %filenamey=[path_dm,atype,'_',datestr(tvec(1,1),'YYYY'),'.mat'];
			  	display(['Saving: ',filenamey])
			  	save(filenamey,'hs_m_mean','hs_m_maxi','hs_m_meanf','hs_m_maxif',... 
		       'hs_mean_rmse','hs_maxi_rmse','hs_mean_bias','hs_maxi_bias','hs_mean_corr','hs_maxi_corr',... 
		       'hs_mean_frmse','hs_maxi_frmse','hs_mean_fbias','hs_maxi_fbias','hs_mean_fcorr','hs_maxi_fcorr')
        end

      end % if ck_mod_mat==1 && exist(filenamey)==2


      if plot_amap
        icopi=find(time_cop==time_lima(1)); 
        icopf=find(time_cop==time_lima(end)); 
        for i=1:size(hs_cop_mean,1)
          for ii=1:size(hs_cop_mean,2)
            hs_mean_mrmse(i,ii)=nanrmse(hs_m_mean(i,ii,:),hs_cop_mean(i,ii,icopi:icopf));
            hs_maxi_mrmse(i,ii)=nanrmse(hs_m_maxi(i,ii,:),hs_cop_maxi(i,ii,icopi:icopf));
            hs_mean_mbias(i,ii)=nanmean(hs_m_mean(i,ii,:)-hs_cop_mean(i,ii,icopi:icopf));
            hs_maxi_mbias(i,ii)=nanmean(hs_m_maxi(i,ii,:)-hs_cop_maxi(i,ii,icopi:icopf));
            hs_mean_mcorr(i,ii)=nancorr(hs_m_mean(i,ii,:),hs_cop_mean(i,ii,icopi:icopf));
            hs_maxi_mcorr(i,ii)=nancorr(hs_m_maxi(i,ii,:),hs_cop_maxi(i,ii,icopi:icopf));

            hs_mean_fmrmse(i,ii)=nanrmse(hs_m_meanf(i,ii,:),hs_cop_mean(i,ii,icopi:icopf));
            hs_maxi_fmrmse(i,ii)=nanrmse(hs_m_maxif(i,ii,:),hs_cop_maxi(i,ii,icopi:icopf));
            hs_mean_fmbias(i,ii)=nanmean(hs_m_meanf(i,ii,:)-hs_cop_mean(i,ii,icopi:icopf));
            hs_maxi_fmbias(i,ii)=nanmean(hs_m_maxif(i,ii,:)-hs_cop_maxi(i,ii,icopi:icopf));
            hs_mean_fmcorr(i,ii)=nancorr(hs_m_meanf(i,ii,:),hs_cop_mean(i,ii,icopi:icopf));
            hs_maxi_fmcorr(i,ii)=nancorr(hs_m_maxif(i,ii,:),hs_cop_maxi(i,ii,icopi:icopf));

            % compute obs stats 
            hs_mean_ostd(i,ii)=nanstd(hs_cop_mean(i,ii,icopi:icopf));
						hs_maxi_ostd(i,ii)=nanstd(hs_cop_maxi(i,ii,icopi:icopf));
	          hs_mean_omean(i,ii)=nanmean(hs_cop_mean(i,ii,icopi:icopf));
	          hs_maxi_omean(i,ii)=nanmean(hs_cop_maxi(i,ii,icopi:icopf));
	          hs_mean_omin(i,ii)=nanmin(hs_cop_mean(i,ii,icopi:icopf));
	          hs_maxi_omin(i,ii)=nanmin(hs_cop_maxi(i,ii,icopi:icopf));
	          hs_mean_omax(i,ii)=nanmax(hs_cop_mean(i,ii,icopi:icopf));
	          hs_maxi_omax(i,ii)=nanmax(hs_cop_maxi(i,ii,icopi:icopf));

          end
        end
      end
    
      if plot_aserie

         figure('position',scrsz,'color',[1 1 1],'visible','on');
         hold on
         set(gca,'fontsize',12,'fontweight','bold')

        for i=1:length(acol)

          display(['data=',avals{acol(i)},';']);
          eval(['data=',avals{acol(i)},';']);

          subplot(length(acol),1,i);
          set(gca,'fontsize',12,'fontweight','bold')
          hold on
          if ke==1 & plot_obs==1
            %plot(time_obs,obs(:,dcol(i))','.','color',[0 .7 0],'linewidth',2)
          end
          plot(time_lima,data,'-','color',colors{ke},'linewidth',2)
          
          if set_xlim==1;
            xlim([time_lima(1)-2 time_lima(end)+2])
          end
        
          %if i==3
            datetick('x','dd/mmm/yyyy','keeplimits')
          %end
          title([replace(atype,'_',' '),' ',replace(avals{acol(i)},'_',' ')])
          grid('on')
          
          %modeli=interp1(time_mod,model.data(:,dcol(i)),time_obs);
          
          %if dcol(i)==6
          %  mod_rmse=nanrmse(obs(:,dcol(i)),modeli');
          %  mod_corr=nancorr(obs(:,dcol(i)),modeli');
          %  mod_bias=nanmean(modeli')-nanmean(obs(:,dcol(i)));
          %  legh=[legh,{[expt,' rmse=',num2str(mod_rmse,'%.2f'),' corr=',num2str(mod_corr,'%.2f'),', bias=',num2str(mod_bias,'%.2f')]}];
          %end
          legh=[expt,' mean=',num2str(nanmean(data),'%.2f')];
          legend([legh],'location','best')
           
        end

      end % if plot_aserie

      if plot_amap

        if ka==0
          scrsz=[1 41 1920 962];
          %scrsz=[1 1 1366 768];
          scrsz=[2    42   958   953];
          scrsz=get(0,'screensize');

          figure('position',scrsz,'color',[1 1 1],'visible','on');
          hold on
          set(gca,'fontsize',12,'fontweight','bold')
        end

        for i=1:length(acol)
          ka=ka+1;
          
          display(['data=',avalm{acol(i)},';']);
          eval(['data=',avalm{acol(i)},';']);

          if length(runs)==1 
            lsub=1; csub=1;
            if length(acol)==1
              lsub=1; csub=1;
            elseif mod(length(acol),2)==1
              lsub=round(length(acol)/2); csub=round(length(acol)/2);
            elseif mod(length(acol),2)==0 
              lsub=round(length(acol)/2); csub=round(length(acol)/2)-1;
            end
          else
            if mod(length(runs),2)==1 
              lsub=round(length(runs)/2); csub=round(length(runs)/2);
            elseif mod(length(runs),2)==0 
              lsub=round(length(runs)/2); csub=round(length(runs)/2)+1;
            end
          end

          % plotting obs if it is the first experiment
          % plot map of selected variable
          if plot_obs==1 & ka==1

            display(['obsdata=',ovalm{acol(i)},';']); % ovalm = contains the name of the Obs VAriabLes for Maps 
            eval(['obsdata=',ovalm{acol(i)},';']);

            ax=subplot(lsub,csub,ka);
            %subplot(length(acol),1,i);
            set(gca,'fontsize',10,'fontweight','bold')
            hold on
            m_proj('lambert','long', [min(lon_cop)+1 max(lon_cop)-1],'lat',[min(lat_cop) max(lat_cop)-1]);
            loni=nanmin(lon_cop(:)); lonf=nanmax(lon_cop(:));

            m_pcolor(lon_cop-1,lat_cop-1,obsdata')
            shading flat

            avalname=ovalm{acol(i)};
            if strcmp(avalname(end-2:end),'std')
              colormap(ax,parula)%cmocean('haline'))
              %caxis([0 1.5])
            elseif strcmp(avalname(end-3:end),'mean')
              colormap(ax,cmocean('thermal'))
              %caxis([0 6])
            elseif strcmp(avalname(end-3:end),'corr')
              colormap(ax,cmocean('matter'))
              [cs,h]=contour(lon_cop,lat_cop,data','w');%[.5:.2:1],'w');
              clabel(cs,h,'color',[1 1 1],'fontsize',8,'fontweight','bold','LabelSpacing',2000)
              caxis([0 1])
            end

            cb=colorbar; %('southoutside');
            %set(get(cb,'ylabel'),'string','','fontsize',12,'fontweight','bold');
            %set(cb,'fontsize',12,'fontweight','bold');
            %h=colorbar;
            %set(h,'position',[0.915 0.24 0.011 0.44]) % [left bottom width height]
            %set(get(h,'ylabel'),'string','Probability density estimate','fontsize',10,'fontweight','bold');


            %[cs,h]=m_contour(lon_mod,lat_mod,-depth_mod',[-1000 -1000],'color',[1 1 1],'linewidth',2);
            %clabel(cs,h,'color',[0 0 0],'fontsize',10,'fontweight','bold','LabelSpacing',2000)

            %[cs,h]=m_contour(lon_aviso,lat_aviso,z_aviso_avg,'color','k');
            %[cs,h]=m_contour(lon_aviso,lat_aviso,z_aviso_avg,[.76:.04:.88],'color','k');
            %clabel(cs,h,'color',[0 0 0],'fontsize',14,'fontweight','bold','LabelSpacing',2000) % each contour line represents 0.02 increase/decrease

            m_gshhs_i('patch',[.7 .7 .7],'edgecolor','k','linewidth',.5);
            m_grid('xtick',[round(loni):14:round(lonf)],'xticklabels',[round(loni):14:round(lonf)],'fontname','helvetica','fontsize',12,'fontweight','bold');

            tavalm=replace(ovalm{acol(i)},'_',' '); tavalm=replace(tavalm,'hs',''); 
            title(['Obs. Hs ',tavalm,' = ',num2str(nanmean(obsdata(:)),'%.2f'),' (m)']);

            if strncmp(expt,'NZWAVE-HR',9)~=1
              if crop_l4==0
                datacrop=obsdata(ilonc:flonc,ilatc:flatc);
                m_text(165,-29,['mean=',num2str(nanmean(datacrop(:)),'%.2f')],'fontweight','bold','fontsize',8);
                [lon_modm,lat_modm]=meshgrid(lon_hr,lat_hr);
                m_plot([lon_modm(:,1)],[lat_modm(:,1)],'k','linewidth',2);
                m_plot([lon_modm(:,end)],[lat_modm(:,end)],'k','linewidth',2);
                m_plot([lon_modm(1,:)],[lat_modm(1,:)],'k','linewidth',2);
                m_plot([lon_modm(end,:)],[lat_modm(end,:)],'k','linewidth',2);
              end
            end

            pause(0.5)

            ka=ka+1;
          end % if plot_obs

          ax=subplot(lsub,csub,ka);
          %subplot(length(acol),1,i);
          set(gca,'fontsize',10,'fontweight','bold')
          hold on
          if strncmp(expt,'NZWAVE-HR',9)~=1
            m_proj('lambert','long', [min(lon_cop)+1 max(lon_cop)-1],'lat',[min(lat_cop) max(lat_cop)-1]);
            loni=nanmin(lon_cop(:)); lonf=nanmax(lon_cop(:));
          else
            m_proj('lambert','long', [min(lon_hr)+.5 max(lon_hr)],'lat',[min(lat_hr)+1 max(lat_hr)]);
            loni=nanmin(lon_hr(:)); lonf=nanmax(lon_hr(:));
          end
          if ke==1 & plot_obs==1
            %plot(time_obs,obs(:,dcol(i))','.','color',[0 .7 0],'linewidth',2)
          end

          if strncmp(expt,'NZWAVE-HR',9)~=1
            m_pcolor(lon_cop-1,lat_cop-1,data')
          else
            m_pcolor(lon_cop-1,lat_cop-1,data')
					end
          shading flat
          avalname=avalm{acol(i)};
          if strcmp(avalname(end-3:end),'rmse')
            colormap(ax,avhrr)
            caxis([0.2 0.8])
          elseif strcmp(avalname(end-3:end),'bias')
            colormap(ax,cmocean('balance'))
            caxis([-1 1])
          elseif strcmp(avalname(end-3:end),'corr')
            colormap(ax,cmocean('matter'))
            [cs,h]=contour(lon_cop,lat_cop,data','w');%[.5:.2:1],'w');
            clabel(cs,h,'color',[1 1 1],'fontsize',8,'fontweight','bold','LabelSpacing',2000)
            caxis([0 1])
          end

          cb=colorbar; %('southoutside');
          %set(get(cb,'ylabel'),'string','','fontsize',12,'fontweight','bold');
          %set(cb,'fontsize',12,'fontweight','bold');

          %[cs,h]=m_contour(lon_mod,lat_mod,-depth_mod',[-1000 -1000],'color',[1 1 1],'linewidth',2);
          %clabel(cs,h,'color',[0 0 0],'fontsize',10,'fontweight','bold','LabelSpacing',2000)

          %[cs,h]=m_contour(lon_aviso,lat_aviso,z_aviso_avg,'color','k');
          %[cs,h]=m_contour(lon_aviso,lat_aviso,z_aviso_avg,[.76:.04:.88],'color','k');
          %clabel(cs,h,'color',[0 0 0],'fontsize',14,'fontweight','bold','LabelSpacing',2000) % each contour line represents 0.02 increase/decrease

          m_gshhs_i('patch',[.7 .7 .7],'edgecolor','k','linewidth',.5);
          if strncmp(expt,'NZWAVE-HR',9)~=1
            m_grid('xtick',[round(loni):14:round(lonf)],'xticklabels',[round(loni):14:round(lonf)],'fontname','helvetica','fontsize',12,'fontweight','bold');
          else
            m_grid('xtick',[round(loni):8:round(lonf)],'xticklabels',[round(loni):8:round(lonf)],'fontname','helvetica','fontsize'  ,12,'fontweight','bold');
          end
          %m_grid('fontname','helvetica','fontsize',14,'fontweight','bold');

          tavalm=replace(avalm{acol(i)},'_',' '); tavalm=replace(tavalm,'hs',''); 
          title([expt,' ',tavalm,' = ',num2str(nanmean(data(:)),'%.2f')]);

          if strncmp(expt,'NZWAVE-HR',9)~=1
            if crop_l4==0
              datacrop=data(ilonc:flonc,ilatc:flatc);
              m_text(165,-29,['mean=',num2str(nanmean(datacrop(:)),'%.2f')],'fontweight','bold','fontsize',8);
              [lon_modm,lat_modm]=meshgrid(lon_hr,lat_hr);
              m_plot([lon_modm(:,1)],[lat_modm(:,1)],'k','linewidth',2);
              m_plot([lon_modm(:,end)],[lat_modm(:,end)],'k','linewidth',2);
              m_plot([lon_modm(1,:)],[lat_modm(1,:)],'k','linewidth',2);
              m_plot([lon_modm(end,:)],[lat_modm(end,:)],'k','linewidth',2);
            end
          end

          %modeli=interp1(time_mod,model.data(:,dcol(i)),time_obs);
          
          %if dcol(i)==6
          %  mod_rmse=nanrmse(obs(:,dcol(i)),modeli');
          %  mod_corr=nancorr(obs(:,dcol(i)),modeli');
          %  mod_bias=nanmean(modeli')-nanmean(obs(:,dcol(i)));
          %  legh=[legh,{[expt,' rmse=',num2str(mod_rmse,'%.2f'),' corr=',num2str(mod_corr,'%.2f'),', bias=',num2str(mod_bias,'%.2f')]}];
          %end
 
          pause(0.5)

        end

      end % if plot_amap
    
    end % for expt=expts
    
    if ~isempty(find(dcol==6))
      %subplot(length(dcol),1,find(dcol==6));
      %legend([legh],'location','best')
    end

    path_dm=strcat(path_fig,expt,'/');
    system(['mkdir -p ',path_dm]);
    if plot_amap==1
      figname=[path_dm,avalm{acol(1)},'_',datestr(time_lima(1),'YYYYmmDD'),'_',datestr(time_lima(end),'YYYYmmDD'),'.png'];
    elseif plot_aserie==1
      figname=[path_dm,avals{acol(1)},'_',datestr(time_lima(1),'YYYYmmDD'),'_',datestr(time_lima(end),'YYYYmmDD'),'.png'];
		end
    display(['Plotting: ',figname]);
    if save_fig==1
      display(['Saving: ',figname]);
      export_fig(gcf,figname,'-png','-r150');
      %close
      %clf('reset')
    end

    
  end % if plot_aserie and plot_amap
  
  
end % fobs=fstations


% compute maps of annual average wave height and plot them


toc
