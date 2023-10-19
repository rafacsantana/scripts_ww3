clear
close all

time=round(now)-4; % look for files generated 4 days ago
time=datenum(2021,6,30,0,0,0);

%Loading OBS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/';
path_obs=[path,'data/obs/'];
file='SteepHead_SP_FullRecord_QC';
file_obs=file;
display(['Loading: ',path,file,'.mat']);
load([path_obs,file,'.mat'])%,'time','obs')


tstation={'Sinclair Head','Ohau Head','Mana Island','Makara'  ,'Karori Rock','Wairewa Lake Forysyth'};
stations={'sinclair_head','ohau_head','mana_island','makara'  ,'karori_rock','wairewa_lake_forsyth'};
lat_obss=[-41.366885     ,-42.256063 ,-41.083327   ,-41.203393, -41.344623  ,-43.84110];
lon_obss=[174.716761     , 173.853341,174.763960   ,174.702391, 174.650430  ,172.71835];

tstation={tstation{end}};
stations={stations{end}};
lon_obss=lon_obss(end); lat_obss=lat_obss(end);
lon_obss=lon_obss(end); lat_obss=lat_obss(end);

lon_obs=lon_obss; lat_obs=lat_obss;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path='/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/';
%path='/scale_wlg_devoper/filesets/archive/ecoconnect/EFS/';
ptime=datestr(time,'/YYYY/mm/DD/HH/');
ftime=datestr(time,'YYYYmmDDHH');

grids={'GLOBALWAVE','NZWAVE','NZWAVE-HR','TONGAWAVE'};
gnames={'globalwave+globalum','nzwave+nzlam','nzwave_hr+nzcsm','tongawave+globalum'};

% Loading GRIDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scrsz=[1 1 1366 768];
scrsz=get(0,'screensize');
figure('position',scrsz,'color',[1 1 1],'visible','on')
hold on
set(gca,'fontsize',12,'fontweight','bold')

colors={'b','r','k','k'};

for i=[2]

  pname=[path,grids{i},ptime];
  if i==2;
    fname=[pname,'ww3g_',ftime,'-utc_',gnames{i},'.nc'];
    depth=ncread(fname,'depth',[1 1 1],[Inf Inf 1]); 
  else;
    fname=[pname,'ww3g_',ftime,'-utc_',gnames{i},'.nc'];
  end

  lon_mod=double(ncread(fname,'lon'));
  lat_mod=double(ncread(fname,'lat'));
  [lonm,latm]=meshgrid(lon_mod,lat_mod);

  [dif ilon]=nanmin(abs(lon_mod-lon_obs));
  [dif ilat]=nanmin(abs(lat_mod-lat_obs));
  
  
  %if i==2; return; end

  if i==2;
    pcolor(lon_mod,lat_mod,depth')
    colormap(cmocean('deep'))
    shading flat;
    colorbar
	else;
    %grid contour
    min_lon=min(lon_mod(:)); max_lon=max(lon_mod(:));
    min_lat=min(lat_mod(:)); max_lat=max(lat_mod(:));

    plot([min_lon min_lon],[min_lat max_lat],'color',colors{i},'linewidth',2)
    plot([max_lon max_lon],[min_lat max_lat],'color',colors{i},'linewidth',2)
    plot([min_lon max_lon],[min_lat min_lat],'color',colors{i},'linewidth',2)
    plot([min_lon max_lon],[max_lat max_lat],'color',colors{i},'linewidth',2)
	end

  if i~=2;
    fname=[pname,'ww3p_interp_',ftime,'-utc_',gnames{i},'.nc'];
    lon=ncread(fname,'lon');
  	lat=ncread(fname,'lat');
		plot(lon,lat,'.','color',colors{i},'markersize',4)
  end

end

xlim([170 185]); ylim([-60 -15]) % NZ and AUS
xlim([min(lon_obss)-.9 max(lon_obss)+.9]); % NZ and AUS
ylim([min(lat_obss)-.9 max(lat_obss)+.9]); % NZ and AUS
%xlim([140 210]); ylim([-60 10])  % GLOBAL and TONGA
%xlim([184.1 186.5]); ylim([-22 -15])   % GLOBAL and TONGA
caxis([0 200])

plot_aus=0;
if plot_aus==1;
    fname=['~/ww3/grids/ww3.acs.g3.201612.2017','.nc'];
    lon=ncread(fname,'longitude');
  	lat=ncread(fname,'latitude');
		plot(lon,lat,'.','color','m','markersize',4)
end


plot_obs=1;
if plot_obs==1;
    for i=1:length(lon_obss)
  		plot(lon_obss(i),lat_obss(i),'s','color','r','markersize',4)
  		text(lon_obss(i),lat_obss(i),tstation{i},'fontsize',12,'color','k')%,'markersize',4)
    end
end

%title('Bathymetry = GLOBALWAVE, Red = NZWAVE, and Black = NZWAVE-HR. Points are from p files')

%return

%title('Bathymetry = GLOBALWAVE, Cyan = TONGAWAVE. Points are from ww3 p files')

save_fig=0;
if save_fig==1;
  path_fig=['/scale_wlg_nobackup/filesets/nobackup/niwa03150/santanarc/figures/']; % GLOBALWAVE/'];%2018/01/05/00'];
  export_fig(gcf,[path_fig,'south_wellington_stations'],'-png','-r150');
end


