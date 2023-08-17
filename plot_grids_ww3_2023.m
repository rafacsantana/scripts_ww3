clear
close all

time=round(now)-4; % look for files generated 4 days ago
                  % datenum(2023,7,18,0,0,0);

%Loading OBS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/';
path_obs=[path,'data/obs/'];
file='SteepHead_SP_FullRecord_QC';
file_obs=file;
display(['Loading: ',path,file,'.mat']);
load([path_obs,file,'.mat'])%,'time','obs')


path='/scale_wlg_devoper/filesets/archive/ecoconnect/EFS/';
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
  if i==1;
    fname=[pname,'ww3g_',ftime,'-utc_',gnames{i},'.nc'];
    depth=ncread(fname,'depth',[1 1 1],[Inf Inf 1]); 
  else;
    fname=[pname,'ww3g_',ftime,'-utc_',gnames{i},'.nc'];
  end

  lon_mod=double(ncread(fname,'lon'));
  lat_mod=double(ncread(fname,'lat'));
  [lonm,latm]=meshgrid(lon_mod,lat_mod);
  %if i==2; return; end

  if i==1;
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

  if i~=1;
    fname=[pname,'ww3p_interp_',ftime,'-utc_',gnames{i},'.nc'];
    lon=ncread(fname,'lon');
  	lat=ncread(fname,'lat');
		plot(lon,lat,'.','color',colors{i},'markersize',4)
  end

end

xlim([170 185]); ylim([-60 -15]) % NZ and AUS
%xlim([140 210]); ylim([-60 10])  % GLOBAL and TONGA
%xlim([184.1 186.5]); ylim([-22 -15])   % GLOBAL and TONGA

plot_aus=0;
if plot_aus==1;
    fname=['~/ww3/grids/ww3.acs.g3.201612.2017','.nc'];
    lon=ncread(fname,'longitude');
  	lat=ncread(fname,'latitude');
		plot(lon,lat,'.','color','m','markersize',4)
end

plot_obs=1;
if plot_obs==1;
		plot(lon_obs,lat_obs,'s','color','g','markersize',4)
end

%title('Bathymetry = GLOBALWAVE, Red = NZWAVE, and Black = NZWAVE-HR. Points are from p files')

%return

title('Bathymetry = GLOBALWAVE, Cyan = TONGAWAVE. Points are from ww3 p files')

save_fig=1;
if save_fig==1;
  export_fig(gcf,'ww3_2023_tonga_grid','-png','-r150');
end












