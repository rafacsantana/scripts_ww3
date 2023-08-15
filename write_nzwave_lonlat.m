clear
close all

time=round(now)-4; % look for files from 4 days in the past
                  % datenum(2023,7,18,0,0,0);

path='/scale_wlg_devoper/filesets/archive/ecoconnect/EFS/';

ptime=datestr(time,'/YYYY/mm/DD/HH/');
ftime=datestr(time,'YYYYmmDDHH');

grids={'GLOBALWAVE','NZWAVE','NZWAVE-HR'};
gnames={'globalwave+globalum','nzwave+nzlam','nzwave_hr+nzcsm'};

%ww3g_2023071800-utc_globalwave+globalum.nc
%ww3g_2023071800-utc_nzwave+nzlam.nc
%ww3g_2023071800-utc_nzwave_hr+nzcsm.nc

scrsz=[1 1 1366 768];
%scrsz=get(0,'screensize');
figure('position',scrsz,'color',[1 1 1],'visible','on')
hold on
set(gca,'fontsize',12,'fontweight','bold')

colors={'b','r','k'};

for i=2

  pname=[path,grids{i},ptime]
  if i==1;
    fname=[pname,'ww3g_',ftime,'-utc_',gnames{i},'.nc'];
    depth=ncread(fname,'depth',[1 1 1],[Inf Inf 1]); 
  else;
    fname=[pname,'ww3g_',ftime,'-utc_',gnames{i},'.nc'];
  end
  lon=ncread(fname,'lon');
  lat=ncread(fname,'lat');

  fda=['WW3_grid_',grids{i},'.nc'];
  if exist(fda)==2
     system(['rm -rf ',fda]);
  end

  display(['Writing: lonlat ',fda])
  % creating nc and variables
  nccreate(fda,'lon','Dimensions', {'lon',  length(lon)}); 
  nccreate(fda,'lat','Dimensions', {'lat',  length(lat)}); 
  ncwrite(fda,'lon',lon)
  ncwrite(fda,'lat',lat)

  close all
  if i==2; return; end

  if i==1;
    pcolor(lon,lat,depth')
    colormap(cmocean('deep'))
    shading flat;
    colorbar
    xlim([140 185])
    ylim([-60 -15])
	else;
    %grid contour
    min_lon=min(lon(:)); max_lon=max(lon(:));
    min_lat=min(lat(:)); max_lat=max(lat(:));

    plot([min_lon min_lon],[min_lat max_lat],'color',colors{i},'linewidth',2)
    plot([max_lon max_lon],[min_lat max_lat],'color',colors{i},'linewidth',2)
    plot([min_lon max_lon],[min_lat min_lat],'color',colors{i},'linewidth',2)
    plot([min_lon max_lon],[max_lat max_lat],'color',colors{i},'linewidth',2)
	end

    fname=[pname,'ww3p_',ftime,'-utc_',gnames{i},'.nc'];
    lon=ncread(fname,'lon');
  	lat=ncread(fname,'lat');
		plot(lon,lat,'.','color',colors{i},'markersize',4)

end

title('Bathymetry = GLOBALWAVE, Black = NZWAVE, and Red = NZWAVE-HR. Points are from p files')

export_fig(gcf,'ww3_2023_grids','-png','-r150');

return














