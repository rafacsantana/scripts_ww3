clear
close all

time=round(now)-4; % look for files generated 4 days ago
time=datenum(2021,1,1,0,0,0);


%Loading OBS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/';
path_obs=[path,'data/obs/'];
%file='SteepHead_SP_FullRecord_QC';
%file_obs=file;
%display(['Loading: ',path,file,'.mat']);
%load([path_obs,file,'.mat'])%,'time','obs')


tstation={'Sinclair Head','Ohau Head','Mana Island','Makara'  ,'Karori Rock','Wairewa Lake Forsyth','Taumutu','Gulf Harbour','Mangawhai'};
stations={'sinclair_head','ohau_head','mana_island','makara'  ,'karori_rock','wairewa_lake_forsyth','Taumutu','Gulf_Harbour','Mangawhai'};
lat_obss=[-41.366885     ,-42.256063 ,-41.083327   ,-41.203393, -41.344623  ,-43.84110,-43.866527,-36.675007,-36.114621];
lon_obss=[174.716761     , 173.853341,174.763960   ,174.702391, 174.650430  ,172.71835,172.381380,174.842037,174.651937];

tstation={tstation{end}};
stations={stations{end}};
scase='manukau_offshore'; %'tairua' % 'port-waikato'; % 'dixon-anderson';
scase='kaitorete'; % 'graveyard_hills' % 'port-waikato'; % 'dixon-anderson';
[tstation,lat_obss,lon_obss]=read_wave_stations(scase);

% grid selection
gsel=3; % grid selected (1=GLOBALWAVE, 2=NZWAVE, 3=NZWAVE-HR, 4=TONGAWAVE)

%tstation={tstation{14}};
%lon_obss= lon_obss(14); 
%lat_obss= lat_obss(14);

%lon_obss=lon_obss(end); lat_obss=lat_obss(end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Model path  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path='/scale_wlg_devoper/filesets/archive/ecoconnect/EFS/';
path='/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/';

ptime=datestr(time,'/YYYY/mm/DD/HH/');
ftime=datestr(time,'YYYYmmDDHH');

grids={'GLOBALWAVE','NZWAVE','NZWAVE-HR','TONGAWAVE'};
gnames={'globalwave+globalum','nzwave+nzlam','nzwave_hr+nzcsm','tongawave+globalum'};

% Loading GRIDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scrsz=[1 1 1366 768];
scrsz=[1 1 1920 1080];
scrsz=[1 1 1910 990];
%scrsz=get(0,'screensize');
figure('position',scrsz,'color',[1 1 1],'visible','on')
hold on
set(gca,'fontsize',12,'fontweight','bold')

colors={'b','r','r','k'};



for i=gsel

  pname=[path,grids{i},ptime];
  %if i==3;
    fname=[pname,'ww3g_',ftime,'-utc_',gnames{i},'.nc'];
    disp(['Readin: ',fname])
    depth=ncread(fname,'depth',[1 1 1],[Inf Inf 1]); 
  %else;
    %fname=[pname,'ww3g_',ftime,'-utc_',gnames{i},'.nc'];
  %end

  lon_mod=double(ncread(fname,'lon'));
  lat_mod=double(ncread(fname,'lat'));
  [lonm,latm]=meshgrid(lon_mod,lat_mod);

  plot_raw_grid=1;

  if plot_raw_grid==1;
   
    if i==gsel;
      inan=isnan(depth); depth(inan)=99999;
      pcolor(lon_mod,lat_mod,depth')
      %[cs,h]=contour(lon_mod,lat_mod,depth',[0:10:50],'k');
      %clabel(cs,h,'color',[0 0 0],'fontsize',12,'fontweight','bold','LabelSpacing',2000)

      colormap(cmocean('deep'))
      shading flat;
      cb=colorbar;
      set(get(cb,'ylabel'),'string','Depth (m)','fontsize',12,'fontweight','bold');
	  else;
      %grid contour
      min_lon=min(lon_mod(:)); max_lon=max(lon_mod(:));
      min_lat=min(lat_mod(:)); max_lat=max(lat_mod(:));

      plot([min_lon min_lon],[min_lat max_lat],'color',colors{i},'linewidth',2)
      plot([max_lon max_lon],[min_lat max_lat],'color',colors{i},'linewidth',2)
      plot([min_lon max_lon],[min_lat min_lat],'color',colors{i},'linewidth',2)
      plot([min_lon max_lon],[max_lat max_lat],'color',colors{i},'linewidth',2)
	  end

    plot_ww3p=0;
    if plot_ww3p==1;
      fname=[pname,'ww3p_',ftime,'-utc_',gnames{i},'.nc'];
      lon=ncread(fname,'lon');
    	lat=ncread(fname,'lat');
	  	plot(lon,lat,'.','color',colors{i},'markersize',4)
    end

    plot_obs=1;
    if plot_obs==1;
       for i=1:length(lon_obss)
         [dif ilon]=nanmin(abs(lon_mod-lon_obss(i)));
         [dif ilat]=nanmin(abs(lat_mod-lat_obss(i)));
       	 %plot(lon_mod(ilon),lat_mod(ilat),'s','color','b','markersize',10)
       	 plot(lon_obss(i),lat_obss(i)    ,'.','color','r','markersize',20)
       	 text(lon_obss(i)+.05,lat_obss(i),[replace(tstation{i},'_',' '),' ',num2str(depth(ilon,ilat),'%.2f'),' m'],'fontsize',12,'color','k')%,'markersize',4)
       end
       dx=.5; dy=.5;
       xlim([min(lon_obss)-dx max(lon_obss)+dx]); ylim([min(lat_obss)-dy max(lat_obss)+dy]); % NZ and AUS
       xlim([min(lon_obss)-dx max(lon_obss)+dx*4]); ylim([min(lat_obss)-dy max(lat_obss)+dy]); % NZ and AUS
    end

  else

    if i==gsel;

      gx = geoaxes;
      geobasemap(gx,'none')% 'grayland')%'streets-light'); %
      pcolor(lon_mod,lat_mod,depth')

return
      %'darkwater', 'grayland', 'bluegreen', 'grayterrain', 'colorterrain', 'landcover', 'streets', 'streets-light', 'streets-dark', 'satellite', 'topographic', 'none'
      geolimits(gx,[-45.2 -41.6],[171 174.2])
      set(gx,'fontsize',12,'fontweight','bold')
      hold on
      geoscatter(gx,clatlon(:,1),clatlon(:,2),'Marker','o','MarkerFaceColor',[.7 0 .7], 'MarkerEdgeColor','black')%,'MarkerSize',10);
      %geoplot(clatlon(:,1),clatlon(:,2),'Marker','o','MarkerFaceColor','magenta', 'MarkerEdgeColor','black','MarkerSize',10);
      geoscatter(gx,lat_mod(points(:,2)),lon_mod(points(:,1)),'Marker','o','MarkerFaceColor',[0 .7 0], 'MarkerEdgeColor','black');


      pcolor(lon_mod,lat_mod,depth')
      [cs,h]=contour(lon_mod,lat_mod,depth',[0:10:50],'k');
      clabel(cs,h,'color',[0 0 0],'fontsize',12,'fontweight','bold','LabelSpacing',2000)

      colormap(cmocean('deep'))
      shading flat;
      cb=colorbar;
      set(get(cb,'ylabel'),'string','Depth (m)','fontsize',12,'fontweight','bold');
	  else;
      %grid contour
      min_lon=min(lon_mod(:)); max_lon=max(lon_mod(:));
      min_lat=min(lat_mod(:)); max_lat=max(lat_mod(:));

      plot([min_lon min_lon],[min_lat max_lat],'color',colors{i},'linewidth',2)
      plot([max_lon max_lon],[min_lat max_lat],'color',colors{i},'linewidth',2)
      plot([min_lon max_lon],[min_lat min_lat],'color',colors{i},'linewidth',2)
      plot([min_lon max_lon],[max_lat max_lat],'color',colors{i},'linewidth',2)
	  end

 
  end % if plot_raw_grid==1;

end

%xlim([170 185]); ylim([-60 -15]) % NZ and AUS
%xlim([140 210]); ylim([-60 10])  % GLOBAL and TONGA
%xlim([184.1 186.5]); ylim([-22 -15])   % GLOBAL and TONGA
%xlim([171.5 173.3]); ylim([-44.7 -43])   % focus on lake ellesmere
caxis([0 3000])

plot_aus=0;
if plot_aus==1;
    fname=['~/ww3/grids/ww3.acs.g3.201612.2017','.nc'];
    lon=ncread(fname,'longitude');
  	lat=ncread(fname,'latitude');
		plot(lon,lat,'.','color','m','markersize',4)
end



%title('Bathymetry = GLOBALWAVE, Red = NZWAVE, and Black = NZWAVE-HR. Points are from p files')

%return

title('Locations of analysis and model bathymetry')
ylabel('Latitude')
xlabel('Longitude')

save_fig=0;
if save_fig==1;
  path_fig=['/scale_wlg_nobackup/filesets/nobackup/niwa03150/santanarc/figures/']; % GLOBALWAVE/'];%2018/01/05/00'];
  display(['Saving: ',path_fig,scase]);
  export_fig(gcf,[path_fig,scase],'-png','-r150');
  saveas(gcf,[path_fig,scase],'fig');
end


