clear 
close all



path_bath='/scale_wlg_persistent/filesets/project/niwa00020/data/topography/dtm/';
file_bath='nzbath250_2016.nc';

lat_obs=-41.4022; lon_obs=174.84669;

% full doamin map 
loni=lon_obs-.1; 
lonf=lon_obs+.1; 
lati=lat_obs-.1; 
latf=lat_obs+.1; 

bath_nc=[path_bath,file_bath];
lon_bath=double(ncread(bath_nc,'lon'));
lat_bath=double(ncread(bath_nc,'lat'));
[difs ilonb]=nanmin(abs(lon_bath-loni)); 
[difs flonb]=nanmin(abs(lon_bath-lonf)); 
[difs ilatb]=nanmin(abs(lat_bath-lati)); 
[difs flatb]=nanmin(abs(lat_bath-latf));
lon_bath=lon_bath(ilonb:flonb);
lat_bath=lat_bath(ilatb:flatb);
bath=ncread(bath_nc,'height',[ilonb ilatb],[flonb-ilonb+1 flatb-ilatb+1]);
%bath=double(ncread(bath_nc,'height'));
%bath=bath(ilonb:flonb,ilatb:flatb);
[lon_bathm,lat_bathm]=meshgrid(lon_bath,lat_bath);
lon_bathm=lon_bathm'; lat_bathm=lat_bathm';
      


figure
hold on
pcolor(lon_bath,lat_bath,bath')
plot(lon_obs,lat_obs,'.r')
contour(lon_bath,lat_bath,bath',[0 -20 -100],'k')
shading flat
colorbar
caxis([-100 10])
