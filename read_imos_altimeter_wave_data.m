%function [lon,lat,time,swh]=read_imos_altimeter_wave_data(fnumber);

clear; close all


path_alt='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/swan/projects/VANUATU/data/satellite_wave_data/IMOS_-_SRS_Surface_Waves_Sub-Facility_-_altimeter_wave_wind_source_files_smaller_area';

falt=dir(path_alt);

for i=[5,9] % length(falt)-1]
  
  fname=[path_alt,'/',falt(i).name];

  display(['Reading: ',falt(i).name])

  clear swh

  swh(:,1)=ncread([fname],'SWH_C');
  swh(:,2)=ncread([fname],'SWH_KU');
  swh(:,3)=ncread([fname],'SWH_KU_CAL');
  
  time=ncread([fname],'TIME')+datenum('1985-01-01 00:00:00');
  lon=ncread([fname],'LONGITUDE');
  lat=ncread([fname],'LATITUDE');


  scrsz=[1 10 1900 990];
  figure('position',scrsz,'color',[1 1 1],'visible','on')
  hold on
  set(gca,'fontsize',12,'fontweight','bold')
  
  plot(time,swh,'.','markersize',20)

  legend('SWH C','SWH KU','SWH KU/KA CAL')

  xlim([datenum(2020,4,1) datenum(2020,4,11)])

  datetick('x','dd-mm','keeplimits')

  title(replace(falt(i).name,'_','-'))

  [dif loci]=nanmin(abs(time-datenum(2020,4,4)))
  [dif locf]=nanmin(abs(time-datenum(2020,4,6)))

end

%end

