clear
close all

fwave='/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/NZWAVE-HR/2023/01/01/00/ww3g_2023010100-utc_nzwave_hr+nzcsm.nc';

fccam='/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/projections/historical/hourly/uas/uas_historical_GFDL-ESM4_CCAM_hourly_NZ12km_raw.nc';

fgccam='/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/projections/historical/6hrly_global/ua/ua_historical_GFDL-ESM4_CCAM_6hrly_Global_raw.nc';

lonw=ncread(fwave,'lon');
latw=ncread(fwave,'lat');

lonc=ncread(fccam,'lon');
latc=ncread(fccam,'lat');

long=ncread(fgccam,'lon');
latg=ncread(fgccam,'lat');

display(['NZWAVEHR corners: minlon=',num2str(min(lonw)),' minlat=',num2str(min(latw)),' maxlon=',num2str(max(lonw)),' maxlat=',num2str(max(latw))])

display(['NZCCAM corners: minlon=',num2str(min(lonc)),' minlat=',num2str(min(latc)),' maxlon=',num2str(max(lonc)),' maxlat=',num2str(max(latc))])

