clear
%close all
%clf('reset')

time=datenum(2020,4,1,0,0,0):1/24:datenum(2020,4,11,0,0,0);

%output
expt='lv_era5_tc_harold_corr/'; % 3';
file='era5_tc_harold_corr.dat';
swanInputFilePath=['/scale_wlg_persistent/filesets/project/niwa03150/santanarc/swan/projects/VANUATU/',expt,file];  

% input winds
filewnd=('/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/ERA5/ERA5_u10_y2020.nc');
lonwnd =double(ncread(filewnd,'longitude'));
latwnd =double(ncread(filewnd,'latitude'));
timewnd=double(ncread(filewnd,'time'))./(24)+datenum(1900,1,1); %:1/48:datenum(2000,01,18,0,0,0);

loni=163; lonf=171;
lati=-19; latf=-12;

[dif loi]=nanmin(abs(lonwnd-loni));
[dif lof]=nanmin(abs(lonwnd-lonf));
[dif lai]=nanmin(abs(latwnd-latf));
[dif laf]=nanmin(abs(latwnd-lati));

[dif ti]=nanmin(abs(timewnd-time(1)));
[dif tf]=nanmin(abs(timewnd-time(end)));

uwnd=squeeze((ncread(filewnd,'u10', [loi lai ti],[lof-loi+1 laf-lai+1 tf-ti+1])));

filewnd=('/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/ERA5/ERA5_v10_y2020.nc');
vwnd=squeeze((ncread(filewnd,'v10', [loi lai ti],[lof-loi+1 laf-lai+1 tf-ti+1])));

latwnd=latwnd(lai:laf);
lonwnd=lonwnd(loi:lof);
timewnd=timewnd(ti:tf);

return

% flipping latitude so it starts at -90
latwnd=flipud(latwnd);
% flipping u and v according to latitude
uwnd=fliplr(uwnd);
vwnd=fliplr(vwnd);

% correct era5 winds to match max WMO wind speed
% ibtracks/TC_Harold_track_max_speed.csv
ibfile=['/scale_wlg_persistent/filesets/project/niwa03150/santanarc/swan/projects/VANUATU/data/ibtracks/TC_Harold_track_max_speed.csv'];
ob=importdata(ibfile);
for i=3:length(ob.textdata)
  time_tc(i-2)=datenum(ob.textdata{i,1},'dd/mm/YYYY HH:MM');
end
maxwnd=ob.data./1.94384; % converting knots to m/s

for i=1:length(timewnd)

  display(['Processing wind',' at ',datestr(timewnd(i),'yyyymmddHHMM')])
  magwnd=sqrt(squeeze(uwnd(:,:,i)).^2+squeeze(vwnd(:,:,i)).^2);
  [maxera5]=nanmax(magwnd(:));
  [dif itc]=nanmin(abs(time_tc-timewnd(i)));
  cfactor=maxwnd(itc)./maxera5;

  if cfactor>1 
    display(['Correcting winds by factor: ',num2str(cfactor),' at ',datestr(timewnd(i),'yyyymmddHHMM')])
    magwnd=magwnd.*cfactor;
  end

  dwnd=atand(vwnd(:,:,i)./uwnd(:,:,i));
  uwnd(:,:,i)=magwnd.*cosd(dwnd);
  vwnd(:,:,i)=magwnd.*sind(dwnd);

end

% Open the SWAN input file for writing
fid = fopen(swanInputFilePath, 'w');

for k=1:length(timewnd)
  display(['Writing: ',datestr(timewnd(k),'yyyymmddHHMM')])
  for j=1:length(latwnd)
    % uwnd
    for i=1:length(lonwnd)
      if i==length(lonwnd)
       fprintf(fid, '%.2f\n', squeeze(uwnd(i,j,k)));
      else
        fprintf(fid, '%.2f\t', squeeze(uwnd(i,j,k)));
      end
    end
  end
  for j=1:length(latwnd)
    % vwnd
    for i=1:length(lonwnd)
      if i==length(lonwnd)
       fprintf(fid, '%.2f\n', squeeze(vwnd(i,j,k)));
      else
        fprintf(fid, '%.2f\t', squeeze(vwnd(i,j,k)));
      end
    end
  end
end

% Close the SWAN input file
fclose(fid);

fprintf('SWAN input file created: %s\n', swanInputFilePath);


return
