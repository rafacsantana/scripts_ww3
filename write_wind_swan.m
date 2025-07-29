clear
close all
%clf('reset')


expt='luganville_3'; % 3';
expt='lv_tcwindgen_tc_harold_cfac'; % 3';
expt='lv_tcwindgen_tc_harold_xuz'; % 3';


if strcmp(expt,'luganville_3')
  time=datenum(2000,1,1,0,0,0):.5/24:datenum(2000,01,02,0,0,0);
elseif strncmp(expt,'lv_',3)
  time=datenum(2020,4,1,0,0,0):.5/24:datenum(2020,4,11,0,0,0);
end

% winds
if strcmp(expt,'luganville_3')
  filewnd=('/scale_wlg_persistent/filesets/project/niwa03150/santanarc/swan/projects/VANUATU/data/TCwindgen/test_LV.nc');
elseif strncmp(expt,'lv_',3)
  %filewnd=strcat('/scale_wlg_persistent/filesets/project/niwa03150/santanarc/TCwindgen_zhonghou/projects/',expt,'/tc_harold_windgen.nc');
  filewnd=strcat('/scale_wlg_persistent/filesets/project/niwa03150/santanarc/swan/projects/VANUATU/',expt,'/Luganville_harold.nc');
end
lonwnd =double(ncread(filewnd,'x'));
latwnd =double(ncread(filewnd,'y'));
timewnd=double(ncread(filewnd,'time'))./(24*3600)+time(1); %:1/48:datenum(2000,01,18,0,0,0);

%return

uwnd=double(ncread(filewnd,'U')); % ,[1 1 iw],[Inf Inf 1])');
vwnd=double(ncread(filewnd,'V')); % ,[1 1 iw],[Inf Inf 1])');

uwnd=permute(uwnd,[2 1 3]);
vwnd=permute(vwnd,[2 1 3]);

%return

corr_wind=1;
if corr_wind==1;
  cfactor=0.93;
  magwnd=sqrt(uwnd.^2+vwnd.^2);
  magwnd=magwnd.*cfactor;
  dwnd=atan2d(vwnd,uwnd);
  %dwnd=mod(-90-dwnd,360);

  uwnd=magwnd.*cosd(dwnd);
  vwnd=magwnd.*sind(dwnd);
end

uwnd=permute(uwnd,[2 1 3]);
vwnd=permute(vwnd,[2 1 3]);

%it=find(timewnd==datenum(2020,4,4,6,30,0)); it=it(1);
%vs=40;
%figure; quiver(lonwnd(1:vs:end),latwnd(1:vs:end),uwnd((1:vs:end),(1:vs:end),it),vwnd((1:vs:end),(1:vs:end),it));
%return

if strcmp(expt,'luganville_3')
  swanInputFilePath='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/swan/projects/VANUATU/luganville_3/test_LV_matlab5.dat';
elseif strncmp(expt,'lv_',3)
  swanInputFilePath=strcat('/scale_wlg_persistent/filesets/project/niwa03150/santanarc/swan/projects/VANUATU/',expt,'/tc_harold_windgen_matlab.dat');
end

if exist(swanInputFilePath)==2
   system(['rm -rf ',swanInputFilePath]);
end


% Open the SWAN input file for writing
fid = fopen(swanInputFilePath, 'w');

it=find(timewnd==timewnd(end));
for k=1:it(1)
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
