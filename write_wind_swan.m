clear
%close all
%clf('reset')


expt='luganville_3'; % 3';
expt='lv_tcwindgen_tc_harold'; % 3';


if strcmp(expt,'luganville_3')
  time=datenum(2000,1,1,0,0,0):.5/24:datenum(2000,01,02,0,0,0);
elseif strcmp(expt,'lv_tcwindgen_tc_harold')
  time=datenum(2020,4,1,0,0,0):.5/24:datenum(2020,4,11,0,0,0);
end

% winds
if strcmp(expt,'luganville_3')
  filewnd=('/scale_wlg_persistent/filesets/project/niwa03150/santanarc/swan/projects/VANUATU/data/TCwindgen/test_LV.nc');
elseif strcmp(expt,'lv_tcwindgen_tc_harold')
  filewnd=('/scale_wlg_persistent/filesets/project/niwa03150/santanarc/TCwindgen_zhonghou/projects/tc_harold_windgen_smaller/tc_harold_windgen.nc');
end
lonwnd =double(ncread(filewnd,'x'));
latwnd =double(ncread(filewnd,'y'));
timewnd=double(ncread(filewnd,'time'))./(24*3600)+time(1); %:1/48:datenum(2000,01,18,0,0,0);

%return

uwnd=double(ncread(filewnd,'U')); % ,[1 1 iw],[Inf Inf 1])');
vwnd=double(ncread(filewnd,'V')); % ,[1 1 iw],[Inf Inf 1])');

if strcmp(expt,'luganville_3')
  swanInputFilePath='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/swan/projects/VANUATU/luganville_3/test_LV_matlab5.dat';
elseif strcmp(expt,'lv_tcwindgen_tc_harold')
  swanInputFilePath='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/swan/projects/VANUATU/lv_tcwindgen_tc_harold/tc_harold_windgen_matlab.dat';
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
