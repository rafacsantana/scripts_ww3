clear
%close all
%clf('reset')

time=datenum(2000,1,1,0,0,0):.5/24:datenum(2000,01,02,0,0,0);

expt='luganville_3'; % 3';



% winds
filewnd=('/scale_wlg_persistent/filesets/project/niwa03150/santanarc/delft3d/projects/VANUATU/data/TCwindgen/test_LV.nc');
lonwnd =double(ncread(filewnd,'x'));
latwnd =double(ncread(filewnd,'y'));
timewnd=double(ncread(filewnd,'time'))./(24*3600)+datenum(2000,01,01,0,0,0); %:1/48:datenum(2000,01,18,0,0,0);

uwnd=double(ncread(filewnd,'U')); % ,[1 1 iw],[Inf Inf 1])');
vwnd=double(ncread(filewnd,'V')); % ,[1 1 iw],[Inf Inf 1])');

swanInputFilePath='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/delft3d/projects/VANUATU/luganville_3/test_LV_matlab5.dat';

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
