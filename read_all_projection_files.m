clear 
close all
%bsose  cesm  cesm2  glorys  isas20_argo

path={
'/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/projections'
};


p1={'historical','ssp370'};
p2={'6hrly_global','hourly'};
p3={'psl','ua','va','uas','vas','sic'};


vars={
'time','time','time','time','time',...
'time','time','time','time',...
'time','time',...
'time_bound','time_bound','time','time'};

for i=1:length(p1)
for ii=1:length(p2)
for iii=1:length(p3)

  paths=strcat(path,'/',p1{i},'/',p2{ii},'/',p3{iii},'/');

  paths=paths{1};
  display(['Checking: ',paths]);
  files=dir(paths);

  for j=1:length(files)

    fname=files(j).name;
    fbytes=files(j).bytes;
    if length(fname)>5

      if strcmp(fname(end-2:end),'.nc') & fbytes > 1
        display(['Reading: ',paths,fname])
        t=ncread([paths,fname],'time'); %,vars{i});
        %ncdisp([paths,fname]) %,vars{i});
        fclose('all');
      end

    end

  end

end
end
end

