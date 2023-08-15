clear
close all

%Time
start_day  =5;
start_month=1;
start_year =2018;
end_day    =20;
end_month  =1;
end_year   =2018;

% Daily time
time_lima=datenum(start_year,start_month,start_day,0,0,0):1:datenum(end_year,end_month,end_day,0,0,0);

path='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/';
path_mod=[path,'model/'];
path=[path,'data/obs/'];


path_expt=['/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/GLOBALWAVE/'];%2018/01/05/00'];
file_prepath_expt=['ww3g_'];%2018010500-utc_globalwave+globalum.nc']
file_pos=['-utc_globalwave+globalum.nc'];

gnames={'globalwave+globalum','nzwave+nzlam','nzwave_hr+nzcsm','tongawave+globalum'};


proc_data=0;
save_fig =0;

% OBS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file='SteepHead_SP_FullRecord_QC';
file_obs=file;
%    {'date'}    {'PeakPerid'}    {'PeakDirection'}    {'Directional spread'}    {'Tm01'}    {'Tm02'}    {'Hs'}    {'Qp'}

if proc_data==1
  
  display(['Processing: ',path,file,'.txt']);
  ob=importdata([path,file,'.txt']);
  ob.textdata(1,:)
  for i=2:length(ob.textdata)
    time_obs(i-1)=datenum(ob.textdata{i,1},'YYYY-mm-ddTHH:MM:SS');
  end
  obs=ob.data;
  lat_obs=-43+(45/60); lon_obs=173+(20/60);
  save([path,file,'.mat'],'time_obs','obs','lat_obs','lon_obs')

else

  display(['Loading: ',path,file,'.mat']);
  load([path,file,'.mat'])%,'time','obs')
  lat_obs=-43+(45/60); lon_obs=173+(20/60);

end


path_dm=[path_expt,'matlab/'];
if exist([path_dm,datestr(time_lima(1),'YYYYmmDDHH'),'_',datestr(time_lima(end),'YYYYmmDDHH'),'.mat'])~=2

  hs_mod=[];tp_mod=[];
  pd_mod=[];time_mod=[];
  
  % Loop in model time
  for t=time_lima
  
    ptime=datestr(t,'YYYY/mm/DD/HH/');
    ftime=datestr(t,'YYYYmmDDHH');
  
    pname=[path_expt,ptime];
  
    %if strncmp('GLOBAL')
    %end
  
    fname=[pname,'ww3g_',ftime,'-utc_',gnames{1},'.nc'];
    display(['Loading: ',fname]);
    if t==time_lima(1)
      lon_mod=ncread(fname,'lon');
      lat_mod=ncread(fname,'lat');
      [dif ilon]=nanmin(abs(lon_mod-lon_obs));
      [dif ilat]=nanmin(abs(lat_mod-lat_obs));
      if strcmp(file_obs,'SteepHead_SP_FullRecord_QC')
        ilon=ilon+3;
      end
    end
  
    time=squeeze(double(ncread(fname,'time')))./24+t; time=time(1:end-1);
    hs=squeeze(double(ncread(fname,'hsig',[ilon ilat 1],[1 1 Inf]))); hs=hs(1:end-1);
    tp=squeeze(double(ncread(fname,'tpeak',[ilon ilat 1],[1 1 Inf]))); tp=tp(1:end-1);
    pd=squeeze(double(ncread(fname,'peak_direction',[ilon ilat 1],[1 1 Inf]))); pd=pd(1:end-1);
  
    time_mod=[time_mod;time];
    hs_mod=[hs_mod;hs];
    tp_mod=[tp_mod;tp];
    pd_mod=[pd_mod;pd];
  
  end
  
  model.data(:,6)=hs_mod;
  model.data(:,1)=tp_mod;
  model.data(:,2)=pd_mod;
  
  path_dm=[path_expt,'matlab/'];
  system(['mkdir -p ',path_dm]);
  save([path_dm,datestr(time_lima(1),'YYYYmmDDHH'),'_',datestr(time_lima(end),'YYYYmmDDHH'),'.mat'],'time_mod','model')

else

  load([path_dm,datestr(time_lima(1),'YYYYmmDDHH'),'_',datestr(time_lima(end),'YYYYmmDDHH'),'.mat'])%,'time_mod','model')

end



scrsz=[1 1 1366 768];
%scrsz=get(0,'screensize');
figure('position',scrsz,'color',[1 1 1],'visible','on')
hold on
set(gca,'fontsize',12,'fontweight','bold')

dcol=[1,4,5]; %[6,1,2]; % Hs, Tp, dir
dcol=[6,1,2]; % Hs, Tp, dir
tnames={'Peak Period',    'Peak Direction',    'Directional spread',   'Tm01','Tm02','Hs','Qp'};
for i=1:length(dcol)

  subplot(3,1,i);
  set(gca,'fontsize',12,'fontweight','bold')
  hold on
  plot(time_obs,obs(:,dcol(i))','.g','linewidth',2)
  plot(time_mod,model.data(:,dcol(i))','.k','linewidth',2)

  xlim([time_lima(1)-2 time_lima(end)+2])
  %if i==3
    datetick('x','dd/mmm/yyyy','keeplimits')
  %end
  title(tnames{dcol(i)})
  grid

end


if save_fig==1
  export_fig(gcf,'steephead_period','-png','-r150');
end


