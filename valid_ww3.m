clear
close all

% Daily time
time_lima=datenum(2014,7,9,0,0,0):1:datenum(2015,5,21,0,0,0); % initial date
%time_lima=datenum(2014,8,5,0,0,0):1:datenum(2014,8,17,0,0,0); % 6-m wave event
%time_lima=datenum(2014,8,5,0,0,0):1:datenum(2014,8,6,0,0,0); % Firts two days of output

runs=[1,2];

expt_names={'GLOBALWAVE','NZWAVE'};

set_xlim  =1;
proc_obs  =0;
save_fig  =0;
check_file=0;

ww3pre={'ww3p_interp','ww3g'};
ww3pre=ww3pre{2};



for i=1:length(runs)
  expts{i}=expt_names{runs(i)};
end

path_source=['/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/']; % GLOBALWAVE/'];%2018/01/05/00'];
path_santanarc='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/';
path_obs=[path_santanarc,'data/obs/'];


% OBS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file='SteepHead_SP_FullRecord_QC';
file_obs=[path_obs,file];
%    {'date'}    {'PeakPerid'}    {'PeakDirection'}    {'Directional spread'}    {'Tm01'}    {'Tm02'}    {'Hs'}    {'Qp'}

if proc_obs==1
  
  display(['Processing: ',file_obs,'.txt']);
  ob=importdata([file_obs,'.txt']);
  ob.textdata(1,:)
  for i=2:length(ob.textdata)
    time_obs(i-1)=datenum(ob.textdata{i,1},'YYYY-mm-ddTHH:MM:SS');
  end
  time_obs=time_obs-.5;
  obs=ob.data;
  lat_obs=-43.7567 % from ECAN website -43+(45/60);
  lon_obs=173.3358 % 173+(20/60);
  save([file_obs,'.mat'],'time_obs','obs','lat_obs','lon_obs')

else

  display(['Loading: ',file_obs,'.mat']);
  load([file_obs,'.mat'])%,'time','obs')

end

scrsz=[1 1 1366 768];
scrsz=get(0,'screensize');

figure('position',scrsz,'color',[1 1 1],'visible','on')
hold on
set(gca,'fontsize',12,'fontweight','bold')
  
colors={'r','k'};
legh={'Obs'}; legt={'Obs'};

ke=0;
for expt=expts
  ke=ke+1;

  expt=expt{1};

  gnames={'globalwave+globalum','nzwave+nzlam','nzwave_hr+nzcsm','tongawave+globalum'};

  if strncmp(expt,'GLOBAL',6)
    ig=1;
  elseif strncmp(expt,'NZWAVE-HR',9)
    ig=3;
  elseif strncmp(expt,'NZWAVE',6)
    ig=2;
  end

  path_expt=[path_source,expt,'/']; % GLOBALWAVE/'];%2018/01/05/00'];
  
  path_dm=[path_expt,'matlab/'];

  if check_file==1 && exist([path_dm,datestr(time_lima(1),'YYYYmmDDHH'),'_',datestr(time_lima(end),'YYYYmmDDHH'),'.mat'])==2
  
    display(['Loading: ',path_dm,datestr(time_lima(1),'YYYYmmDDHH'),'_',datestr(time_lima(end),'YYYYmmDDHH'),'.mat'])%,'time_mod','model')
    load([path_dm,datestr(time_lima(1),'YYYYmmDDHH'),'_',datestr(time_lima(end),'YYYYmmDDHH'),'.mat'])%,'time_mod','model')
  
  else
  
    hs_mod=[];tp_mod=[];
    pd_mod=[];time_mod=[];
    
    % Loop in model time
    for t=time_lima
    
      ptime=datestr(t,'YYYY/mm/DD/HH/');
      ftime=datestr(t,'YYYYmmDDHH');
    
      pname=[path_expt,ptime];
    
      fname=[pname,ww3pre,'_',ftime,'-utc_',gnames{ig},'.nc'];
      %fname=[pname,'ww3g_',ftime,'-utc_',gnames{1},'.nc'];
      display(['Loading: ',fname]);

      if t==time_lima(1)
        lon_mod=double(ncread(fname,'lon'));
        lat_mod=double(ncread(fname,'lat'));
        [dif ilon]=nanmin(abs(lon_mod-lon_obs));
        [dif ilat]=nanmin(abs(lat_mod-lat_obs));
        if strcmp(file_obs,'SteepHead_SP_FullRecord_QC')
          %ilon=ilon+3;
        end
        dis=sqrt((lon_mod(ilon)-lon_obs).^2 + (lat_mod(ilat)-lat_obs).^2)*110;
        display(['Distance between obs and grid point is: ',num2str(dis),' km'])
      end
      
      time=squeeze(double(ncread(fname,'time')))./24+t; time=time(1:end-1);
      if strncmp(ww3pre,'ww3p',4)
        hs=squeeze(double(ncread(fname,'hsig',          [ilon 1],[1 Inf]))); hs=hs(1:end-1)';
        tp=squeeze(double(ncread(fname,'tpeak',         [ilon 1],[1 Inf]))); tp=tp(1:end-1)';
        pd=squeeze(double(ncread(fname,'peak_direction',[ilon 1],[1 Inf]))); pd=pd(1:end-1)';
      elseif strncmp(ww3pre,'ww3g',4)
        hs=squeeze(double(ncread(fname,'hsig',          [ilon ilat 1],[1 1 Inf]))); hs=hs(1:end-1);
        tp=squeeze(double(ncread(fname,'tpeak',         [ilon ilat 1],[1 1 Inf]))); tp=tp(1:end-1);
        pd=squeeze(double(ncread(fname,'peak_direction',[ilon ilat 1],[1 1 Inf]))); pd=pd(1:end-1);
      end

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
  
  end



  dcol=[1,4,5]; %[6,1,2]; % Hs, Tp, dir
  dcol=[6,1,2]; % Hs, Tp, dir
  tnames={'Peak Period',    'Peak Direction',    'Directional spread',   'Tm01','Tm02','Significant Wave Height (Hs)','Qp'};
  for i=1:length(dcol)
  
    subplot(3,1,i);
    set(gca,'fontsize',12,'fontweight','bold')
    hold on
    if ke==1
      plot(time_obs,obs(:,dcol(i))','.','color',[0 .7 0],'linewidth',2)
    end
    plot(time_mod,model.data(:,dcol(i))','.','color',colors{ig},'linewidth',2)
  
    if set_xlim==1;
      xlim([time_lima(1)-2 time_lima(end)+2])
    end
  
    %if i==3
      datetick('x','dd/mmm/yyyy','keeplimits')
    %end
    title(tnames{dcol(i)})
    grid('on')

    modeli=interp1(time_mod,model.data(:,dcol(i)),time_obs);

    if dcol(i)==6
      legh=[legh,{[expt,' rmse = ',num2str(nanrmse(obs(:,dcol(i)),modeli'),'%.2f')]}];
    end

  end

end


subplot(3,1,find(dcol==6));
legend([legh],'location','best')


if save_fig==1
  export_fig(gcf,'steephead_period','-png','-r150');
end


