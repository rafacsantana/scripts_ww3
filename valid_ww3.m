clear
close all
warning off

% Daily time
time_lima=datenum(2014,7,9,0,0,0):1:datenum(2014,7,14,0,0,0); % initial date
time_lima=datenum(2014,7,9,0,0,0):1:datenum(2018,1,6,0,0,0); % Graham Harrington ECAN request on 02/10/2023
%time_lima=datenum(2014,8,5,0,0,0):1:datenum(2014,8,17,0,0,0); % 9-m wave event
%time_lima=datenum(2018,1,4,0,0,0):1:datenum(2018,1,10,0,0,0); % crossing between hindcast and ecoconnect
time_lima=datenum(2021,5,15,0,0,0):1:datenum(2021,7,15,0,0,0); % Firts two days of output
%time_lima=datenum(2021,5,15,0,0,0):1:datenum(2021,5,16,0,0,0); % Firts two days of output
%time_lima=datenum(2021,5,15,0,0,0):1:datenum(2021,7,2,0,0,0); % Firts two days of output

runs=[3,4,5];
%runs=[3];

expt_names={'GLOBALWAVE','NZWAVE','NZWAVE-ST4','NZWAVE-ST6','NZWAVE-HR'};

file='Baring_Head';% 'SteepHead_SP_FullRecord_QC';
%file='wairewa_lake_forsyth';% 'SteepHead_SP_FullRecord_QC';
%file='SteepHead_SP_FullRecord_QC';
%file='Banks_Peninsula';

check_file=1;
proc_obs  =1;
save_csv  =0;
set_xlim  =1;
save_fig  =0;

ww3pre={'ww3p_interp','ww3g'};
ww3pre=ww3pre{2};

colors={'r','b','k','m','y'};

for i=1:length(runs)
  expts{i}=expt_names{runs(i)};
end

path_source=['/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/']; % GLOBALWAVE/'];%2018/01/05/00'];
path_santanarc='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/';
path_obs=[path_santanarc,'data/obs/'];

% OBS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_obs=[path_obs,file];
%    {'date'}    {'PeakPerid'}    {'PeakDirection'}    {'Directional spread'}    {'Tm01'}    {'Tm02'}    {'Hs'}    {'Qp'}

if proc_obs==1
  if strcmp(file,'SteepHead_SP_FullRecord_QC') || strcmp(file,'Banks_Peninsula')
    display(['Processing: ',file_obs,'.txt']);
    ob=importdata([file_obs,'.txt']);
    ob.textdata(1,:)
    for i=2:length(ob.textdata)
      time_obs(i-1)=datenum(ob.textdata{i,1},'YYYY-mm-ddTHH:MM:SS');
    end
    %time_obs=time_obs-.5;
    obs=ob.data;
    lat_obs=-43.7567 % from ECAN website -43+(45/60);
    lon_obs=173.3358 % 173+(20/60);
  elseif strcmp(file,'Baring_Head') 
    display(['Processing: ',file_obs,'.csv']);
    ob=importdata([file_obs,'.csv']);
    ob.textdata(1,:)
    for i=2:length(ob.textdata)
      time_obs(i-1)=datenum(ob.textdata{i,1},'dd/mm/YYYY HH:MM:SS');
    end
    time_obs=time_obs-.5;
    obs(:,6)=ob.data(:,1);
    obs(obs==0)=nan;
    %lat_obs=-41.434334; lon_obs=174.853727; 
    lat_obs=-41.416667; lon_obs=174.866667;

  elseif strcmp(file,'wairewa_lake_forsyth') 
    lat_obs=-43.84110; lon_obs=172.71835; % NZWAVE-HR 
    lat_obs=-43.84110; lon_obs=172.71835; % 
    time_obs=nan; obs(1,1:6)=nan;
  end
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
  
legh={'Obs'}; legt={'Obs'};

ke=0;
for expt=expts
  ke=ke+1;

  expt=expt{1};

  gnames={'globalwave+globalum','nzwave+nzlam','nzwave_hr+nzcsm','tongawave+globalum'};

  if strncmp(expt,'GLOBAL',6)
    ig=1;
    ltime=25;
  elseif strncmp(expt,'NZWAVE-HR',9)
    ig=3;
    ltime=49;
  elseif strncmp(expt,'NZWAVE',6)
    ig=2;
    ltime=25;
  end

  path_expt=[path_source,expt,'/']; % GLOBALWAVE/'];%2018/01/05/00'];
  
  path_dm=[path_expt,'matlab/'];
  system(['mkdir -p ',path_dm]);
  %filename=[path_dm,file,'_',datestr(time_lima(1),'YYYYmmDDHH'),'_',datestr(time_lima(end),'YYYYmmDDHH'),'.mat'];
  hs_mod=[];tp_mod=[]; pd_mod=[];time_mod=[];
  checkin=0;  
  % Loop in model time
  for t=time_lima
    filename=[path_dm,file,'_',datestr(t,'YYYYmmDDHH'),'.mat'];

    if check_file==1 && exist(filename)==2
  
      display(['Loading: ',filename])
      load(filename)%,'time_mod','model')
  
    else
      ptime=datestr(t,'YYYY/mm/DD/HH/');
      ftime=datestr(t,'YYYYmmDDHH');
    
      pname=[path_expt,ptime];
    
      fname=[pname,ww3pre,'_',ftime,'-utc_',gnames{ig},'.nc'];
      %fname=[pname,'ww3g_',ftime,'-utc_',gnames{1},'.nc'];
      display(['Loading: ',fname]);

      if checkin==0; % t==time_lima(1)
        checkin=1; 
        lon_mod=double(ncread(fname,'lon'));
        lat_mod=double(ncread(fname,'lat'));
        [dif ilon]=nanmin(abs(lon_mod-lon_obs));
        [dif ilat]=nanmin(abs(lat_mod-lat_obs));

        if strcmp(file,'SteepHead_SP_FullRecord_QC')
          %ilon=ilon+3;
        elseif strcmp(file,'wairewa_lake_forsyth') & ig==2 % NZWAVE
          ilat=ilat-2;
        end

        dis=sqrt((lon_mod(ilon)-lon_obs).^2 + (lat_mod(ilat)-lat_obs).^2)*110;
        display(['Distance between obs and grid point is: ',num2str(dis),' km'])
      end
      
      time=squeeze(double(ncread(fname,'time')))./24+t; time=time(1:ltime-1);
      if strncmp(ww3pre,'ww3p',4)
        hs=squeeze(double(ncread(fname,'hsig',          [ilon 1],[1 ltime]))); hs=hs(1:end-1)';
        tp=squeeze(double(ncread(fname,'tpeak',         [ilon 1],[1 ltime]))); tp=tp(1:end-1)';
        pd=squeeze(double(ncread(fname,'peak_direction',[ilon 1],[1 ltime]))); pd=pd(1:end-1)';
      elseif strncmp(ww3pre,'ww3g',4)
        hs=squeeze(double(ncread(fname,'hsig',          [ilon ilat 1],[1 1 ltime]))); hs=hs(1:end-1);
        tp=squeeze(double(ncread(fname,'tpeak',         [ilon ilat 1],[1 1 ltime]))); tp=tp(1:end-1);
        pd=squeeze(double(ncread(fname,'peak_direction',[ilon ilat 1],[1 1 ltime]))); pd=pd(1:end-1);
      end

      display(['Saving: ',filename])
      save(filename,'time','hs','tp','pd')

    end

      time_mod=[time_mod;time];
      hs_mod=[hs_mod;hs];
      tp_mod=[tp_mod;tp];
      pd_mod=[pd_mod;pd];
    
    clear model
 
    model.data(:,6)=hs_mod;
    model.data(:,1)=tp_mod;
    model.data(:,2)=pd_mod;
    
    %path_dm=[path_expt,'matlab/'];
    %system(['mkdir -p ',path_dm]);
    %display(['Saving: ',filename])
    %save(filename,'time_mod','model')
  
  end

  if save_csv==1

    fda=[filename(1:end-4),'.csv'];
    display(['Saving: ',fda]);

    if exist(fda)==2
       system(['rm -rf ',fda]);
    end

    model.data(:,3)=nan;
    model.data(:,4)=nan;
    model.data(:,5)=nan;
    model.data(:,7)=nan;
    csvdata=model.data;
    csvdata=num2cell(csvdata);
    for i=1:length(time_mod); csvtime{i,1}=datestr(time_mod(i),'yyyy-mm-ddTHH:MM:SS'); end
    csvhead={'date','PeakPeriod','PeakDirection','Directional spread','Tm01','Tm02','Hs','Qp'};

    % concatenating str matrices
    csvall=[csvtime,csvdata];
    %csvall=[csvhead;csvall];
    %csvwrite(fda,csvall);
    T = cell2table(csvall,'VariableNames',csvhead);
    writetable(T,fda)
    fclose('all');
  end


  dcol=[1,4,5]; %[6,1,2]; % Hs, Tp, dir
  dcol=[6,1,2]; % Hs, Tp, dir
  tnames={'Peak Period',    'Peak Direction',    'Directional spread',   'Tm01','Tm02','Significant Wave Height (Hs)','Qp'};
  for i=1:length(dcol)
  
    subplot(length(dcol),1,i);
    set(gca,'fontsize',12,'fontweight','bold')
    hold on
    if ke==1
      plot(time_obs,obs(:,dcol(i))','.','color',[0 .7 0],'linewidth',2)
    end
    plot(time_mod,model.data(:,dcol(i))','.','color',colors{ke},'linewidth',2)
  
    if set_xlim==1;
      xlim([time_lima(1)-2 time_lima(end)+2])
    end
  
    %if i==3
      datetick('x','dd/mmm/yyyy','keeplimits')
    %end
    title([replace(file,'_',' '),' ',tnames{dcol(i)}])
    grid('on')

    modeli=interp1(time_mod,model.data(:,dcol(i)),time_obs);

    if dcol(i)==6
      mod_rmse=nanrmse(obs(:,dcol(i)),modeli');
      mod_bias=nanmean(modeli')-nanmean(obs(:,dcol(i)));
      legh=[legh,{[expt,' rmse=',num2str(mod_rmse,'%.2f'),', bias=',num2str(mod_bias,'%.2f')]}];
    end

  end

end


subplot(length(dcol),1,find(dcol==6));
legend([legh],'location','best')


if save_fig==1
  export_fig(gcf,'steephead_period','-png','-r150');
end


