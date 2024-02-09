function [obs,time_obs,lon_obs,lat_obs]=proc_wave_obs(file,proc_obs) % ,plot_obs)


  path_santanarc='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/';
  path_obs=[path_santanarc,'data/obs/'];

  % OBS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  file_obs=[path_obs,file];
  % {'date'}    {'PeakPerid'}    {'PeakDirection'}    {'Directional spread'}    {'Tm01'}    {'Tm02'}    {'Hs'}    {'Qp'}
  
  if proc_obs==1 % & plot_obs==1
    if strcmp(file,'Banks_Peninsula') || strcmp(file,'SteepHead_SP_FullRecord_QC') 
      display(['Processing: ',file_obs,'.txt']);
      ob=importdata([file_obs,'.txt']);
      ob.textdata(1,:)
      for i=2:length(ob.textdata)
        time_obs(i-1)=datenum(ob.textdata{i,1},'YYYY-mm-ddTHH:MM:SS');
      end
      [dif iutc]=nanmin(abs(time_obs-datenum(2021,05,15)));
      time_obs(1:iutc)=time_obs(1:iutc)-.5;
      obs=ob.data;
      lat_obs=-43.7567 % from ECAN website -43+(45/60);
      lon_obs=173.3358 % 173+(20/60);

      % removing similar dates
      [time_obs,iorg,~]=unique(time_obs);
      for i=1:size(obs,2)
        obsn(:,i)=obs(iorg,i);
      end
      obs=obsn;

      % hourly moving mean
      % decomposing direction before moving average
      pd_obs=obs(:,2); hs_obs=obs(:,6);
      pd=mod(-90-pd_obs,360);
      u_obs=cosd(pd).*hs_obs; v_obs=sind(pd).*hs_obs;
      % hourly moving mean
      u_obs=movmean(u_obs,4); v_obs=movmean(v_obs,4);
      %hs_obs=sqrt(u_obs.^2+v_obs.^2); % squaring after moving mean doesnt work. it creates very small values sometimes
      dp_obs=270-atan2(v_obs,u_obs)*180/pi;
      dp_obs(dp_obs>360)=dp_obs(dp_obs>360)-360;

      %figure; hold on; plot(time_obs,obs(:,2),'.r');
      %plot(time_obs,dp_obs,'.b')
      %figure; hold on; plot(time_obs,obs(:,6),'.r');
      % hourly moving mean
      obs(:,6)=movmean(obs(:,6),4);
      obs(:,2)=dp_obs;
      %plot(time_obs,obs(:,6),'.g');
      %plot(time_obs,hs_obs,'.b')

      % 1 {'PeakPerid'}   2 {'PeakDirection'}   3 {'Directional spread'}   4 {'Tm01'}   5 {'Tm02'}   6 {'Hs'}   7 {'Qp'}  8 Tm
      % 9 wlv
      obs(:,8)=nan;
      obs(:,9)=nan;
  
    elseif strcmp(file,'Baring_Head') 
      %lat_obs=-41.434334; lon_obs=174.853727; 
      %lat_obs=-41.416667; lon_obs=174.866667;
      lat_obs=-41.4022; lon_obs=174.84669;
      % hs
      file_obsn=[file_obs,'_hs_20210101_20211231.csv'];
      display(['Processing: ',file_obsn]);
      ob=importdata([file_obsn]);
      ob.textdata(1,:)
      % time
      for i=2:length(ob.textdata)
        time_obs(i-1)=datenum(ob.textdata{i,1},'dd/mm/YYYY HH:MM:SS');
      end
      time_obs=time_obs-.5;
      obs(:,6)=ob.data(:,1);
      obs(obs==0)=nan;
      
      %figure; plot(time_obs,obs(:,6),'.k'); datetick('x','keeplimits')
      %return

      % tp
      file_obsn=[file_obs,'_tp_20210101_20211231.csv']; % Tsig
      display(['Processing: ',file_obsn]);
      ob=importdata([file_obsn]);
      ob.textdata(1,:)
      obs(:,1)=ob.data(:,1);
      inan=find(obs(:,1)==0); obs(inan,1)=nan;
      obs(:,4)=obs(:,1);  obs(:,5)=obs(:,1);
      % dp or pd
      file_obsn=[file_obs,'_dp_20210101_20211231.csv']; % peak direction 
      display(['Processing: ',file_obsn]);
      ob=importdata([file_obsn]);
      ob.textdata(1,:)
      obs(:,2)=ob.data(:,1);
      obs(inan,2)=nan;
 
      %bad dates
      time_ibad=[datenum(2021,09,2),datenum(2021,10,13)];
      time_fbad=[datenum(2021,10,1),datenum(2021,10,15)];
      for i=1:length(time_ibad)
        [dif ibad]=nanmin(abs(time_obs-time_ibad(i)));
        [dif fbad]=nanmin(abs(time_obs-time_fbad(i)));
        obs(ibad:fbad,:)=nan;
      end     

      % hourly moving mean
      % decomposing direction before moving average
      pd_obs=obs(:,2); hs_obs=obs(:,6);
      pd=mod(-90-pd_obs,360);
      u_obs=cosd(pd).*hs_obs; v_obs=sind(pd).*hs_obs;
      % hourly moving mean
      u_obs=movmean(u_obs,4); v_obs=movmean(v_obs,4);
      %hs_obs=sqrt(u_obs.^2+v_obs.^2); % squaring after moving mean doesnt work. it creates very small values sometimes
      dp_obs=270-atan2(v_obs,u_obs)*180/pi;
      dp_obs(dp_obs>360)=dp_obs(dp_obs>360)-360;

      %figure; hold on; plot(time_obs,obs(:,2),'.r');
      %plot(time_obs,dp_obs,'.b')
      %figure; hold on; plot(time_obs,obs(:,6),'.r');
      % hourly moving mean
      obs(:,6)=movmean(obs(:,6),4);
      obs(:,2)=dp_obs;
      %plot(time_obs,obs(:,6),'.g');
      %plot(time_obs,hs_obs,'.b')

      % 1 {'PeakPerid'}   2 {'PeakDirection'}   3 {'Directional spread'}   4 {'Tm01'}   5 {'Tm02'}   6 {'Hs'}   7 {'Qp'}  8 Tm
      % 9 wlv
      obs(:,8)=obs(:,1);
      obs(:,9)=nan;
 
    elseif strcmp(file,'wairewa_lake_forsyth') 
      lat_obs=-43.84110; lon_obs=172.71835; % NZWAVE-HR 
      time_obs=nan; obs(1,1:9)=nan;
 
    elseif strcmp(file,'Pelorus_Sound') 
      lat_obs=-40.955196; lon_obs=174.038144;
      time_obs=time_lima; obs(1:length(time_lima),1:9)=0;

    elseif strcmp(file,'Taumutu') 
      lat_obs=-43.866527; lon_obs=172.381380;
      time_obs=time_lima; obs(1:length(time_lima),1:9)=0;

    else %if strcmp(file,'Taumutu') 
      lat_obs=-43.866527; lon_obs=172.381380;
      time_obs=time_lima; obs(1:length(time_lima),1:9)=0;

    end

    save([file_obs,'.mat'],'time_obs','obs','lat_obs','lon_obs')
  
  elseif proc_obs==0 %& plot_obs==1
  
    display(['Loading: ',file_obs,'.mat']);
    load([file_obs,'.mat'])%,'time','obs')

  %else

  %  lon_obs=lon_obss(fobs);
  %  lat_obs=lat_obss(fobs);
  %  time_obs=time_lima; obs(1:length(time_lima),1:9)=0;

  end





