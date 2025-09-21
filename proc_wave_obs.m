function [obs,time_obs,lon_obs,lat_obs]=proc_wave_obs(file,proc_obs,time_lima,scase) % ,plot_obs)


  path_santanarc='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/';
  path_obs=[path_santanarc,'data/obs/'];

  % OBS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  file_obs=[path_obs,file];
  % {'date'}    {'PeakPerid'}    {'PeakDirection'}    {'Directional spread'}    {'Tm01'}    {'Tm02'}    {'Hs'}    {'Qp'}

  display([' ']);
  display([' ']);
  display(['Remember to QAQC wave data from stations:']);
  display(['Mokohinau']);
  display(['Pukehina']);
  display([' ']);
  display([' ']);
  
  if proc_obs==1 % & plot_obs==1

    display(['Working to save: ',file_obs,'.mat']);

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
      pd_obs=obs(:,2); 
      hs_obs=obs(:,6);
      pd=mod(-90-pd_obs,360);
      u_obs=cosd(pd).*hs_obs; v_obs=sind(pd).*hs_obs;
      % hourly moving mean
      u_obs=movmean(u_obs,4); v_obs=movmean(v_obs,4);
      %hs_obs=sqrt((u_obs.^2)+(v_obs.^2)); % squaring after moving mean doesnt work. it creates very small values sometimes
      dp_obs=270-atan2(v_obs,u_obs)*180/pi;
      dp_obs(dp_obs>360)=dp_obs(dp_obs>360)-360;

      %figure; hold on; plot(time_obs,obs(:,2),'.r');
      %plot(time_obs,dp_obs,'.b')
      %figure; hold on; plot(time_obs,obs(:,6),'.r');
      %plot(time_obs,hs_obs,'.b')

      % hourly moving mean
      obs(:,6)=movmean(obs(:,6),4);
      obs(:,6)=movmean(obs(:,6),4);
      obs(:,2)=dp_obs;
      %plot(time_obs,obs(:,6),'.g');
      %plot(time_obs,hs_obs,'.b')

      % 1 {'PeakPerid'}   2 {'PeakDirection'}   3 {'Directional spread'}   4 {'Tm01'}   5 {'Tm02'}   6 {'Hs'}   7 {'Qp'}  8 Tm
      % 9 wlv
      % hourly moving mean
      obs(:,5)=movmean(obs(:,5),4);
      obs(:,4)=obs(:,5);
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
      obs(:,1)=movmean(obs(:,1),4);

      inan=find(obs(:,1)==0); obs(inan,1)=nan;
      obs(:,4)=obs(:,1);  obs(:,5)=obs(:,1); % assigning Tp as Tm01 and Tm02
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
      obs(:,6)=movmean(obs(:,6),4);
      obs(:,2)=dp_obs;
      %plot(time_obs,obs(:,6),'.g');
      %plot(time_obs,hs_obs,'.b')

      % 1 {'PeakPerid'}   2 {'PeakDirection'}   3 {'Directional spread'}   4 {'Tm01'}   5 {'Tm02'}   6 {'Hs'}   7 {'Qp'}  8 Tm
      % 9 wlv
      obs(:,8)=obs(:,4);
      obs(:,9)=nan;

    elseif strcmp(file,'Baring_Head_2025') 
      %lat_obs=-41.434334; lon_obs=174.853727; 
      %lat_obs=-41.416667; lon_obs=174.866667;
      lat_obs=-41.4022; lon_obs=174.84669;
      % hs
      file_obsn=[file_obs,'_hs_north_20250701_20250718.csv'];
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
      
      obs(:,1:5)=nan; 
      obs(:,7:9)=nan; 

      %obs(:,6)=movmean(obs(:,6),4);
      %obs(:,6)=movmean(obs(:,6),4);


    elseif strcmp(file,'Mokohinau') 

      %clear 
      %close all

      file_obsn='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/data/obs/all_nz_bouys/buoy_1998043012_utc_2004061000_utc_arc_NZMokohinau_e17508s3588.nc';
      display(['Processing: ',file_obsn]);

      time_obs=ncread(file_obsn,'time')./24+datenum('1970-01-01 00:00'); 
      lat_obs=ncread(file_obsn,'lat'); 
      lon_obs=ncread(file_obsn,'lon'); 
      hs_obs=ncread(file_obsn,'hm0'); 
      tm02_obs=ncread(file_obsn,'tmean02'); 
      tp_obs=ncread(file_obsn,'tpeak'); 
      pd_obs=ncread(file_obsn,'peak_direction'); 
      ds_obs=ncread(file_obsn,'directional_spread'); 
      md_obs=ncread(file_obsn,'mean_direction'); 
      tm_obs=ncread(file_obsn,'tave'); 
      f_obs=ncread(file_obsn,'flag'); 
      
			%QAQC
      hs_obs=qaqc_extremes(hs_obs,0.23,8);

      % removing repeated times from wave data
      [time_obs,ia,ic]=unique(time_obs);
      hs_obs=hs_obs(ia); tm02_obs=tm02_obs(ia); tp_obs=tp_obs(ia); pd_obs=pd_obs(ia); ds_obs=ds_obs(ia); md_obs=md_obs(ia); tm_obs=tm_obs(ia); f_obs=f_obs(ia);

      time_int=[datenum(1998,4,30,12,30,0):1/24:time_obs(end)]'; 
      hs_int=interp1(time_obs,hs_obs,time_int);

			tm02_int=interp1(time_obs,tm02_obs,time_int);
			tp_int=interp1(time_obs,tp_obs,time_int);
			pd_int=interp1(time_obs,pd_obs,time_int);
			ds_int=interp1(time_obs,ds_obs,time_int);
			md_int=interp1(time_obs,md_obs,time_int);
			tm_int=interp1(time_obs,tm_obs,time_int);
			f_int=interp1(time_obs,f_obs,time_int);


      % inserting nans in long gaps longer than 3 hours
      hs_int=qaqc_time_gap(time_obs,time_int,hs_int,3/24);
			tm02_int=qaqc_time_gap(time_obs,time_int,tm02_int,3/24);
			tp_int=qaqc_time_gap(time_obs,time_int,tp_int,3/24);
			pd_int=qaqc_time_gap(time_obs,time_int,pd_int,3/24);
			ds_int=qaqc_time_gap(time_obs,time_int,ds_int,3/24);
			md_int=qaqc_time_gap(time_obs,time_int,md_int,3/24);
			tm_int=qaqc_time_gap(time_obs,time_int,tm_int,3/24);
			f_int=qaqc_time_gap(time_obs,time_int,f_int,3/24);

      %data=qaqc_remove_lin_interp(data,1.000002);

      %figure; hold on; 
      %plot(datetime(datevec(time_obs)),hs_obs,'.')
      %plot(datetime(datevec(time_int)),hs_int,'.k')

      time_obs=time_int; % time_obs is the time of the interpolated data
			hs_obs=hs_int;	tm02_obs=tm02_int;	tp_obs=tp_int;	pd_obs=pd_int; ds_obs=ds_int;	md_obs=md_int;	tm_obs=tm_int;	f_obs=f_int;


      %bad dates
      %time_ibad=[datenum(2021,09,2),datenum(2021,10,13)];
      %time_fbad=[datenum(2021,10,1),datenum(2021,10,15)];
      %for i=1:length(time_ibad)
        %[dif ibad]=nanmin(abs(time_obs-time_ibad(i)));
        %[dif fbad]=nanmin(abs(time_obs-time_fbad(i)));
        %obs(ibad:fbad,:)=nan;
      %end     

      % hourly moving mean
      % decomposing direction before moving average
      %pd=mod(-90-pd_obs,360);
      %u_obs=cosd(pd).*hs_obs; v_obs=sind(pd).*hs_obs;
      %% hourly moving mean
      %u_obs=movmean(u_obs,4); v_obs=movmean(v_obs,4);
      %%hs_obs=sqrt(u_obs.^2+v_obs.^2); % squaring after moving mean doesnt work. it creates very small values sometimes
      %dp_obs=270-atan2(v_obs,u_obs)*180/pi;
      %dp_obs(dp_obs>360)=dp_obs(dp_obs>360)-360;

      %figure; hold on; plot(time_obs,obs(:,2),'.r');
      %plot(time_obs,dp_obs,'.b')
      %figure; hold on; plot(time_obs,obs(:,6),'.r');
      % hourly moving mean
      %obs(:,6)=movmean(obs(:,6),4);
      %obs(:,6)=movmean(obs(:,6),4);

      %plot(time_obs,obs(:,6),'.g');
      %plot(time_obs,hs_obs,'.b')

      % 1 {'PeakPerid'}   2 {'PeakDirection'}   3 {'Directional spread'}   4 {'Tm01'}   5 {'Tm02'}   6 {'Hs'}   7 {'Qp'}  8 Tm
      % 9 wlv
      obs(:,8)=tm_obs';
      obs(:,9)=nan;
      obs(:,1)=tp_obs';
      obs(:,2)=pd_obs'; 
      obs(:,3)=ds_obs'; 
      obs(:,4)=nan; 
      obs(:,5)=tm02_obs'; 
      obs(:,6)=hs_obs';
      obs(:,7)=nan';


    elseif strcmp(file,'Pukehina') 

      %clear 
      %close all

      %buoy_2003092214_utc_2016101401_utc_ebop_NZPukehina_e17661s3769.nc

      file_obsn='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/data/obs/all_nz_bouys/buoy_2003092214_utc_2016101401_utc_ebop_NZPukehina_e17661s3769.nc';
      %file_obsn='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/data/obs/all_nz_bouys/buoy_2003092214_utc_2013013112_utc_ebop_NZPukehina_e17661s3769.nc';
      display(['Processing: ',file_obsn]);

      time_obs=ncread(file_obsn,'time')./24+datenum('1970-01-01 00:00'); 
      lat_obs=ncread(file_obsn,'lat'); 
      lon_obs=ncread(file_obsn,'lon'); 
      hs_obs=ncread(file_obsn,'hm0'); 
      tm02_obs=ncread(file_obsn,'tmean02'); 
      tp_obs=ncread(file_obsn,'tpeak'); 
      %pd_obs=ncread(file_obsn,'peak_direction'); 
      ds_obs=ncread(file_obsn,'directional_spread'); 
      md_obs=ncread(file_obsn,'mean_direction'); 
      tm_obs=ncread(file_obsn,'tave'); 
      %f_obs=ncread(file_obsn,'flag'); 
      
			%QAQC
      %hs_obs=qaqc_extremes(hs_obs,0.23,8);

      % removing repeated times from wave data
      [time_obs,ia,ic]=unique(time_obs);
      hs_obs=hs_obs(ia); tm02_obs=tm02_obs(ia); tp_obs=tp_obs(ia); 
			%pd_obs=pd_obs(ia); 
			ds_obs=ds_obs(ia); md_obs=md_obs(ia); tm_obs=tm_obs(ia); 

      time_int=[time_obs(1):1/24:time_obs(end)]'; 
      hs_int=interp1(time_obs,hs_obs,time_int);

			tm02_int=interp1(time_obs,tm02_obs,time_int);
			tp_int=interp1(time_obs,tp_obs,time_int);
			%pd_int=interp1(time_obs,pd_obs,time_int);
			ds_int=interp1(time_obs,ds_obs,time_int);
			md_int=interp1(time_obs,md_obs,time_int);
			tm_int=interp1(time_obs,tm_obs,time_int);
			%f_int=interp1(time_obs,f_obs,time_int);


      % inserting nans in long gaps longer than 3 hours
      hs_int=qaqc_time_gap(time_obs,time_int,hs_int,3/24);
			tm02_int=qaqc_time_gap(time_obs,time_int,tm02_int,3/24);
			tp_int=qaqc_time_gap(time_obs,time_int,tp_int,3/24);
			%pd_int=qaqc_time_gap(time_obs,time_int,pd_int,3/24);
			ds_int=qaqc_time_gap(time_obs,time_int,ds_int,3/24);
			md_int=qaqc_time_gap(time_obs,time_int,md_int,3/24);
			tm_int=qaqc_time_gap(time_obs,time_int,tm_int,3/24);
			%f_int=qaqc_time_gap(time_obs,time_int,f_int,3/24);

      %data=qaqc_remove_lin_interp(data,1.000002);

      %figure; hold on; 
      %plot(datetime(datevec(time_obs)),hs_obs,'.')
      %plot(datetime(datevec(time_int)),hs_int,'.k')

      %figure; hold on; 
      %plot(datetime(datevec(time_obs)),md_obs,'.')
      %plot(datetime(datevec(time_int)),md_int,'.k')

      time_obs=time_int; % time_obs is the time of the interpolated data
			hs_obs=hs_int;	tm02_obs=tm02_int;	tp_obs=tp_int;	
			%pd_obs=pd_int; 
			ds_obs=ds_int;	md_obs=md_int;	tm_obs=tm_int;	


      %bad dates
      %time_ibad=[datenum(2021,09,2),datenum(2021,10,13)];
      %time_fbad=[datenum(2021,10,1),datenum(2021,10,15)];
      %for i=1:length(time_ibad)
        %[dif ibad]=nanmin(abs(time_obs-time_ibad(i)));
        %[dif fbad]=nanmin(abs(time_obs-time_fbad(i)));
        %obs(ibad:fbad,:)=nan;
      %end     

      % hourly moving mean
      % decomposing direction before moving average
      %pd=mod(-90-pd_obs,360);
      %u_obs=cosd(pd).*hs_obs; v_obs=sind(pd).*hs_obs;
      %% hourly moving mean
      %u_obs=movmean(u_obs,4); v_obs=movmean(v_obs,4);
      %%hs_obs=sqrt(u_obs.^2+v_obs.^2); % squaring after moving mean doesnt work. it creates very small values sometimes
      %dp_obs=270-atan2(v_obs,u_obs)*180/pi;
      %dp_obs(dp_obs>360)=dp_obs(dp_obs>360)-360;

      %figure; hold on; plot(time_obs,obs(:,2),'.r');
      %plot(time_obs,dp_obs,'.b')
      %figure; hold on; plot(time_obs,obs(:,6),'.r');
      % hourly moving mean
      %obs(:,6)=movmean(obs(:,6),4);
      %obs(:,6)=movmean(obs(:,6),4);

      %plot(time_obs,obs(:,6),'.g');
      %plot(time_obs,hs_obs,'.b')

      % 1 {'PeakPerid'}   2 {'PeakDirection'}   3 {'Directional spread'}   4 {'Tm01'}   5 {'Tm02'}   6 {'Hs'}   7 {'Qp'}  8 Tm
      % 9 wlv
      obs(:,8)=tm_obs';
      obs(:,9)=nan;
      obs(:,1)=tp_obs';
      obs(:,2)=nan;% pd_obs'; 
      obs(:,3)=ds_obs'; 
      obs(:,4)=nan; 
      obs(:,5)=tm02_obs'; 
      obs(:,6)=hs_obs';
      obs(:,7)=nan';


    elseif strcmp(file,'Taharoa') 

      %clear; close all

      file_obsn='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/data/obs/all_nz_bouys/buoy_1996030704_utc_2000123112_utc_nzsteel_NZTaharoa_e17467s3816.nc';
      display(['Processing: ',file_obsn]);

      time_obs=ncread(file_obsn,'time')./24+datenum('1970-01-01 00:00') ; 
      lat_obs=ncread(file_obsn,'lat'); 
      lon_obs=ncread(file_obsn,'lon') -(4/100); 
      hs_obs=ncread(file_obsn,'hsig'); 
      %tm02_obs=ncread(file_obsn,'tmean02'); 
      tp_obs=ncread(file_obsn,'tpeak'); 
      %pd_obs=ncread(file_obsn,'peak_direction'); 
      %ds_obs=ncread(file_obsn,'directional_spread'); 
      %md_obs=ncread(file_obsn,'mean_direction'); 
      tm_obs=ncread(file_obsn,'tave'); 
      %f_obs=ncread(file_obsn,'flag'); 
      
			%QAQC
      hs_obs=qaqc_extremes(hs_obs,0.23,10);

      % removing repeated times from wave data
      [time_obs,ia,ic]=unique(time_obs);
      hs_obs=hs_obs(ia); 
      %tm02_obs=tm02_obs(ia); 
      tp_obs=tp_obs(ia); 
			%pd_obs=pd_obs(ia); 
			%ds_obs=ds_obs(ia); md_obs=md_obs(ia); 
      tm_obs=tm_obs(ia); 

      time_int=[time_obs(1):1/24:time_obs(end)]'; 
      hs_int=interp1(time_obs,hs_obs,time_int);
			%tm02_int=interp1(time_obs,tm02_obs,time_int);
			tp_int=interp1(time_obs,tp_obs,time_int);
			%pd_int=interp1(time_obs,pd_obs,time_int);
			%ds_int=interp1(time_obs,ds_obs,time_int);
			%md_int=interp1(time_obs,md_obs,time_int);
			tm_int=interp1(time_obs,tm_obs,time_int);
			%f_int=interp1(time_obs,f_obs,time_int);


      % inserting nans in long gaps longer than 3 hours
      hs_int=qaqc_time_gap(time_obs,time_int,hs_int,3/24);
			%tm02_int=qaqc_time_gap(time_obs,time_int,tm02_int,3/24);
			tp_int=qaqc_time_gap(time_obs,time_int,tp_int,3/24);
			%pd_int=qaqc_time_gap(time_obs,time_int,pd_int,3/24);
			%ds_int=qaqc_time_gap(time_obs,time_int,ds_int,3/24);
			%md_int=qaqc_time_gap(time_obs,time_int,md_int,3/24);
			tm_int=qaqc_time_gap(time_obs,time_int,tm_int,3/24);
			%f_int=qaqc_time_gap(time_obs,time_int,f_int,3/24);

      data=hs_int;
      %data=qaqc_remove_lin_interp(data,0.000002);
      hs_int=data;

      %figure; hold on; 
      %plot(datetime(datevec(time_obs)),hs_obs,'.')
      %plot(datetime(datevec(time_int)),hs_int,'.k')

      %figure; hold on; 
      %plot(datetime(datevec(time_obs)),md_obs,'.')
      %plot(datetime(datevec(time_int)),md_int,'.k')

      time_obs=time_int; % time_obs is the time of the interpolated data
			hs_obs=hs_int;	%tm02_obs=tm02_int;	
      tp_obs=tp_int;	
			%pd_obs=pd_int; 
			%ds_obs=ds_int;	md_obs=md_int;	
			tm_obs=tm_int;	


      %bad dates
      %time_ibad=[datenum(2021,09,2),datenum(2021,10,13)];
      %time_fbad=[datenum(2021,10,1),datenum(2021,10,15)];
      %for i=1:length(time_ibad)
        %[dif ibad]=nanmin(abs(time_obs-time_ibad(i)));
        %[dif fbad]=nanmin(abs(time_obs-time_fbad(i)));
        %obs(ibad:fbad,:)=nan;
      %end     

      % hourly moving mean
      % decomposing direction before moving average
      %pd=mod(-90-pd_obs,360);
      %u_obs=cosd(pd).*hs_obs; v_obs=sind(pd).*hs_obs;
      %% hourly moving mean
      %u_obs=movmean(u_obs,4); v_obs=movmean(v_obs,4);
      %%hs_obs=sqrt(u_obs.^2+v_obs.^2); % squaring after moving mean doesnt work. it creates very small values sometimes
      %dp_obs=270-atan2(v_obs,u_obs)*180/pi;
      %dp_obs(dp_obs>360)=dp_obs(dp_obs>360)-360;

      %figure; hold on; plot(time_obs,obs(:,2),'.r');
      %plot(time_obs,dp_obs,'.b')
      %figure; hold on; plot(time_obs,obs(:,6),'.r');
      % hourly moving mean
      %obs(:,6)=movmean(obs(:,6),4);
      %obs(:,6)=movmean(obs(:,6),4);

      %plot(time_obs,obs(:,6),'.g');
      %plot(time_obs,hs_obs,'.b')

      % 1 {'PeakPerid'}   2 {'PeakDirection'}   3 {'Directional spread'}   4 {'Tm01'}   5 {'Tm02'}   6 {'Hs'}   7 {'Qp'}  8 Tm
      % 9 wlv
      obs(:,8)=tm_obs';
      obs(:,9)=nan;
      obs(:,1)=tp_obs';
      obs(:,2)=nan;% pd_obs'; 
      obs(:,3)=nan;%ds_obs'; 
      obs(:,4)=nan; 
      obs(:,5)=nan;%tm02_obs'; 
      obs(:,6)=hs_obs';
      obs(:,7)=nan';

 
    elseif strcmp(file,'wairewa_lake_forsyth') 
      lat_obs=-43.84110; lon_obs=172.71835; % NZWAVE-HR 
      time_obs=nan; obs(1,1:9)=nan;
 
    elseif strcmp(file,'Pelorus_Sound') 
      lat_obs=-40.955196; lon_obs=174.038144;
      time_obs=time_lima; obs(1:length(time_lima),1:9)=0;

    elseif strcmp(file,'Taumutu') 
      lat_obs=-43.866527; lon_obs=172.381380;
      time_obs=time_lima; obs(1:length(time_lima),1:9)=0;

    elseif strcmp(file,'Mangawhai') 
      lat_obs=-36.114621; lon_obs=174.651937;
      time_obs=nan; obs(1,1:9)=nan;

    elseif strcmp(file,'BoP_leatherback_turtles') 
      lat_obs=-37.26; lon_obs=177.59;
      time_obs=nan; obs(1,1:9)=nan;

    elseif strcmp(file,'Hjort_Trench')
      lat_obs=-58.497015; lon_obs=157.814017;
      time_obs=nan; obs(1,1:9)=nan;

    elseif strcmp(file,'Puysegur_Trench')
      lat_obs=-47.767419; lon_obs=164.576396;
      time_obs=nan; obs(1,1:9)=nan;

    elseif strcmp(file,'Port_Waikato')
      lat_obs=-37.403755; lon_obs=174.595955;
      time_obs=nan; obs(1,1:9)=nan;

    elseif strcmp(file,'Blenheim')
      lat_obs=-41.4360; lon_obs=174.23765;
      time_obs=nan; obs(1,1:9)=nan;

    %elseif strcmp(file,'Tairua')
    %  lat_obs=-36.885410; lon_obs=176.069924; 
    %  time_obs=nan; obs(1,1:9)=nan;

    else
 
      disp('NO OBSERVATIONS AVAILABLE. GETTING LON and LAT FROM read_wave_stations.m')
      %lat_obs=-43.866527; lon_obs=172.381380;
      [tstation,lat_obss,lon_obss]=read_wave_stations(scase); % used in data request
      for i=1:length(tstation);
        if strcmp(tstation{i},file); break; end
      end
      lon_obs=lon_obss(i); lat_obs=lat_obss(i);
      time_obs=time_lima; obs(1:length(time_lima),1:9)=0;

    end


    % 
    save([file_obs,'.mat'],'time_obs','obs','lat_obs','lon_obs')



  
  elseif proc_obs==0 % & plot_obs==1
     
    display(['Loading: ',file_obs,'.mat']);
    load([file_obs,'.mat'])%,'time','obs')
    
  %else
    
  %  lon_obs=lon_obss(fobs);
  %  lat_obs=lat_obss(fobs);
  %  time_obs=time_lima; obs(1:length(time_lima),1:9)=0;

  end





