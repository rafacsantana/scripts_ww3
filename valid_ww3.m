clear
close all
warning off

% Daily time
time_lima=datenum(2014,7,9,0,0,0):1:datenum(2014,7,14,0,0,0); % initial date
%time_lima=datenum(2014,8,5,0,0,0):1:datenum(2014,8,17,0,0,0); % 9-m wave event
%time_lima=datenum(2018,1,4,0,0,0):1:datenum(2018,1,10,0,0,0); % crossing between hindcast and ecoconnect
time_lima=datenum(2021,1,9,0,0,0):1:datenum(2021,1,16,0,0,0); % Firts two days of output
%time_lima=datenum(2021,5,20,0,0,0):1:datenum(2021,6,1,0,0,0); % Baring Head tide and swell
%time_lima=datenum(2021,10,19,0,0,0):1:datenum(2021,10,24,0,0,0); % Banks Peninsula head tide and swell
%time_lima=datenum(2021,1,1,0,0,0):1:datenum(2021,12,2,0,0,0); % Firts two days of output
%time_lima=datenum(2021,1,1,0,0,0):1:datenum(2021,12,31,0,0,0); % Firts two days of output

runs=[1:6];
runs=6;
expt_names={'GLOBALWAVE','NZWAVE','NZWAVE-ST4','NZWAVE-ST6','NZWAVE-HR-NOTIDES','NZWAVE-HR'};

stations=[2];
files={'Banks_Peninsula','Baring_Head'};
%file='wairewa_lake_forsyth';% 'SteepHead_SP_FullRecord_QC';

plot_series=0;
plot_map   =1;

ck_mod_mat =1;
proc_obs   =1;
save_csv   =0;
set_xlim   =1;
save_fig   =1;

% 1 {'PeakPerid'}   2 {'PeakDirection'}   3 {'Directional spread'}   4 {'Tm01'}   5 {'Tm02'}   6 {'Hs'}   7 {'Qp'}  8 Tm
% 9 wlv
dcol=[1,4,5,9]; %[6,1,2]; % Tp, Tm01, Tm02
dcol=[4,5,8,9]; %[6,1,2]; % Tp, Tm01, Tm02
%dcol=[6,1,2]; % Hs, Tp, dir
%dcol=[2,9]; % Hs, Tp, dir

pcname ='wlv'; % m_pcolor:
vecname='curr'; % m_vec: hs, curr

colors={'r','b','k','m','y','c'};

for i=1:length(runs)
  expts{i}=expt_names{runs(i)};
end

path_source=['/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/']; % GLOBALWAVE/'];%2018/01/05/00'];
path_santanarc='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/';
path_obs=[path_santanarc,'data/obs/'];
path_fig='/scale_wlg_nobackup/filesets/nobackup/niwa03150/santanarc/figures/';

% stations or grids - stations are bad in Baring_Head
ww3pre={'ww3p_interp','ww3g'};
ww3pre=ww3pre{2};


for fobs=stations

  file=files{fobs};
  
  % OBS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  file_obs=[path_obs,file];
  %    {'date'}    {'PeakPerid'}    {'PeakDirection'}    {'Directional spread'}    {'Tm01'}    {'Tm02'}    {'Hs'}    {'Qp'}
  
  if proc_obs==1
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
      lat_obs=-41.40; lon_obs=174.85;
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
      lat_obs=-43.84110; lon_obs=172.71835; % 
      time_obs=nan; obs(1,1:9)=nan;
  
    end
    save([file_obs,'.mat'],'time_obs','obs','lat_obs','lon_obs')
  
  else
  
    display(['Loading: ',file_obs,'.mat']);
    load([file_obs,'.mat'])%,'time','obs')
  
  end

  if plot_series
  
    scrsz=[1 1 1366 768];
    %scrsz=get(0,'screensize');
    
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
      tm01_mod=[];tm_mod=[]; pd_mod=[];time_mod=[];
      tm02_mod=[];ds_mod=[]; wlv_mod=[];time_mod=[];
    
      checkin=0;  
      % Loop in model time
      for t=time_lima
        filename=[path_dm,file,'_',datestr(t,'YYYYmmDDHH'),'.mat'];
    
        if ck_mod_mat==1 && exist(filename)==2
      
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
            tm01=squeeze(double(ncread(fname,'tmean01',     [ilon 1],[1 ltime]))); tm01=tm01(1:end-1)';
            tm=squeeze(double(ncread(fname,'tmean',         [ilon 1],[1 ltime]))); tm=tm(1:end-1)';
            pd=squeeze(double(ncread(fname,'peak_direction',[ilon 1],[1 ltime]))); pd=pd(1:end-1)';
          elseif strncmp(ww3pre,'ww3g',4)
            hs=squeeze(double(ncread(fname,'hsig',          [ilon ilat 1],[1 1 ltime]))); hs=hs(1:end-1);
            tp=squeeze(double(ncread(fname,'tpeak',         [ilon ilat 1],[1 1 ltime]))); tp=tp(1:end-1);
            tm01=squeeze(double(ncread(fname,'tmean01',     [ilon ilat 1],[1 1 ltime]))); tm01=tm01(1:end-1);
            tm02=squeeze(double(ncread(fname,'tmean02',     [ilon ilat 1],[1 1 ltime]))); tm02=tm02(1:end-1);
            tm=squeeze(double(ncread(fname,'tmean',         [ilon ilat 1],[1 1 ltime]))); tm=tm(1:end-1);
            pd=squeeze(double(ncread(fname,'peak_direction',[ilon ilat 1],[1 1 ltime]))); pd=pd(1:end-1);
            ds=squeeze(double(ncread(fname,'directional_spread',[ilon ilat 1],[1 1 ltime]))); ds=ds(1:end-1);
            wlv=squeeze(double(ncread(fname,'wlv',[ilon ilat 1],[1 1 ltime]))); wlv=wlv(1:end-1);
          end
    
          display(['Saving: ',filename])
          save(filename,'time','tp','pd','ds','tm01','tm02','hs','tm','wlv')
  
        end % if ck_mod_mat==1 && exist(filename)==2
    
        time_mod=[time_mod;time];
        hs_mod=[hs_mod;hs];
        tp_mod=[tp_mod;tp];
        tm01_mod=[tm01_mod;tm01];
        tm02_mod=[tm02_mod;tm02];
        tm_mod=[tm_mod;tm];
        pd_mod=[pd_mod;pd];
        ds_mod=[ds_mod;ds];
        wlv_mod=[wlv_mod;wlv];
        
        clear model
     
        % 1 {'PeakPerid'}   2 {'PeakDirection'}   3 {'Directional spread'}   4 {'Tm01'}   5 {'Tm02'}   6 {'Hs'}   7 {'Qp'}  8 Tm
        model.data(:,9)=wlv_mod;
        model.data(:,8)=tm_mod;
        model.data(:,1)=tp_mod;
        model.data(:,2)=pd_mod;
        model.data(:,3)=ds_mod;
        model.data(:,4)=tm01_mod;
        model.data(:,5)=tm02_mod;
        model.data(:,6)=hs_mod;
        
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
    
    
      tnames={'Peak Period',    'Peak Direction',    'Directional spread',   'Tm01','Tm02','Significant Wave Height (Hs)','Qp','Tm','Water Level'};
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
          mod_corr=nancorr(obs(:,dcol(i)),modeli');
          mod_bias=nanmean(modeli')-nanmean(obs(:,dcol(i)));
          legh=[legh,{[expt,' rmse=',num2str(mod_rmse,'%.2f'),' corr=',num2str(mod_corr,'%.2f'),', bias=',num2str(mod_bias,'%.2f')]}];
        end
        
      end
    
    end
    
    if ~isempty(find(dcol==6))
      subplot(length(dcol),1,find(dcol==6));
      legend([legh],'location','best')
    end
    
    if save_fig==1
      %export_fig(gcf,'steephead_period','-png','-r150');
    end

  end % if plot_series

  if plot_map

    scrsz=[1 1 1366 768];
    scrsz=[1 1 1920 1080];
    %scrsz=get(0,'screensize');
    
    figure('position',scrsz,'color',[1 1 1],'visible','on')
    hold on
    set(gca,'fontsize',12,'fontweight','bold')
      
    legh={'Obs'}; legt={'Obs'};

    checkin=0;  
    % Loop in desired time
    for t=time_lima
    
      ptime=datestr(t,'YYYY/mm/DD/HH/');
      ftime=datestr(t,'YYYYmmDDHH');
    
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
    
        pname=[path_expt,ptime];
        
        fname=[pname,ww3pre,'_',ftime,'-utc_',gnames{ig},'.nc'];
        %fname=[pname,'ww3g_',ftime,'-utc_',gnames{1},'.nc'];
        display(['Loading: ',fname]);
    
        if checkin==0; % t==time_lima(1)
          checkin=1; 
          depth_mod=double(ncread(fname,'depth',[1 1 1],[Inf Inf 1]));
          wlv_mod=double(ncread(fname,'wlv',[1 1 1],[Inf Inf 1]));
          depth_mod=depth_mod-wlv_mod;

          lon_mod=double(ncread(fname,'lon'));
          lat_mod=double(ncread(fname,'lat'));
          [dif ilon]=nanmin(abs(lon_mod-lon_obs));
          [dif ilat]=nanmin(abs(lat_mod-lat_obs));

          [lon_modm,lat_modm]=meshgrid(lon_mod,lat_mod);
 
          if strcmp(file,'SteepHead_SP_FullRecord_QC')
            %ilon=ilon+3;
          elseif strcmp(file,'wairewa_lake_forsyth') & ig==2 % NZWAVE
            ilat=ilat-2;
          end
    
          dis=sqrt((lon_mod(ilon)-lon_obs).^2 + (lat_mod(ilat)-lat_obs).^2)*110;
          display(['Distance between obs and grid point is: ',num2str(dis),' km'])
        end
        
        display(['Reading nc variables']); tic
        time=squeeze(double(ncread(fname,'time')))./24+t; time=time(1:ltime-1);
        hs=squeeze(double(ncread(fname,'hsig',          [1 1 1],[Inf Inf ltime-1]))); 
        tp=squeeze(double(ncread(fname,'tpeak',         [1 1 1],[Inf Inf ltime-1]))); 
        %tm01=squeeze(double(ncread(fname,'tmean01',     [1 1 1],[Inf Inf ltime-1]))); 
        %tm=squeeze(double(ncread(fname,'tmean',         [1 1 1],[Inf Inf ltime-1]))); 
        pd=squeeze(double(ncread(fname,'peak_direction',[1 1 1],[Inf Inf ltime-1]))); 
        if strcmp(expt,'NZWAVE-HR')
          wlv=squeeze(double(ncread(fname,'wlv',          [1 1 1],[Inf Inf ltime-1]))); 
          model(ke).wlv=wlv;
          ucur=squeeze(double(ncread(fname,'ucur',        [1 1 1],[Inf Inf ltime-1]))); 
          model(ke).ucur=ucur;
          vcur=squeeze(double(ncread(fname,'vcur',        [1 1 1],[Inf Inf ltime-1]))); 
          model(ke).vcur=vcur;
        else
          model(ke).wlv=nan(size(hs));
          model(ke).ucur=nan(size(hs));
          model(ke).vcur=nan(size(hs));
        end
        display(['Finished reading nc variables']); toc

        model(ke).time=time;
        model(ke).hs=hs;
        model(ke).tp=tp;
        model(ke).pd=pd;
        %model(ke).=;

      end % for expt=expts

      time=model(1).time;
      tnames={'Peak Period',    'Peak Direction',    'Directional spread',   'Tm01','Tm02 or Tm','Significant Wave Height (Hs)','Qp'};

      % converting to u and v before interpolating
      pd_obs=obs(:,2); hs_obs=obs(:,6);
      pd=mod(-90-pd_obs,360); u_obs=cosd(pd).*hs_obs; v_obs=sind(pd).*hs_obs;
      u_obsl=interp1(time_obs,u_obs,time); v_obsl=interp1(time_obs,v_obs,time);
      u_obsc=interp1(time_obs,u_obs,time,'cubic'); v_obsc=interp1(time_obs,v_obs,time,'cubic');
      dp_obs=270-atan2(v_obs,u_obs)*180/pi;
      dp_obs(dp_obs>360)=dp_obs(dp_obs>360)-360;

      pd_obs=dp_obs;
      hs_obs=interp1(time_obs,obs(:,6),time);
      tp_obs=interp1(time_obs,obs(:,1),time);

        
      for i=1:length(time)
 
        ke=0;
        for expt=expts
          ke=ke+1;

          display(['Ploting map']); tic
          ax(ke)=subplot(1,length(expts),ke);
          hold on
          m_proj('lambert','long', [lon_obs-.3 lon_obs+.3],'lat',[lat_obs-.3 lat_obs+.3]);

          if strcmp(pcname,'wlv')
            data=model(ke).wlv(:,:,i);
          %elseif pcname='wlv'
          elseif strcmp(pcname,'curr')
            u_mod=model(ke).ucur(:,:,i);
            v_mod=model(ke).vcur(:,:,i);
            data=sqrt(u_mod.^2+v_mod.^2);
          end
          %colormap(ax(1),gwbreal)
          %m_pcolor(lon_ccmp,lat_ccmp,wsc_ccmp');
          %set(get(cs,'ylabel'),'string','N/m^3','fontsize',20,'fontweight','bold');
          %caxis([-2E-06 2E-06])
          cmap=cmocean('balance');
          colormap(ax(ke),cmap) % parula, avhrr
          m_pcolor(lon_mod,lat_mod,data');
          caxis([-2 2])
          %caxis([0 .13])
          cb=colorbar;%('southoutside');
          set(get(cb,'ylabel'),'string','Water level (m)','fontsize',12,'fontweight','bold');
          set(cb,'fontsize',12,'fontweight','bold');

          %%set(cs,'position','south');
          %shading interp
          %%caxis([nanmin(t_avhrr(:)) nanmax(t_avhrr(:))]); % 12 25

          %[cs,h]=m_contour(lon_bath,lat_bath,bath',[-1000 -1000],'color',[1 1 1],'linewidth',2);
          %clabel(cs,h,'color',[1 1 1],'fontsize',14,'fontweight','bold','LabelSpacing',2000)
          [cs,h]=m_contour(lon_mod,lat_mod,-depth_mod',[-2 -20 -50 -100 -200],'color',[1 1 1],'linewidth',2);
          clabel(cs,h,'color',[0 0 0],'fontsize',10,'fontweight','bold','LabelSpacing',2000)

          %[cs,h]=m_contour(lon_aviso,lat_aviso,z_aviso_avg,'color','k');
          %[cs,h]=m_contour(lon_aviso,lat_aviso,z_aviso_avg,[.76:.04:.88],'color','k');
          %clabel(cs,h,'color',[0 0 0],'fontsize',14,'fontweight','bold','LabelSpacing',2000) % each contour line represents 0.02 increase/decrease

          m_gshhs_f('patch',[.7 .7 .7],'edgecolor','k','linewidth',.5);
          %%m_grid('fontname','helvetica','fontsize',22,'fontweight','bold');
          %%m_grid('xtick',[round(loni):2:round(lonf)],'xticklabels',[round(loni):2:round(lonf)],'fontname','helvetica','fontsize',20,'fontweight','bold');
          m_grid('fontname','helvetica','fontsize',14,'fontweight','bold');
          title([replace(file,'_',' '),' Hs and water level at ',datestr(time(i),'HH:MM dd/mmm/yyyy')],'fontsize',14,'fontweight','bold')

          %set(gca,'fontsize',14,'fontname','helvetica','fontweight','bold')
          %%m_contour(plon,plat,bathy',[200 2000],'w','linewidth',1);
          %%caxis([0 .5])

          if strcmp(vecname,'wlv')
          elseif strcmp(vecname,'hs')
            data=model(ke).hs(:,:,i);
            pd_mod=model(ke).pd(:,:,i);
            pd=mod(-90-pd_mod,360);
            u_mod=cosd(pd).*data; v_mod=sind(pd).*data;
          elseif strcmp(vecname,'curr')
            u_mod=model(ke).ucur(:,:,i);
            v_mod=model(ke).vcur(:,:,i);
            data=sqrt(u_mod.^2+v_mod.^2);
          end

          v_spa=1;
          if strcmp(vecname,'wlv')
            m_vec(10,lon_modm(1:v_spa:end,1:v_spa:end),lat_modm(1:v_spa:end,1:v_spa:end),u_mod(1:v_spa:end,1:v_spa:end)',v_mod(1:v_spa:end,1:v_spa:end)','k','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5)
            [hp ht]=m_vec(10,lon_obs+.1,lat_obs+.14,2,0,'k','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5,'key',[num2str(2.0),' m']);
            set(ht,'fontsize',14,'fontweight','bold')
          elseif strcmp(vecname,'curr')
            m_vec(5,lon_modm(1:v_spa:end,1:v_spa:end),lat_modm(1:v_spa:end,1:v_spa:end),u_mod(1:v_spa:end,1:v_spa:end)',v_mod(1:v_spa:end,1:v_spa:end)','k','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5)
            [hp ht]=m_vec(5,lon_obs+.1,lat_obs+.14,2,0,'k','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5,'key',[num2str(2.0),' m/s']);
            set(ht,'fontsize',14,'fontweight','bold')
          end

          pd_mod=model(ke).pd(:,:,i);
          m_text(lon_obs+.1,lat_obs+.2,['Mod dp: ',num2str(pd_mod(ilon,ilat),'%.1f'),'^o'],'fontsize',14,'fontweight','bold') 
          m_text(lon_obs+.1,lat_obs+.18,['Obs dp: ',num2str(pd_obs(i),'%.2f'),'^o'],'fontsize',14,'fontweight','bold') 

          % model single vector 
          data=model(ke).hs(:,:,i);
          pd_mod=model(ke).pd(:,:,i);
          pd=mod(-90-pd_mod,360);
          u_mod=cosd(pd).*data; v_mod=sind(pd).*data;
          m_vec(10,lon_mod(ilon),lat_mod(ilat),u_mod(ilon,ilat),v_mod(ilon,ilat),'b','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5)
          [hp ht]=m_vec(10,lon_obs+.1,lat_obs+.12,2,0,'b','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5,'key',[num2str(2.0),' m']);
          set(ht,'fontsize',14,'fontweight','bold')
          
          % obs single vector 
          %pd=mod(-90-pd_obs,360);
          %u_obs=cosd(pd).*hs_obs; v_obs=sind(pd).*hs_obs;
          m_vec(10,lon_obs,lat_obs,u_obsc(i),v_obsc(i),'m','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5)
          %m_vec(10,lon_obs,lat_obs,u_obsl(i),v_obsl(i),'c','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5)

          display(['Finished ploting map']); toc

          %if dcol(i)==6
          %  mod_rmse=nanrmse(obs(:,dcol(i)),modeli');
          %  mod_bias=nanmean(modeli')-nanmean(obs(:,dcol(i)));
          %  legh=[legh,{[expt,' rmse=',num2str(mod_rmse,'%.2f'),', bias=',num2str(mod_bias,'%.2f')]}];
          %end

        end % for expt=expts

        path_dm=strcat(path_fig,expt{1},'/');
        system(['mkdir -p ',path_dm]);
        figname=[path_dm,file,'_',pcname,'_',vecname,'_',datestr(time(i),'YYYY_mm_DD_HH_MM'),'.png'];
        display(['Plotting: ',figname]);
        if save_fig==1
          display(['Saving: ',figname]);
          export_fig(gcf,figname,'-png','-r150');
          %close
          clf('reset')
        else
          pause
          clf('reset')
        end
          
      end % for i=1:length(time)

    end % for t=time_lima
    
    close 

  end % if plot_map


end % fobs=fstations
