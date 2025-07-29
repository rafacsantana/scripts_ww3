clear
close all
%run('/scale_wlg_persistent/filesets/home/santanarc/scripts/niwa/matlab/startup.m')

% Daily time
time_lima=datenum(2023,8,9,0,0,0):1:datenum(2023,8,13,0,0,0); % MPI request on 28/08/2023
time_lima=datenum(2018,1,7,0,0,0):1:datenum(2023,9,17,0,0,0); % Graham Harrington ECAN request on 18/09/2023
time_lima=datenum(2014,7,9,0,0,0):1:datenum(2018,1,6,0,0,0); % Graham Harrington ECAN request on 02/10/2023
time_lima=datenum(2018,1,6,0,0,0):1:datenum(2023,8,13,0,0,0); % Graham Harrington ECAN request on 18/09/2023
time_lima=datenum(2019,1,1,0,0,0):1:datenum(2023,12,31,0,0,0); % Graham Harrington ECAN request on 18/09/2023

time_lima=datenum(2019,1,1,0,0,0):1:datenum(2019,12,31,0,0,0); % Graham Harrington ECAN request on 18/09/2023
%time_lima=datenum(2020,1,1,0,0,0):1:datenum(2020,12,31,0,0,0); % Graham Harrington ECAN request on 18/09/2023
%time_lima=datenum(2021,1,1,0,0,0):1:datenum(2021,12,31,0,0,0); % Graham Harrington ECAN request on 18/09/2023
time_lima=datenum(2022,1,1,0,0,0):1:datenum(2022,12,31,0,0,0); % Graham Harrington ECAN request on 18/09/2023
time_lima=datenum(2023,1,1,0,0,0):1:datenum(2023,12,31,0,0,0); % Graham Harrington ECAN request on 18/09/2023
time_lima=datenum(2024,1,1,0,0,0):1:datenum(2024,8,20,0,0,0); % Graham Harrington ECAN request on 18/09/2023

time_lima=datenum(2019,1,1,0,0,0):1:datenum(2024,8,20,0,0,0); % Graham Harrington ECAN request on 18/09/2023

%time_lima=datenum(2019,1,1,0,0,0):1:datenum(2019,1,1,0,0,0); % Graham Harrington ECAN request on 18/09/2023

%time_lima=datenum(2021,1,1,0,0,0):1:datenum(2021,1,4,0,0,0); % Graham Harrington ECAN request on 18/09/2023
%time_lima=datenum(2023,2,8,0,0,0):1:datenum(2023,2,16,0,0,0); % cyclone gabrielle

% switches
check_file=1; % check if there is a saved matlab file
ck_l_file =1; % check if the whole time series is available
load_atm  =1;
set_xlim  =0;
proc_obs  =0;

plot_fig  =0; % plot short time seires of wind and waves
plot_rose =0;
plot_hist =1;

save_fig  =0;
save_nc   =0;
save_csv  =0;
portrait  =0;

for dcol=[6]
if dcol==3
   rlevels=0:5:55;
elseif dcol==6
   rlevels=0:.5:10.5;
end


path_fig=['/scale_wlg_nobackup/filesets/nobackup/niwa03150/santanarc/figures/']; % GLOBALWAVE/'];%2018/01/05/00'];

path_efs='/niwa/archive/ecoconnect/EFS/';
path_efs='/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/';

prefixes={'tn_','ww3g_'};
gnames  ={'nzcsm','nzwave_hr+nzcsm'};
%gnames  ={'nzcsm','nzwave+nzlam'};

% always atm    wave  models
expts={'NZCSM','NZWAVE-HR'};

if portrait
  scrsz=[1 1 1366 768];
  scrsz=get(0,'screensize');
else
  scrsz=[2 42 958 953];
end
  
colors={'k','k'};

tstation={'Sinclair Head','Ohau Head','Mana Island','Makara'  ,'Karori Rock','Wairewa Lake Forsyth','Taumutu' ,'Gulf Harbour','BoP_leatherback_turtles'};
stations={'sinclair_head','ohau_head','mana_island','makara'  ,'karori_rock','wairewa_lake_forsyth','Taumutu' ,'Gulf_Harbour','BoP_leatherback_turtles'};
lat_obss=[-41.366885     ,-42.256063 ,-41.083327   ,-41.203393, -41.344623  ,-43.84110             ,-43.866527,-36.675007    ,-37.26];
lon_obss=[174.716761     , 173.853341,174.763960   ,174.702391, 174.650430  ,172.71835             ,172.381380,174.842037    ,177.59];

scase='kaitorete'; % 'graham-harrington'; %dixon-anderson';
[tstation,lat_obss,lon_obss]=read_wave_stations(scase);
stations=tstation;

k=0;

for s=1:length(stations)

  k=k+1;
  station=stations{s};

  if plot_fig==1
    figure('position',scrsz,'color',[1 1 1],'visible','on')
    hold on
    set(gca,'fontsize',12,'fontweight','bold')
  end

  lon_obs=lon_obss(s);
  lat_obs=lat_obss(s);

  expt=expts{2};
  path_expt=[path_efs,'/',expt,'/']; % GLOBALWAVE/'];%2018/01/05/00'];
  path_matlab=['/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/',expt,'/']; % GLOBALWAVE/'];%2018/01/05/00'];
  path_dm=[path_matlab,'matlab/'];
  if ck_l_file==1
    filename=[path_dm,station,'_',datestr(time_lima(1),'YYYYmmDDHH'),'_',datestr(time_lima(end),'YYYYmmDDHH'),'.mat'];
  else
    filename='bla'; system(['mkdir -p ',filename]);
  end

  display(['Looking: ',filename])%,'time_mod','model')

  if check_file==1 && exist([filename])==2
  
    display(['Loading: ',filename])%,'time_mod','model')
    load([filename])%,'time_mod','model')
  
  else
  
    wlv_mod=[];vw_mod=[];
    ucur_mod=[];vcur_mod=[];
    uw_mod=[];vw_mod=[];
    time_csm=[]; 
    hs_mod=[];tp_mod=[];
    pd_mod=[];time_mod=[];
    tm01_mod=[];tm02_mod=[];
    dm_mod=[];tm_mod=[];
    
    % Loop in model time
    for t=time_lima

      fileday=[path_dm,station,'_',datestr(t,'YYYYmmDDHH'),'.mat'];
      display(['Looking: ',fileday])%,'time_mod','model')

      if check_file==1 && exist([fileday])==2
      
        display(['Loading: ',fileday])%,'time_mod','model')
        load([fileday])%,'time_mod','model')
      
      else
    
        ptime=datestr(t,'YYYY/mm/DD/HH/');
        ftime=datestr(t,'YYYYmmDDHH');
    
        expt=expts{2};
        path_expt=[path_efs,'/',expt,'/']; % GLOBALWAVE/'];%2018/01/05/00'];
        pname=[path_expt,ptime];
        fname=[pname,prefixes{2},ftime,'-utc_',gnames{2},'.nc'];
        display(['Loading wave : ',fname]);
        if t==time_lima(1)
          lon_mod=double(ncread(fname,'lon'));
          lat_mod=double(ncread(fname,'lat'));
          [dif ilon]=nanmin(abs(lon_mod-lon_obs));
          [dif ilat]=nanmin(abs(lat_mod-lat_obs));
          dis=sqrt((lon_mod(ilon)-lon_obs).^2 + (lat_mod(ilat)-lat_obs).^2)*111;
          display(['Lon obs = ',num2str(lon_obs),' and ','Lat obs = ',num2str(lat_obs)])
          display(['Lon mod = ',num2str(lon_mod(ilon)),' and ','Lat mod = ',num2str(lat_mod(ilat))])
          display(['Distance between obs and grid point is: ',num2str(dis),' km'])
        end
        %display(['Reading time']); tic;
        timee=squeeze(double(ncread(fname,'time')))./24+t; timee=timee(1:48); %toc
        %display(['Reading Hs']); tic;
        hs=squeeze(double(ncread(fname,'hsig',          [ilon ilat 1],[1 1 49]))); hs=hs(1:end-1); %toc
        tp=squeeze(double(ncread(fname,'tpeak',         [ilon ilat 1],[1 1 49]))); tp=tp(1:end-1);
        ltime=49;
        tm01=squeeze(double(ncread(fname,'tmean01',     [ilon ilat 1],[1 1 ltime]))); tm01=tm01(1:end-1);
        %tm02=squeeze(double(ncread(fname,'tmean02',     [ilon ilat 1],[1 1 ltime]))); tm02=tm02(1:end-1);
        tm=squeeze(double(ncread(fname,'tmean',         [ilon ilat 1],[1 1 ltime]))); tm=tm(1:end-1);
        pd=squeeze(double(ncread(fname,'peak_direction',[ilon ilat 1],[1 1 ltime]))); pd=pd(1:end-1);
        dm=squeeze(double(ncread(fname,'mean_direction',[ilon ilat 1],[1 1 ltime]))); dm=dm(1:end-1);
        %ds=squeeze(double(ncread(fname,'directional_spread',[ilon ilat 1],[1 1 ltime]))); ds=ds(1:end-1);
        try
          vinfo = ncinfo(fname, 'wlv'); wlv_exist=1;
        catch; wlv_exist=0; disp('wlv not found');
        end
        if wlv_exist==1; %strcmp(expt,'NZWAVE-HR')
          wlv=squeeze(double(ncread(fname,'wlv',          [ilon ilat 1],[1 1 ltime-1])));
          ucur=squeeze(double(ncread(fname,'ucur',        [ilon ilat 1],[1 1 ltime-1])));
          vcur=squeeze(double(ncread(fname,'vcur',        [ilon ilat 1],[1 1 ltime-1])));
        else
          wlv=nan(size(hs));
          ucur=nan(size(hs));
          vcur=nan(size(hs));
        end


        
        if load_atm==1
          % Wind data
          get_nzcsm=0;
          if get_nzcsm==1
            expt=expts{1};
            path_expt=[path_efs,'/',expt,'/']; % GLOBALWAVE/'];%2018/01/05/00'];
            pname=[path_expt,ptime];
            fname=[pname,prefixes{1},ftime,'-utc_',gnames{1},'.nc']; % NZCSM
          else
            expt=expts{2};
            path_expt=[path_efs,'/',expt,'/']; % GLOBALWAVE/'];%2018/01/05/00'];
            pname=[path_expt,ptime];
            fname=[pname,prefixes{2},ftime,'-utc_',gnames{2},'.nc'];
          end
          display(['Loading weather: ',fname]);
          if t==time_lima(1)
            if get_nzcsm==1
              lon_mod=double(ncread(fname,'longitude'));
              lat_mod=double(ncread(fname,'latitude'));
              dis=sqrt((lon_mod-lon_obs).^2 + (lat_mod-lat_obs).^2)*110;
              [ilon,ilat]=find(dis==min(dis(:)));
              display(['Distance between obs and grid point is: ',num2str(dis(ilon,ilat)),' km'])
            end
            %[dif ilon]=nanmin(abs(lon_mod-lon_obs));
            %[dif ilat]=nanmin(abs(lat_mod-lat_obs));
            %dis=sqrt((lon_mod(ilon)-lon_obs).^2 + (lat_mod(ilat)-lat_obs).^2)*110;
            %lon_m=squeeze(double(ncread(fname,'longitude',          [ilon ilat],[1 1]))); 
          end

          if get_nzcsm==1
            timew=squeeze(double(ncread(fname,'time0')))./24+t; timew=timew(1:25-1);
            uw=squeeze(double(ncread(fname,'sfc_zonal_wind',          [ilonc ilatc 1],[1 1 25]))); uw=uw(1:end-1);
            vw=squeeze(double(ncread(fname,'sfc_merid_wind',          [ilonc ilatc 1],[1 1 25]))); vw=vw(1:end-1);
          else
            timew=timee; % squeeze(double(ncread(fname,'time')))./24+t; timew=timew(1:48);
            uw=squeeze(double(ncread(fname,'uwnd',          [ilon ilat 1],[1 1 49]))); uw=uw(1:end-1);
            vw=squeeze(double(ncread(fname,'vwnd',          [ilon ilat 1],[1 1 49]))); vw=vw(1:end-1);
          end
 
        end

        % saving daily file
        display(['Saving: ',fileday])
        save(fileday,'timee','hs','tp','tm','tm01','pd','dm','timew','uw','vw','wlv','ucur','vcur','ilon','ilat')

        %model.data(:,1)=tp_mod;    model.data(:,2)=pd_mod;    model.data(:,3)=uw_mod;    model.data(:,4)=vw_mod;    model.data(:,5)=nan;    model.data(:,6)=hs_mod;    model.data(:,7)=nan;    model.data(:,8)=tm01_mod; % model.data(:,8)=dm_mod;
        %time_mod=[time_mod;timee];
        %%display(['Concatenating']); tic;
        %hs_mod=[hs_mod;hs];
        %tp_mod=[tp_mod;tp];
        %tm_mod=[tm_mod;tm];
        %tm01_mod=[tm01_mod;tm01];
        %%tm02_mod=[tm02_mod;tm02];
        %pd_mod=[pd_mod;pd];
        %dm_mod=[dm_mod;dm];

      end % fileday

      time_mod=[time_mod;timee];
      %display(['Concatenating']); tic;
      hs_mod=[hs_mod;hs];
      tp_mod=[tp_mod;tp];
      tm_mod=[tm_mod;tm];
      tm01_mod=[tm01_mod;tm01];
      %tm02_mod=[tm02_mod;tm02];
      pd_mod=[pd_mod;pd];
      dm_mod=[dm_mod;dm];
      %toc;
      if exist('ucur','var')==1
        wlv_mod=[wlv_mod;wlv];
        ucur_mod=[ucur_mod;ucur];
        vcur_mod=[vcur_mod;vcur];
      end
      time_csm=[time_csm;timew];
      uw_mod=[uw_mod;uw];
      vw_mod=[vw_mod;vw];

      fclose('all');

    end %     for t=time_lima

    %csvhead={'date','PeakPeriod','PeakDirection','Directional spread','Tm01','Tm02','Hs','Wave Energy'};
    csvhead={'date','PeakPeriod','PeakDirection' ,'Uwind(m/s)'       ,'Vwind(m/s)','None1','Hs','None2','MeanDirection'};
    model.data(:,1)=tp_mod;    model.data(:,2)=pd_mod;    model.data(:,3)=uw_mod;    model.data(:,4)=vw_mod;    model.data(:,5)=nan;    model.data(:,6)=hs_mod;    model.data(:,7)=nan; model.data(:,8)=dm_mod; %  model.data(:,8)=tm01_mod; % 
    if load_atm==1
      csm.data(:,3)=uw_mod;
      csm.data(:,4)=vw_mod;
    else
      csm=[];
    end
    %lon_wave=lon_mod(ilon); lat_wave=lat_mod(ilat);

    expt=expts{2};
    path_matlab=['/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/',expt,'/']; % GLOBALWAVE/'];%2018/01/05/00'];
    path_dm=[path_matlab,'matlab/'];
    system(['mkdir -p ',path_dm]);
    display(['Saving: ',filename]);
    save([filename],'time_mod','model','time_csm','csm')%,'lon_wave','lat_wave')

  end


  if save_nc==1
    fda=[filename(1:end-4),'.nc'];
    display(['Saving: ',fda]);
    if exist(fda)==2
       system(['rm -rf ',fda]);
    end
    display(['Writing: lonlat ',fda])
    time=time_mod-datenum(1990,1,1);
    % creating nc and variables
    %nccreate(fda,'lon','Dimensions', {'lon',  length(lon_mod(ilon))});
    %nccreate(fda,'lat','Dimensions', {'lat',  length(lat_mod(ilat))});
    nccreate(fda,'time','Dimensions', {'time',  length(time)});
    nccreate(fda,'hs','Dimensions', {'time',  length(time)}); %,'unit',{ " significant wave height (m)"});
    nccreate(fda,'tp','Dimensions', {'time',  length(time)}); %,'unit',{ " significant wave height (m)"});
    nccreate(fda,'pd','Dimensions', {'time',  length(time)}); %,'unit',{ " significant wave height (m)"});
    ncwriteatt(fda,'time','unit',"days since 1990-01-01 00:00:00 (UTC)");
    ncwriteatt(fda,'hs','Description',"significant wave height (m)");
    ncwriteatt(fda,'tp','Description',"peak wave period (s)");
    ncwriteatt(fda,'pd','Description',"peak wave direction (degrees). 0 = north");
    ncwrite(fda,'time',time)
    ncwrite(fda,'hs',squeeze(model.data(:,6)))
    ncwrite(fda,'tp',squeeze(model.data(:,1)))
    ncwrite(fda,'pd',squeeze(model.data(:,2)))
    ncdisp(fda)
  end

  if save_csv==1
    fda=[filename(1:end-4),'.csv'];
    display(['Working on: ',fda]);
    if exist(fda)==2
       system(['rm -rf ',fda]);
    end

    %model.data(:,3)=nan;
    %model.data(:,4)=nan;
    %model.data(:,5)=nan;
    %model.data(:,7)=nan;

    csvhead={'date','PeakPeriod','PeakDirection','Hs'};
    if time_lima(1)==datenum(2018,1,6,0,0,0);
      time_mo=time_mod:1/24:time_mod(end);
      tp_mod=interp1(time_mo,tp_mod,time_mod);
      pd_mod=interp1(time_mo,pd_mod,time_mod);
      hs_mod=interp1(time_mo,hs_mod,time_mod);
      csvdata=[tp_mod',pd_mod',hs_mod'];

    elseif strcmp(station,'Banks_Peninsula')
      csvdata=[tp_mod,pd_mod,hs_mod];

    else;
      csvdata=[model.data(:,1),model.data(:,2),model.data(:,6)]; 

    end

    csvdata=num2cell(csvdata);
    for i=1:length(time_mod); csvtime{i,1}=datestr(time_mod(i),'yyyy-mm-ddTHH:MM:SS'); end
    %csvhead={'date','PeakPeriod','PeakDirection','DirectionalSpread','Tm01'      ,'Tm02','Hs','MeanDirection'};
    %csvhead={'date','PeakPeriod','PeakDirection' ,'Uwind(m/s)'       ,'Vwind(m/s)','None1','Hs','None2','MeanDirection'};
    %csvhead={'date','PeakPeriod','PeakDirection' ,'Uwind(m/s)'       ,'Vwind(m/s)','None1','Hs','MeanPeriod','MeanDirection'};

    % concatenating str matrices
    csvall=[csvtime,csvdata];
    %csvall=[csvhead;csvall];
    %csvwrite(fda,csvall);
    T = cell2table(csvall,'VariableNames',csvhead);
    display(['Saving: ',fda]);
    writetable(T,fda)
    fclose('all');
  end


  if plot_fig==1

    %model.data(:,1)=tp_mod;    model.data(:,2)=pd_mod;    model.data(:,3)=uw_mod;    model.data(:,4)=vw_mod;    model.data(:,5)=nan;    model.data(:,6)=hs_mod;    model.data(:,7)=nan; model.data(:,8)=dm_mod; %   
    %csm.data(:,3)=uw_mod;
    %csm.data(:,4)=vw_mod;

    dcol=[1,4,5]; %[6,1,2]; % Hs, Tp, dir
    dcol=[3,6,1]; % wind,hs,mean period
%   csvhead={'PeakPeiod'       ,'PeakDirection'      ,'Uwind(m/s)'                  ,'Vwind(m/s)'    ,'None1','Hs'                       ,'MeanDirection'};
    tnames={'wave peak period','wave peak direction',    'wind speed and direction','wind direction','Tm02' ,'wave height and direction','Mean Period (s)'};
    ynames={'(seconds)',       'From direction (^o)',    'Wind speed (knots)',     'Wind direction' ,'Tm02' ,'Significant Wave Height'  ,'Mean Period (s)'};
    time_mod=time_mod+0.5;
    time_csm=time_csm+0.5;
    time=julian(datevec(time_mod)); 

    %hs_mod=model.data(:,dcol(2))';
    %plot(time,hs_mod,'color','k','linewidth',2)
    %return

    for i=1:length(dcol)
 
      subplot(3,1,i);
      %set(gca,'fontsize',12,'fontweight','bold')
      hold on
      if dcol(i)==3

        if load_atm==1
        uw_mod=csm.data(:,3); vw_mod=csm.data(:,4); 
        mag_csm=sqrt(uw_mod.^2+vw_mod.^2).*1.94384; % from m/s to knots
        dir_csm=atan2d(uw_mod,vw_mod); 
        dir_csm=dir_csm+180; % From direction
        dir_csm(dir_csm<0)=dir_csm(dir_csm<0)+360;      

        %yyaxis left
        time=julian(datevec(time_csm)); %time_lima=julian(datevec(time_lima)); 
        %time=time_csm;

        plot(time,mag_csm,'color','k','linewidth',2)
        ylabel(['Wind speed (knots)'])

        yyaxis right
        plot(time,dir_csm,'.','color',[0 .5 0],'linewidth',2)
        ylim([0 360])
        set(gca,'YColor',[0 .5 0]);
        ylabel(['Winds from direction (^o)'])

        yyaxis left
        set(gca,'YColor','k');
        end

      elseif dcol(i)==6
        time=julian(datevec(time_mod)); 
        %time=time_mod;
        hs_mod=model.data(:,dcol(i))';
        plot(time,hs_mod,'color','k','linewidth',2)
        ylabel(['Wave height (m)'])

        yyaxis right
        plot(time,model.data(:,2)','.','color',[0 .5 0],'linewidth',2)
        ylim([0 360])
        set(gca,'YColor',[0 .5 0]);
        ylabel(['Waves from direction (^o)'])

        yyaxis left
        set(gca,'YColor','k');

        %xlim([time_lima(1) time_lima(end)+1])
        %datetick('x','dd/mmm/yyyy','keeplimits')

        %return
      else % if dcol(i)==8
        time=julian(datevec(time_mod)); 
        %time=time_mod;
        pd_mod=model.data(:,dcol(i))';
        plot(time,pd_mod,'.','color','k','linewidth',2)
        ylabel([ynames{dcol(i)}])

        %plot(time_mod,model.data(:,dcol(i))','.','color',colors{ig},'linewidth',2)

      end
      %xlim([time_lima(1) time_lima(end)+1])
      %datetick('x','dd/mmm/yyyy','keeplimits')
      if length(time_lima)<10
        gregaxh(time,12)
      else
        gregaxd(time,(time_lima(end)-time_lima(1))/30) 
      end
    
      if set_xlim==1;
        %xlim([time_lima(1)-2 time_lima(end)+2])
      end
    
      %if i==3
        %datetick('x','dd/mmm/yyyy','keeplimits')
      %end

      title([replace(tstation{s},'_',' '),' ',tnames{dcol(i)}])
      grid('on')
    
    end

    if save_fig==1
      path_dm=[path_fig,expt,'/'];
      system(['mkdir -p ',path_dm]);
      figname=[path_dm,station,'_',datestr(time_lima(1),'YYYYmmDDHH'),'_',datestr(time_lima(end),'YYYYmmDDHH'),''];
      display(['Saving: ',figname]);
      saveas(gcf,[figname,'.fig'],'fig')
      export_fig(gcf,figname,'-png','-r150');
      %close
    end

  end % if plot_fig==1



  if plot_rose

    %model.data(:,1)=tp_mod;    model.data(:,2)=pd_mod;    model.data(:,3)=uw_mod;    model.data(:,4)=vw_mod;    model.data(:,5)=nan;    model.data(:,6)=hs_mod;    model.data(:,7)=nan; model.data(:,8)=dm_mod; %   
    %csm.data(:,3)=uw_mod;
    %csm.data(:,4)=vw_mod;
      
    legh={'Obs'}; legt={'Obs'};


    i=1;

%   csvhead={'PeakPeiod'       ,'PeakDirection'      ,'Uwind(m/s)'                  ,'Vwind(m/s)'    ,'None1','Hs'                       ,'MeanDirection'};
    tnames={'Wave peak period','Wave peak direction',    'wind speed (knots)','wind direction','Tm02' ,'sig. wave height (Hs)','Mean Period (s)'};
    ynames={'(seconds)',       'From direction (^o)',    'Wind speed (knots)',        'Wind direction','Tm02' ,'Significant Wave Height'  ,'Mean Period (s)'};
    %time_mod=time_mod+0.5;
    %time_csm=time_csm+0.5;
    %tnames={'Peak Period',    'Peak Direction',    'Directional spread',   'Tm01','Tm02','Significant Wave Height','Qp','Tm','Water Level'};
    ynames={'Tp',              'Dp',                     'Ws',              'Wind dir'      ,'Tm02' ,'Hs (m)',                       'Tm'             ,'Tm','Wlv'};
    vnames={'Tp',              'Dp',                     'Ws',              'Wind dir'      ,'Tm02' ,'Hs',                       'Tm'             ,'Tm','Wlv'};

    %hs_mod=model.data(:,dcol(2))';
    %plot(time,hs_mod,'color','k','linewidth',2)
    %return
 
    %subplot(length(stations),1,s);

    %set(gca,'fontsize',12,'fontweight','bold')
    
    if dcol(i)==3

      uw_mod=csm.data(:,3); vw_mod=csm.data(:,4); 
      mag_csm=sqrt(uw_mod.^2+vw_mod.^2).*1.94384; % from m/s to knots
      dir_csm=atan2d(uw_mod,vw_mod); 
      dir_csm=dir_csm+180; % From direction
      dir_csm(dir_csm<0)=dir_csm(dir_csm<0)+360;      

      data_dir=dir_csm;
      data_val=mag_csm;


    elseif dcol(i)==6

      %time=time_mod;
      data_val=model.data(:,dcol(i));
      data_dir=model.data(:,2); % ','.','color',[0 .5 0],'linewidth',2)

    end
    

    % wave rose options        
    Option = {'vWinds', [rlevels],'ndirections',12,'anglenorth', 0, 'angleeast', 90, 'labels', {'N (0°)', 'E (90°)', 'S (180°)', 'W (270°)'}, 'freqlabelangle', 45, 'legendtype', 0};    
    %Option={'legendtype', 0, 'vWinds', [0:6]};

    %tnames={'Peak Period',    'Peak Direction',    'Directional spread',   'Tm01','Tm02','Significant Wave Height','Qp','Tm','Water Level'};
    %vnames={'Tp',             'Dp',                'Spr',                  'Tm01','Tm02','Hs',                         'Qp' ,'Tm','Wlv'};
    %for i=1:length(dcol)

    portrait=0;
    if portrait
      scrsz=[2    42   958   953];
    else
      scrsz=[1 41 1920 962];
      scrsz=[459   313   872   501];
    end

    h=figure;
    set(h,'position',scrsz,'color',[1 1 1],'visible','on')
    %set(h,'color',[1 1 1],'visible','on')
     
    ax=subplot(1,1,1);
    set(gca,'fontsize',12,'fontweight','bold')

    %Options = [Option, {'axes', ax, 'cmap', 'invbone','TitleString', {['',replace(station,'_',' '), ' ',tnames{dcol(1)}];''},...
    %'LegendType', 0, 'LegendPosition', 'north','LegendOrientation', 'vertical', 'LabLegend',vnames{dcol(1)} , 'LegendVariable', vnames{dcol(1)}, 'scalefactor', 1.0}];%, 'vWinds', [0:.5 1 1.2 2 5]}];

    %if ke==1
    %  [dif is]=nanmin(abs(time_obs-time_lima(1)));
    %  [dif ie]=nanmin(abs(time_obs-time_lima(end)));
    %  dp_obs=obs(is:ie,2)';
    %  data_obs=obs(is:ie,dcol(i))';
    %  [figure_handle, count, speeds, directions, Table] = WindRose(dp_obs,data_obs,Options);
    %  %title('Obs')
    %end
    %if portrait
    %  ax=subplot(lsub,csub,ke+1);
    %else
    %  ax=subplot(csub,lsub,ke+1);
    %end
    %lett={'(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)','(m)','(n)','(o)','(p)','(q)','(r)','(s)','(t)','(u)','(v)','(w)','(x)','(y)','(z)'};
    ltype=2; % if length(runs)==ke; ltype=2; end

    Options = [Option, {'axes', ax, 'cmap', jet,'TitleString', {[replace(station,'_',' '),' - ',tnames{dcol(1)}];''},... 
    'LegendType', ltype, 'LegendPosition', 'westoutside','LegendOrientation','vertical', 'LabLegend', ynames{dcol(1)} , 'LegendVariable', vnames{dcol(1)}, 'scalefactor', 1.0}];
    [figure_handle, count, speeds, directions, Table] = WindRose(data_dir',data_val',Options);
     
  
    
    path_dm=strcat(path_fig,expt,'/');
    system(['mkdir -p ',path_dm]);
    figname=[path_dm,station,'_rose_',replace(vnames{dcol(1)},' ','_'),'_',datestr(time_lima(1),'YYYY_mm_DD_HH_MM'),'_',datestr(time_lima(end),'YYYY_mm_DD_HH_MM'),'.png'];
    display(['Plotting: ',figname]);
    if save_fig==1
      display(['Saving: ',figname]);
      %display(['Adjust legend and press enter ....']); pause
      export_fig(gcf,figname,'-png','-r150');
      %close
      %clf('reset')
    end

  end % if plot_rose


  if plot_hist

    %model.data(:,1)=tp_mod;    model.data(:,2)=pd_mod;    model.data(:,3)=uw_mod;    model.data(:,4)=vw_mod;    model.data(:,5)=nan;    model.data(:,6)=hs_mod;    model.data(:,7)=nan; model.data(:,8)=dm_mod; %   
    %csm.data(:,3)=uw_mod;
    %csm.data(:,4)=vw_mod;
      
    legh={'Obs'}; legt={'Obs'};

    i=1;

%   csvhead={'PeakPeiod'       ,'PeakDirection'      ,'Uwind(m/s)'                  ,'Vwind(m/s)'    ,'None1','Hs'                       ,'MeanDirection'};
    tnames={'Wave peak period','Wave peak direction',    'Wind speed (knots)','wind direction','Tm02' ,'sig. wave height (m)','Mean Period (s)'};
    ynames={'(seconds)',       'From direction (^o)',    'Wind speed (knots)',        'Wind direction','Tm02' ,'Significant Wave Height (m)'  ,'Mean Period (s)'};
    %time_mod=time_mod+0.5;
    %time_csm=time_csm+0.5;
    %tnames={'Peak Period',    'Peak Direction',    'Directional spread',   'Tm01','Tm02','Significant Wave Height','Qp','Tm','Water Level'};
    vnames={'Tp',              'Dp',                     'Ws',              'Wind dir'      ,'Tm02' ,'Hs',                       'Tm'             ,'Tm','Wlv'};

    %hs_mod=model.data(:,dcol(2))';
    %plot(time,hs_mod,'color','k','linewidth',2)
    %return
 
    %subplot(length(stations),1,s);

    %set(gca,'fontsize',12,'fontweight','bold')
    
    if dcol(i)==3

      uw_mod=csm.data(:,3); vw_mod=csm.data(:,4); 
      mag_csm=sqrt(uw_mod.^2+vw_mod.^2).*1.94384; % from m/s to knots
      dir_csm=atan2d(uw_mod,vw_mod); 
      dir_csm=dir_csm+180; % From direction
      dir_csm(dir_csm<0)=dir_csm(dir_csm<0)+360;      

      data_dir=dir_csm;
      data_val=mag_csm;


    elseif dcol(i)==6

      %time=time_mod;
      data_val=model.data(:,dcol(i));
      data_dir=model.data(:,2); % ','.','color',[0 .5 0],'linewidth',2)


    end
    

    %Option={'legendtype', 0, 'vWinds', [0:6]};

    %tnames={'Peak Period',    'Peak Direction',    'Directional spread',   'Tm01','Tm02','Significant Wave Height','Qp','Tm','Water Level'};
    %vnames={'Tp',             'Dp',                'Spr',                  'Tm01','Tm02','Hs',                         'Qp' ,'Tm','Wlv'};
    %for i=1:length(dcol)

    portrait=0;
    if portrait
      scrsz=[2    42   958   953];
    else
      scrsz=[1 41 1920 962];
    end

    h=figure;
    set(h,'position',scrsz,'color',[1 1 1],'visible','on')
    %set(h,'color',[1 1 1],'visible','on')
     
    ax=subplot(1,1,1);
		set(gca,'fontsize',24,'fontweight','bold')      
    hold on
    h = histogram(data_val,rlevels, 'Normalization', 'probability');
    ylabel('Probability of ocurrence [0,1]')
		xlabel(ynames(dcol)) 
    title([replace(station,'_',' ')])%,' - ',tnames{dcol(1)}])

    nanmax(data_val)

    xlim([rlevels(1) rlevels(end)])
    %set(gca,'xtick',rlevels)

    %end
    
    path_dm=strcat(path_fig,expt,'/');
    system(['mkdir -p ',path_dm]);
    figname=[path_dm,station,'_hist_',replace(vnames{dcol(1)},' ','_'),'_',datestr(time_lima(1),'YYYY_mm_DD_HH_MM'),'_',datestr(time_lima(end),'YYYY_mm_DD_HH_MM'),'.png'];
    display(['Plotting: ',figname]);
    if save_fig==1
      display(['Saving: ',figname]);
      export_fig(gcf,figname,'-png','-r150');
      %close
      %clf('reset')

      fda=[filename(1:end-4),'_hist_data.csv'];
      display(['Working on: ',fda]);
      if exist(fda)==2
         system(['rm -rf ',fda]);
      end
      xlsdata=[h.BinEdges(2:end)',h.Values'];
      writetable(cell2table(num2cell(xlsdata)),fda);
    end

    write_table=1;

    if write_table==1;

      fda=[filename(1:end-4),'_hist_data_',replace(vnames{dcol(1)},' ','_'),'.csv'];
      display(['Working on: ',fda]);
      if exist(fda)==2
         system(['rm -rf ',fda]);
      end
      xlsdata=[h.BinEdges(2:end)',round(100.*h.Values,5)'];
      writetable(cell2table(num2cell(xlsdata)),fda);
    end

  end % if plot_hist


end % for s=1:length(stations)


%legend(['Obs',expts],'location','best')


end % for dcol=[6 3]
