clear
close all

% Daily time
time_lima=datenum(2023,8,9,0,0,0):1:datenum(2023,8,13,0,0,0); % MPI request on 28/08/2023

check_file=1;
set_xlim  =0;
proc_obs  =0;
save_fig  =1;

path_fig=['/scale_wlg_nobackup/filesets/nobackup/niwa03150/santanarc/figures/']; % GLOBALWAVE/'];%2018/01/05/00'];

path_efs='/niwa/archive/ecoconnect/EFS/';

prefixes={'tn_','ww3g_'};
gnames  ={'nzcsm','nzwave_hr+nzcsm'};

%tn_2023081000-utc_nzcsm.nc

%/niwa/archive/ecoconnect/EFS/NZCSM/2023/08/10/00
%/niwa/archive/ecoconnect/EFS/NZWAVE-HR/2023/08/10/00

expts={'NZCSM','NZWAVE-HR'};

scrsz=[1 1 1366 768];
scrsz=[2 42 958 953];
%scrsz=get(0,'screensize');
  
colors={'k','k'};

tstation={'Sinclair Head','Ohau Head','Mana Island','Makara','Karori Rock'};
stations={'sinclair_head','ohau_head','mana_island','makara','karori_rock'};

lat_obss=[-41.366885,-42.256063,  -41.083327, -41.203393, -41.344623];

lon_obss=[174.716761, 173.853341, 174.763960, 174.702391, 174.650430];

k=0;
for s=1:length(stations)
  k=k+1;
  station=stations{s};

  figure('position',scrsz,'color',[1 1 1],'visible','on')
  hold on
  set(gca,'fontsize',12,'fontweight','bold')


  lon_obs=lon_obss(k);
  lat_obs=lat_obss(k);

  expt='NZWAVE-HR';
  path_expt=[path_efs,'/',expt,'/']; % GLOBALWAVE/'];%2018/01/05/00'];
  path_matlab=['/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/',expt,'/']; % GLOBALWAVE/'];%2018/01/05/00'];
  path_dm=[path_matlab,'matlab/'];
  filename=[path_dm,station,'_',datestr(time_lima(1),'YYYYmmDDHH'),'_',datestr(time_lima(end),'YYYYmmDDHH'),'.mat'];

  if check_file==1 && exist([filename])==2
  
    display(['Loading: ',filename])%,'time_mod','model')
    load([filename])%,'time_mod','model')
  
  else
  
    uw_mod=[];vw_mod=[];
    time_csm=[]; 
    hs_mod=[];tp_mod=[];
    pd_mod=[];time_mod=[];
    
    % Loop in model time
    for t=time_lima
    
      ptime=datestr(t,'YYYY/mm/DD/HH/');
      ftime=datestr(t,'YYYYmmDDHH');
    
      expt='NZWAVE-HR';
      path_expt=[path_efs,'/',expt,'/']; % GLOBALWAVE/'];%2018/01/05/00'];
      pname=[path_expt,ptime];
      fname=[pname,prefixes{2},ftime,'-utc_',gnames{2},'.nc'];
      display(['Loading: ',fname]);
      if t==time_lima(1)
        lon_mod=double(ncread(fname,'lon'));
        lat_mod=double(ncread(fname,'lat'));
        [dif ilon]=nanmin(abs(lon_mod-lon_obs));
        [dif ilat]=nanmin(abs(lat_mod-lat_obs));
        dis=sqrt((lon_mod(ilon)-lon_obs).^2 + (lat_mod(ilat)-lat_obs).^2)*110;
        display(['Distance between obs and grid point is: ',num2str(dis),' km'])
      end
      timee=squeeze(double(ncread(fname,'time')))./24+t; timee=timee(1:48);
      hs=squeeze(double(ncread(fname,'hsig',          [ilon ilat 1],[1 1 49]))); hs=hs(1:end-1);
      tp=squeeze(double(ncread(fname,'tpeak',         [ilon ilat 1],[1 1 49]))); tp=tp(1:end-1);
      pd=squeeze(double(ncread(fname,'peak_direction',[ilon ilat 1],[1 1 49]))); pd=pd(1:end-1);
      time_mod=[time_mod;timee];
      hs_mod=[hs_mod;hs];
      tp_mod=[tp_mod;tp];
      pd_mod=[pd_mod;pd];

      % Wind data
      expt='NZCSM';
      path_expt=[path_efs,'/',expt,'/']; % GLOBALWAVE/'];%2018/01/05/00'];
      pname=[path_expt,ptime];
      fname=[pname,prefixes{1},ftime,'-utc_',gnames{1},'.nc'];
      display(['Loading: ',fname]);
      if t==time_lima(1)
        lon_mod=double(ncread(fname,'longitude'));
        lat_mod=double(ncread(fname,'latitude'));
        dis=sqrt((lon_mod-lon_obs).^2 + (lat_mod-lat_obs).^2)*110;
        [ilonc,ilatc]=find(dis==min(dis(:)));
        %[dif ilon]=nanmin(abs(lon_mod-lon_obs));
        %[dif ilat]=nanmin(abs(lat_mod-lat_obs));
        %dis=sqrt((lon_mod(ilon)-lon_obs).^2 + (lat_mod(ilat)-lat_obs).^2)*110;
        display(['Distance between obs and grid point is: ',num2str(dis(ilonc,ilatc)),' km'])
        %lon_m=squeeze(double(ncread(fname,'longitude',          [ilon ilat],[1 1]))); 
      end
      timew=squeeze(double(ncread(fname,'time0')))./24+t; timew=timew(1:25-1);
      uw=squeeze(double(ncread(fname,'sfc_zonal_wind',          [ilonc ilatc 1],[1 1 25]))); uw=uw(1:end-1);
      vw=squeeze(double(ncread(fname,'sfc_merid_wind',          [ilonc ilatc 1],[1 1 25]))); vw=vw(1:end-1);
    
      time_csm=[time_csm;timew];
      uw_mod=[uw_mod;uw];
      vw_mod=[vw_mod;vw];
    
    end
    
    model.data(:,6)=hs_mod;
    model.data(:,1)=tp_mod;
    model.data(:,2)=pd_mod;
    csm.data(:,3)=uw_mod;
    csm.data(:,4)=vw_mod;
    

    expt='NZWAVE-HR';
    path_matlab=['/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/',expt,'/']; % GLOBALWAVE/'];%2018/01/05/00'];
    path_dm=[path_matlab,'matlab/'];
    system(['mkdir -p ',path_dm]);
    display(['Saving: ',filename]);
    save([filename],'time_mod','model','time_csm','csm')
  
  end


  dcol=[1,4,5]; %[6,1,2]; % Hs, Tp, dir
  dcol=[3,6,1]; % Hs, Tp, dir
  tnames={'wave peak period',    'wave peak direction',    'wind speed and direction',   'wind direction','Tm02','wave height and direction','Qp'};
  ynames={'(seconds)',    'From direction (^o)',    'Wind speed (knots)',   'Wind direction','Tm02','Significant Wave Height','Qp'};
  time_mod=time_mod+0.5;
  time_csm=time_csm+0.5;

  for i=1:length(dcol)
 
    subplot(3,1,i);
    %set(gca,'fontsize',12,'fontweight','bold')
    hold on
    if dcol(i)==3
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
    elseif dcol(i)==1
      time=julian(datevec(time_mod)); 
      %time=time_mod;
      pd_mod=model.data(:,dcol(i))';
      plot(time,pd_mod,'.','color','k','linewidth',2)
      ylabel(['Wave period (s)'])

      %plot(time_mod,model.data(:,dcol(i))','.','color',colors{ig},'linewidth',2)

    end
    %xlim([time_lima(1) time_lima(end)+1])
    %datetick('x','dd/mmm/yyyy','keeplimits')
    gregaxh(time,12)
  
    if set_xlim==1;
      %xlim([time_lima(1)-2 time_lima(end)+2])
    end
  
    %if i==3
      %datetick('x','dd/mmm/yyyy','keeplimits')
    %end

    title([tstation{s},' ',tnames{dcol(i)}])
    grid('on')
  
  end

  if save_fig==1
    path_dm=[path_fig,expt,'/'];
    system(['mkdir -p ',path_dm]);
    figname=[path_dm,station,'_',datestr(time_lima(1),'YYYYmmDDHH'),'_',datestr(time_lima(end),'YYYYmmDDHH'),'.png'];
    display(['Saving: ',figname]);
    export_fig(gcf,figname,'-png','-r150');
    close
  end

end % for s=1:length(stations)


%legend(['Obs',expts],'location','best')



