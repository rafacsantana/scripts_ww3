clear
close_all=1;
if close_all==1; close all; end
warning off

% Daily time
time_lima=datenum(1983,1,1,0,0,0):1:datenum(1983,1,3,0,0,0); % 
time_lima=datenum(1985,1,1,0,0,0):1:datenum(1985,1,21,0,0,0); % GLOBALWAVE-ERA5
%time_lima=datenum(2001,1,1,0,0,0):1:datenum(2001,12,31,0,0,0); % GLOBALWAVE-ERA5
%time_lima=datenum(1990,07,01,0,0,0):1:datenum(1990,7,31,0,0,0); % GLOBALWAVE-ERA5
%time_lima=datenum(1993,01,01,0,0,0):1:datenum(2002,12,31,0,0,0); % GLOBALWAVE-ERA5
%time_lima=datenum(2008,1,1,0,0,0):1:datenum(2009,7,30,0,0,0); % 
%time_lima=datenum(2008,1,1,0,0,0):1:datenum(2022,12,31,0,0,0); % 
%time_lima=datenum(1988,6,28,0,0,0):1:datenum(1988,7,2,0,0,0); % NZWAVE-ERA5
%time_lima=datenum(1989,6,1,0,0,0):1:datenum(1989,6,3,0,0,0); % NZWAVE-ERA5
%time_lima=datenum(2017,1,1,0,0,0):1:datenum(2019,12,31,0,0,0); % 
%time_lima=datenum(2021,2,23,0,0,0):1:datenum(2021,12,31,0,0,0); % 
%time_lima=time_lima+3;

%time_lima=datenum(2014,7,9,0,0,0):1:datenum(2014,9,3,0,0,0); % initial date
%time_lima=datenum(2018,1,6,0,0,0):1:datenum(2023,11,15,0,0,0); % Graham Harrington ECAN request on 18/09/2023

%time_lima=datenum(2014,8,5,0,0,0):1:datenum(2014,8,17,0,0,0); % 9-m wave event
%time_lima=datenum(2018,1,20,0,0,0):1:datenum(2020,12,1,0,0,0); % crossing between hindcast and ecoconnect
%time_lima=datenum(2021,1,8,0,0,0):1:datenum(2021,1,15,0,0,0); % Baring head wave direction, time series
%time_lima=datenum(2021,1,13,0,0,0):1:datenum(2021,1,13,0,0,0); % Baring head wave direction and currents, maps
%time_lima=datenum(2021,5,20,0,0,0):1:datenum(2021,5,27,0,0,0); % Baring Head tide and swell
%time_lima=datenum(2021,10,19,0,0,0):1:datenum(2021,10,24,0,0,0); % Banks Peninsula head tide and swell
%time_lima=datenum(2021,10,20,0,0,0):1:datenum(2021,10,22,0,0,0); % Banks Peninsula head tide and direction
%time_lima=datenum(2021,10,1,0,0,0):1:datenum(2021,10,24,0,0,0); % Banks Peninsula data assimilation
%time_lima=datenum(2021,6,27,0,0,0):1:datenum(2021,7,2,0,0,0); % Banks Peninsula data assimilation
%time_lima=datenum(2021,6,26,12,0,0):1/24:datenum(2021,6,30,0,0,0); %
time_lima=datenum(2020,12,1,0,0,0):1:datenum(2020,12,22,0,0,0); %
time_lima=datenum(2021,1,1,0,0,0):1:datenum(2021,12,31,0,0,0); %

%time_lima=datenum(2029,12,1,0,0,0):1:datenum(2030,1,1,0,0,0); % Firts two days of output

stations=[2];
files={'Banks_Peninsula','Baring_Head','wairewa_lake_forsyth','Pelorus_Sound','Taumutu'};% 'SteepHead_SP_FullRecord_QC';

runs=[10];
%runs=[2,3,4,8];
runs=[1:5];

expt_names={'GLOBALWAVE'          ,'NZWAVE'             ,'NZWAVE-ST6'         ,'NZWAVE-HR-NOTIDES','NZWAVE-HR',... 5
            'GLOBALWAVE-GFDL-CCAM','NZWAVE-ST4-GLOBALUM','NZWAVE-ST6-GLOBALUM','NZWAVE-GFDL-CCAM' ,'NZWAVE-ERA5',...10
            'GLOBALWAVE-ERA5'     ,'NZWAVE-ERA5-2021'   ,'NZWAVE-ST6-DEFAULT'};            

plot_series=0;
plot_map   =0;
plot_downs =0;
plot_rose  =1;
plot_disp  =0;
plot_eva   =0;

plot_aserie=0;
plot_amap  =0;

plot_obs   =1;
proc_obs   =0;
ck_mod_mat =1;
save_csv   =0;
save_fig   =0;
save_video =0;
save_gif   =0;

set_xlim   =1;
obs_first  =1;
add_stats  =1;
portrait   =1;

correct_model=1; % -1=scatter_kde, 0=binscatter, 1=nothing, >1 correct it

% series dispersion
% 1 {'PeakPerid'}   2 {'PeakDirection'}   3 {'Directional spread'}   4 {'Tm01'}   5 {'Tm02'}   6 {'Hs'}   7 {'Qp'}  8 Tm
% 9 wlv 10 ucurr 11 vcurr

dcol=[1,4,5,9]; %[6,1,2]; % Tp, Tm01, Tm02
%dcol=[6,5,8,9]; %[6,1,2]; % Tp, Tm01, Tm02
dcol=[6,1,2]; % Hs, Tp, dir
%dcol=[6,1,8]; % Hs, Tp, dir
dcol=[6]; % Hs, wlv
%dcol=[6]; % Hs, wlv

% maps
pcname ='hs'; % m_pcolor: wlv, hs, curr
vecname='hs'; % m_vec: hs, curr

% altimeter stats
acol=[1:6];
acol=[7:12];
avals={'hs_mean_rmse','hs_maxi_rmse','hs_mean_bias','hs_maxi_bias','hs_mean_corr','hs_maxi_corr',...
       'hs_mean_frmse','hs_maxi_frmse','hs_mean_fbias','hs_maxi_fbias','hs_mean_fcorr','hs_maxi_fcorr'};
avalm={'hs_mean_mrmse','hs_maxi_mrmse','hs_mean_mbias','hs_maxi_mbias','hs_mean_mcorr','hs_maxi_mcorr',... 
       'hs_mean_fmrmse','hs_maxi_fmrmse','hs_mean_fmbias','hs_maxi_fmbias','hs_mean_fmcorr','hs_maxi_fmcorr'}; % 

colors={'c','b','r','m','k','c'};

for i=1:length(runs)
  expts{i}=expt_names{runs(i)};
end

%sufixes in netcdf files
gnames={'globalwave+globalum','nzwave+nzlam','nzwave_hr+nzcsm','tongawave+globalum','nzwave_era5'};

path_source=['/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/']; % GLOBALWAVE/'];%2018/01/05/00'];
path_santanarc='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/';
path_obs=[path_santanarc,'data/obs/'];
path_fig='/scale_wlg_nobackup/filesets/nobackup/niwa03150/santanarc/figures/';

% stations or grids - stations are bad in Baring_Head
ww3pre={'ww3p_interp','ww3g'};
ww3pre=ww3pre{2};

% Loading satellite obs
if plot_aserie || plot_amap
  atype='cmems_l4';
  % cmems_l4  
  if strcmp(atype,'cmems_l4')
    path_cop='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/data/copernicus-l4-wave-data/';
    file_cop=[path_cop,'cmems_obs-wave_glo_phy-swh_my_multi-l4-2deg_P1D_multi-vars_143.00E-185.00E_53.00S-21.00S_2002-01-15-2020-12-31.nc'];
    time_cop=ncread(file_cop,'time')+datenum('1950-01-01');
    lon_cop=ncread(file_cop,'longitude');
    lat_cop=ncread(file_cop,'latitude');
    hs_cop_mean=ncread(file_cop,'VAVH_DAILY_MEAN');
    hs_cop_maxi=ncread(file_cop,'VAVH_DAILY_MAX');
  end

  if plot_amap
    %fname=[pname,ww3pre,'_',ftime,'-utc_',gnames{ig},'.nc'];
    fname='/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/NZWAVE/2021/01/01/00/ww3g_2021010100-utc_nzwave+nzlam.nc';
    depth_mod=double(ncread(fname,'depth',[1 1 1],[Inf Inf 1]));
    lon_mod=double(ncread(fname,'lon'));
    lat_mod=double(ncread(fname,'lat'));
  end
end

if plot_map || plot_downs
  %fname=[pname,ww3pre,'_',ftime,'-utc_',gnames{ig},'.nc'];
  fname='/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/NZWAVE/2021/01/01/00/ww3g_2021010100-utc_nzwave+nzlam.nc';
  depth_mod=double(ncread(fname,'depth',[1 1 1],[Inf Inf 1]));
  lon_mod=double(ncread(fname,'lon'));
  lat_mod=double(ncread(fname,'lat'));
  
  %% gebco bathy 
  %file_geb='/scale_wlg_persistent/filesets/project/niwa00020/data/topography/gebco/1-min/gridone.nc';
  %depg    =ncread(file_geb,'z');
  %long    =ncread(file_geb,'x_range');
  %latg    =ncread(file_geb,'y_range');
  %dgeb    =ncread(file_geb,'spacing');
  %lon_geb =long(1):dgeb:long(2);
  %lat_geb =latg(1):dgeb:latg(2);
  %depth_geb=reshape(depg,length(lon_geb),length(lat_geb));
  %%iz=find(lon_geb==0); lon_ge=[lon_geb(iz:end),lon_geb(1:iz-1)]
  %lon_geb=wrapTo360(lon_geb);

  %lon_geb=lon_geb(1:10:end);
  %lat_geb=lat_geb(1:10:end);
  %depth_geb=depth_geb(1:10:end,1:10:end);

  %[dif ilo]=nanmin(abs(min(lon_mod(:))-lon_geb));  [dif flo]=nanmin(abs(max(lon_mod(:))-lon_geb));
  %[dif ila]=nanmin(abs(min(lat_mod(:))-lat_geb));  [dif fla]=nanmin(abs(max(lat_mod(:))-lat_geb));
  %lon_geb=lon_geb(ilo:flo);
  %lat_geb=lat_geb(ila:fla);
  %depth_geb=depth_geb(ilo:flo,ila:fla);

  % gebco bathy
  gebco_nc='/scale_wlg_persistent/filesets/project/niwa00020/data/topography/gebco/1-min/gridone.nc';
  spac=double(ncread(gebco_nc,'spacing'));
  lon_gebco=double(ncread(gebco_nc,'x_range')); lon_gebco=lon_gebco(1):spac(1):lon_gebco(end);
  lat_gebco=double(ncread(gebco_nc,'y_range')); lat_gebco=lat_gebco(end):-spac(2):lat_gebco(1);
  bat_gebco=double(ncread(gebco_nc,'z'));
  bat_gebco=reshape(bat_gebco,length(lon_gebco),length(lat_gebco));
  lon_gebco=lon_gebco(1:10:end);
  lat_gebco=lat_gebco(1:10:end);
  bat_gebco=bat_gebco(1:10:end,1:10:end);
  lon_gebco=wrapTo360(lon_gebco);

  %[dif ilo]=nanmin(abs(min(lon_mod(:))-lon_geb));  [dif flo]=nanmin(abs(max(lon_mod(:))-lon_geb));
  %[dif ila]=nanmin(abs(min(lat_mod(:))-lat_geb));  [dif fla]=nanmin(abs(max(lat_mod(:))-lat_geb));
  %lon_geb=lon_geb(ilo:flo);
  %lat_geb=lat_geb(ila:fla);
  %depth_geb=depth_geb(ilo:flo,ila:fla);

end

lon_obss=[]; lat_obss=[];


% biggest loop on stations
for fobs=stations

  file=files{fobs};
  file=replace(file,' ','_');
  
  % OBS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  file_obs=[path_obs,file];
  [obs,time_obs,lon_obs,lat_obs]=proc_wave_obs(file,proc_obs); %,plot_obs);   
 
  if plot_series
  
    scrsz=[1 1 1366 768];
    %scrsz=get(0,'screensize');
    
    if close_all==1; figure('position',scrsz,'color',[1 1 1],'visible','on'); end
    hold on
    set(gca,'fontsize',12,'fontweight','bold')
    
    if obs_first
      legh={'Obs'}; legt={'Obs'};
    else
      legh=[]; legt=[];
    end

    ke=0;
    for expt=expts
      ke=ke+1;
    
      expt=expt{1};
      [ig,ltime]=grab_nc_sufix(expt);

    
      path_expt=[path_source,expt,'/']; % GLOBALWAVE/'];%2018/01/05/00'];
      
      path_dm=[path_expt,'matlab/'];
      system(['mkdir -p ',path_dm]);
      %filename=[path_dm,file,'_',datestr(time_lima(1),'YYYYmmDDHH'),'_',datestr(time_lima(end),'YYYYmmDDHH'),'.mat'];
      hs_mod=[];tp_mod=[]; pd_mod=[];time_mod=[];
      tm01_mod=[];tm_mod=[]; pd_mod=[]; we_mod=[];
      tm02_mod=[];ds_mod=[]; wlv_mod=[];time_mod=[];
      ucur_mod=[];vcur_mod=[]; 

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
    
            dis=sqrt((lon_mod(ilon)-lon_obs).^2 + (lat_mod(ilat)-lat_obs).^2)*111;
            display(['Distance between obs ',file,' and grid point is: ',num2str(dis),' km'])

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
            if strcmp(expt,'NZWAVE-HR')
              wlv=squeeze(double(ncread(fname,'wlv',          [ilon ilat 1],[1 1 ltime-1]))); 
              model(ke).wlv=wlv;
              ucur=squeeze(double(ncread(fname,'ucur',        [ilon ilat 1],[1 1 ltime-1]))); 
              model(ke).ucur=ucur;
              vcur=squeeze(double(ncread(fname,'vcur',        [ilon ilat 1],[1 1 ltime-1]))); 
              model(ke).vcur=vcur;
            else
              wlv=nan(size(hs));
              ucur=nan(size(hs));
              vcur=nan(size(hs));
            end
          end
          we=(1/16)*1025*9.8*(hs.^2);    

          display(['Saving: ',filename])
          save(filename,'time','tp','pd','ds','tm01','tm02','hs','tm','wlv','ucur','vcur','we')
  
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
        if exist('ucur','var')==1
          ucur_mod=[ucur_mod;ucur];
	  	    vcur_mod=[vcur_mod;vcur];
		      we_mod=[we_mod;we];
        end  
        clear model
     
        % 1 {'PeakPerid'}   2 {'PeakDirection'}   3 {'Directional spread'}   4 {'Tm01'}   5 {'Tm02'}   6 {'Hs'}   7 {'Qp'}  8 Tm
        model.data(:,9)=wlv_mod; % we_mod
        model.data(:,8)=tm_mod;
        model.data(:,1)=tp_mod;
        model.data(:,2)=pd_mod;
        model.data(:,3)=ds_mod;
        model.data(:,4)=tm01_mod;
        model.data(:,5)=tm02_mod;
        model.data(:,6)=hs_mod;
		    model.data(:,7)=nan;
        if exist('ucur','var')==1
		      model.data(:,10)=ucur_mod;
		      model.data(:,11)=vcur_mod;
		      model.data(:,12)=we_mod;
	      end
  
        %path_dm=[path_expt,'matlab/'];
        %system(['mkdir -p ',path_dm]);
        %display(['Saving: ',filename])
        %save(filename,'time_mod','model')
      
      end % for t=time_lima
    
      if save_csv==1
    
        fda=[filename(1:end-14),datestr(time_lima(1),'YYYYmmDDHH'),'_',datestr(time_lima(end),'YYYYmmDDHH'),'.csv'];
        display(['Saving: ',fda]);
    
        if exist(fda)==2
        system(['rm -rf ',fda]);
        end
    
        model.data(:,7)=nan;
        csvdata=model.data;
        csvdata=num2cell(csvdata);
        for i=1:length(time_mod); csvtime{i,1}=datestr(time_mod(i),'yyyy-mm-ddTHH:MM:SS'); end
        csvhead={'date','PeakPeriod','PeakDirection','Directional spread','Tm01','Tm02','Hs','Qp','Tm','Water level','U current','V current','Wave Energy (J/m2)'};
    
        % concatenating str matrices
        csvall=[csvtime,csvdata];
        %csvall=[csvhead;csvall];
        %csvwrite(fda,csvall);
        T = cell2table(csvall,'VariableNames',csvhead);
        writetable(T,fda)
        fclose('all');

      end
    
      % Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      tnames={'Peak Period',    'Peak Direction',    'Directional spread',   'Tm01','Tm02','Significant Wave Height (Hs)','Qp','Tm','Water Level','U current','V current','Wave Energy (J/m2)'};
      for i=1:length(dcol)
        
        subplot(length(dcol),1,i);
        set(gca,'fontsize',12,'fontweight','bold')
        hold on
        if ke==1 & plot_obs==1
          if obs_first==1;
	    if dcol(i)<=size(obs,2) 
              plot(time_obs,obs(:,dcol(i))','.','markersize',12,'color',[0 .7 0],'linewidth',2)
	    end
          end
          for eix=1:length(expts); plot(nan,nan,'.','markersize',20,'color',colors{eix},'linewidth',2); end
          plot(nan,nan,'.','markersize',20,'color',[0 .7 0],'linewidth',2)
        end
        plot(time_mod,model.data(:,dcol(i))','.','markersize',12,'color',colors{ke},'linewidth',2)
        if obs_first==0; 
          if ke==length(expts) & plot_obs==1
	    if dcol(i)<=size(obs,2) 
              plot(time_obs,obs(:,dcol(i))','.','markersize',12,'color',[0 .7 0],'linewidth',2)
	    end
          end
        end

        if set_xlim==1;
          xlim([time_lima(1) time_lima(end)])
        end
      
        %if i==3
          datetick('x','dd/mmm/yyyy','keeplimits')
        %end
        title([replace(file,'_',' '),' ',tnames{dcol(i)}])
        grid('on')
        
        modeli=interp1(time_mod,model.data(:,dcol(i)),time_obs);
        
        if (dcol(i)~=2 || dcol(i)~=1) & i==1
          mod_rmse=nanrmse(obs(:,dcol(i)),modeli');
          mod_corr=nancorr(obs(:,dcol(i)),modeli');
          mod_bias=nanmean(modeli')-nanmean(obs(:,dcol(i)));
          if add_stats
            legh=[legh,{[expt,' rmse=',num2str(mod_rmse,'%.2f'),' corr=',num2str(mod_corr,'%.2f'),', bias=',num2str(mod_bias,'%.2f')]}];
          else
            legh=[legh,{[expt]}];%,' rmse=',num2str(mod_rmse,'%.2f'),' corr=',num2str(mod_corr,'%.2f'),', bias=',num2str(mod_bias,'%.2f')]}];
          end
	elseif dcol(i)==2 || dcol(i)==1
	  ylim([0 360])
        end
        
      end % for i=1:length(dcol)
    
    end % for expt=expts

    if obs_first==0; 
      legh=[legh,{['Obs.']}]; %,' rmse=',num2str(mod_rmse,'%.2f'),' corr=',num2str(mod_corr,'%.2f'),', bias=',num2str(mod_bias,'%.2f')]}];
    end

    %if ~isempty(find(dcol==6))
    %  subplot(length(dcol),1,find(dcol==6));
      subplot(length(dcol),1,1);
      legend([legh],'location','best')
    %end
    
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

        [ig,ltime]=grab_nc_sufix(expt);
 
      
        path_expt=[path_source,expt,'/']; % GLOBALWAVE/'];%2018/01/05/00'];
    
        pname=[path_expt,ptime];
        
        fname=[pname,ww3pre,'_',ftime,'-utc_',gnames{ig},'.nc'];
        %fname=[pname,'ww3g_',ftime,'-utc_',gnames{1},'.nc'];
        display(['Loading: ',fname]);
    
        if checkin==0; % t==time_lima(1)
          checkin=1; 
          if strcmp(expt,'NZWAVE-HR')
            depth_mod=double(ncread(fname,'depth',[1 1 1],[Inf Inf 1]));
            wlv_mod=double(ncread(fname,'wlv',[1 1 1],[Inf Inf 1]));
            depth_mod=depth_mod-wlv_mod;
          end

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
      dp_obs=270-atan2(v_obsc,u_obsc)*180/pi;
      dp_obs(dp_obs>360)=dp_obs(dp_obs>360)-360;

      pd_obs=dp_obs;
      hs_obs=sqrt(u_obsc.^2+v_obsc.^2);
      %hs_obs=interp1(time_obs,obs(:,6),time);
      tp_obs=interp1(time_obs,obs(:,1),time);
        
      for i=1:length(time)
 
        ke=0;
        for expt=expts
          ke=ke+1;

          display(['Ploting map']); tic
          ax(ke)=subplot(1,length(expts),ke);
          hold on
          m_proj('lambert','long', [lon_obs-.3 lon_obs+.25],'lat',[lat_obs-.25 lat_obs+.25]);

          if strcmp(pcname,'hs')
            data=model(ke).hs(:,:,i);
            cmap=cmocean('thermal');
            caxis([-2 2])
          elseif strcmp(pcname,'wlv')
            data=model(ke).wlv(:,:,i);
            pd_mod=model(ke).pd(:,:,i);
            pd=mod(-90-pd_mod,360);
            u_mod=cosd(pd).*data; v_mod=sind(pd).*data;
            cmap=cmocean('balance');
            caxis([-2 2])
          elseif strcmp(pcname,'curr')
            u_mod=model(ke).ucur(:,:,i);
            v_mod=model(ke).vcur(:,:,i);
            data=sqrt(u_mod.^2+v_mod.^2);
            cmap=cmocean('tempo');
          end
          %colormap(ax(1),gwbreal)
          %m_pcolor(lon_ccmp,lat_ccmp,wsc_ccmp');
          %set(get(cs,'ylabel'),'string','N/m^3','fontsize',20,'fontweight','bold');
          %caxis([-2E-06 2E-06])
          colormap(ax(ke),cmap) % parula, avhrr
          m_pcolor(lon_mod,lat_mod,data');
          %caxis([0 .13])
          cb=colorbar;%('southoutside');
          set(get(cb,'ylabel'),'string','','fontsize',12,'fontweight','bold');
          set(cb,'fontsize',12,'fontweight','bold');

          if strcmp(pcname,'hs')
            caxis([0 4])
          elseif strcmp(pcname,'wlv')
            caxis([-2 2])
          elseif strcmp(pcname,'curr')
            caxis([0 2])
          end

          %%set(cs,'position','south');
          %shading interp
          %%caxis([nanmin(t_avhrr(:)) nanmax(t_avhrr(:))]); % 12 25

          %[cs,h]=m_contour(lon_bath,lat_bath,bath',[-1000 -1000],'color',[1 1 1],'linewidth',2);
          %clabel(cs,h,'color',[1 1 1],'fontsize',14,'fontweight','bold','LabelSpacing',2000)
          if strcmp(expt,'NZWAVE-HR')
            [cs,h]=m_contour(lon_mod,lat_mod,-depth_mod',[-2 -20 -50 -100 -200],'color',[1 1 1],'linewidth',2);
            clabel(cs,h,'color',[0 0 0],'fontsize',10,'fontweight','bold','LabelSpacing',2000)
          else
            [cs,h]=m_contour(lon_geb,lat_geb,-depth_geb',[-200 -1000],'color',[1 1 1],'linewidth',2);
            clabel(cs,h,'color',[0 0 0],'fontsize',10,'fontweight','bold','LabelSpacing',2000)
          end

          %[cs,h]=m_contour(lon_aviso,lat_aviso,z_aviso_avg,'color','k');
          %[cs,h]=m_contour(lon_aviso,lat_aviso,z_aviso_avg,[.76:.04:.88],'color','k');
          %clabel(cs,h,'color',[0 0 0],'fontsize',14,'fontweight','bold','LabelSpacing',2000) % each contour line represents 0.02 increase/decrease

          m_gshhs_f('patch',[.7 .7 .7],'edgecolor','k','linewidth',.5);
          %%m_grid('fontname','helvetica','fontsize',22,'fontweight','bold');
          %%m_grid('xtick',[round(loni):2:round(lonf)],'xticklabels',[round(loni):2:round(lonf)],'fontname','helvetica','fontsize',20,'fontweight','bold');
          m_grid('fontname','helvetica','fontsize',14,'fontweight','bold');
          title([replace(file,'_',' '),' Hs and ',pcname,' at ',datestr(time(i),'HH:MM dd/mmm/yyyy')],'fontsize',14,'fontweight','bold')

          %set(gca,'fontsize',14,'fontname','helvetica','fontweight','bold')
          %%m_contour(plon,plat,bathy',[200 2000],'w','linewidth',1);
          %%caxis([0 .5])

          % getting u and v
          if strcmp(vecname,'wlv')
          elseif strcmp(vecname,'hs')
            dat=model(ke).hs(:,:,i);
            pd_mod=model(ke).pd(:,:,i);
            pd=mod(-90-pd_mod,360);
            u_mod=cosd(pd).*dat; v_mod=sind(pd).*dat;
          elseif strcmp(vecname,'curr')
            u_mod=model(ke).ucur(:,:,i);
            v_mod=model(ke).vcur(:,:,i);
          end

          v_spa=1;
          if strcmp(vecname,'hs')
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
      	  vsize=2;
          m_vec(vsize,lon_mod(ilon),lat_mod(ilat),u_mod(ilon,ilat),v_mod(ilon,ilat),'b','shaftwidth',2,'headwidth',6,'headangle',60,'headlength',8)
          [hp ht]=m_vec(vsize,lon_obs+.1,lat_obs+.12,2,0,'b',                           'shaftwidth',2,'headwidth',6,'headangle',60,'headlength',8,'key',[num2str(2.0),' m']);
          set(ht,'fontsize',14,'fontweight','bold')
          
          % obs single vector 
          %pd=mod(-90-pd_obs,360);
          %u_obs=cosd(pd).*hs_obs; v_obs=sind(pd).*hs_obs;
          m_vec(vsize,lon_obs,lat_obs,u_obsc(i),v_obsc(i),'m'                          ,'shaftwidth',2,'headwidth',6,'headangle',60,'headlength',8)
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
          %clf('reset')	
		end

        if save_video==1 

          if i==1 & t==time_lima(1)
            videoName = [path_dm,'/',file,'_',pcname,'_',vecname,'_',datestr(time_lima(1),'yyyymmdd_HHMMSS'),'_',datestr(time_lima(end),'yyyymmdd_HHMMSS'),''];
            display(['Saving video name: ',videoName])
            vidObj = VideoWriter(videoName,'Uncompressed AVI');
            %vidObj = VideoWriter(videoName,'Motion JPEG 2000');
            set(vidObj,'FrameRate',5)
            open(vidObj);
          end
          display(['Printing frame time: ',datestr(time(i),'yyyymmdd_HHMMSS')]);
          pause(0.5)
          M = getframe(gcf);
          writeVideo(vidObj,M);

        end

        if save_gif==1 
          gifName = [path_dm,'/',file,'_',pcname,'_',vecname,'_',datestr(time_lima(1),'yyyymmdd_HHMMSS'),'_',datestr(time_lima(end),'yyyymmdd_HHMMSS'),'.gif'];
          display(['Saving GIF: ',gifName]);
          drawnow
          frame = getframe(1);
          im = frame2im(frame);
          [imind,cm] = rgb2ind(im,256);
          if i == 1 & t==time_lima(1);
            imwrite(imind,cm,gifName,'gif', 'Loopcount',inf,'DelayTime',.2);
          else
            imwrite(imind,cm,gifName,'gif','WriteMode','append','DelayTime',.2);
          end
		end

        if save_fig==1 || save_video==1 || save_gif==1
          clf('reset')
          set(gcf,'color',[1 1 1])
	else
	  pause(0.5)
	end
 
      end % for i=1:length(time)

    end % for t=time_lima
    
    if save_video==1 
      close(vidObj);
    end

    %close 

  end % if plot_map

  if plot_downs

    lon_obss=[lon_obss;lon_obs];
    lat_obss=[lat_obss;lat_obs];

    if fobs==stations(1)

      scrsz=[1 1 1366 768];
      scrsz=[1 1 1920 1080];
      %scrsz=get(0,'screensize');
      
      figure('position',scrsz,'color',[1 1 1],'visible','on')
      hold on
      set(gca,'fontsize',12,'fontweight','bold')

    elseif fobs==stations(end)
        
      legh={'Obs'}; legt={'Obs'};

      checkin=0;  
      % Loop in desired time
      for t=time_lima%(1)
       
        ptime=datestr(t,'YYYY/mm/DD/00/');
        ftime=datestr(t,'YYYYmmDD00');
      
        ax(1)=subplot(1,1,1);
        ke=0;

        for expt=expts
          ke=ke+1;
        
          expt=expt{1};

          [ig,ltime]=grab_nc_sufix(expt);
        
          path_expt=[path_source,expt,'/']; % GLOBALWAVE/'];%2018/01/05/00'];
      
          pname=[path_expt,ptime];
          
          fname=[pname,ww3pre,'_',ftime,'-utc_',gnames{ig},'.nc'];
          %fname=[pname,'ww3g_',ftime,'-utc_',gnames{1},'.nc'];
          display(['Loading: ',fname]);

          % getting time
          time=squeeze(double(ncread(fname,'time')))./24+floor(t); time=time(1:ltime-1);
          % finding the right hour
          [dif i]=nanmin(abs(time-t)); %:length(time)

          if strcmp(expt,'NZWAVE-HR')
            %depth_mod=double(ncread(fname,'depth',[1 1 1],[Inf Inf 1]));
            %wlv_mod=double(ncread(fname,'wlv',[1 1 1],[Inf Inf 1]));
            %depth_mod=depth_mod-wlv_mod;
          end

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
          
          display(['Reading nc variables']); tic
          hs=squeeze(double(ncread(fname,'hsig',          [1 1 i],[Inf Inf 1]))); 
          tp=squeeze(double(ncread(fname,'tpeak',         [1 1 i],[Inf Inf 1]))); 
          %tm01=squeeze(double(ncread(fname,'tmean01',     [1 1 1],[Inf Inf ltime-1]))); 
          %tm=squeeze(double(ncread(fname,'tmean',         [1 1 1],[Inf Inf ltime-1]))); 
          pd=squeeze(double(ncread(fname,'peak_direction',[1 1 i],[Inf Inf 1]))); 
          if strcmp(expt,'NZWAVE-HR')
            wlv=squeeze(double(ncread(fname,'wlv',          [1 1 i],[Inf Inf 1]))); 
            model(ke).wlv=wlv;
            ucur=squeeze(double(ncread(fname,'ucur',        [1 1 i],[Inf Inf 1]))); 
            model(ke).ucur=ucur;
            vcur=squeeze(double(ncread(fname,'vcur',        [1 1 i],[Inf Inf 1]))); 
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

          time=model(ke).time;
          tnames={'Peak Period',    'Peak Direction',    'Directional spread',   'Tm01','Tm02 or Tm','Significant Wave Height (Hs)','Qp'};

          % converting to u and v before interpolating
          pd_obs=obs(:,2); hs_obs=obs(:,6);
          pd=mod(-90-pd_obs,360); u_obs=cosd(pd).*hs_obs; v_obs=sind(pd).*hs_obs;
          u_obsl=interp1(time_obs,u_obs,time); v_obsl=interp1(time_obs,v_obs,time);
          u_obsc=interp1(time_obs,u_obs,time,'cubic'); v_obsc=interp1(time_obs,v_obs,time,'cubic');
          dp_obs=270-atan2(v_obsc,u_obsc)*180/pi;
          dp_obs(dp_obs>360)=dp_obs(dp_obs>360)-360;

          pd_obs=dp_obs;
          hs_obs=sqrt(u_obsc.^2+v_obsc.^2);
          %hs_obs=interp1(time_obs,obs(:,6),time);
          tp_obs=interp1(time_obs,obs(:,1),time);
            
          %for i=1%:length(time)
 
            %ke=0;
            %for expt=expts
            %  ke=ke+1;

          display(['Ploting map']); tic
          if strcmp(pcname,'hs')
            data=model(ke).hs;%(:,:,i);
            cmap=avhrr; % cmocean('matter');
            %cmap=cmocean('ice');
            cmap = mycmap(cmap,1,1,14);
          elseif strcmp(pcname,'wlv')
            data=model(ke).wlv;%(:,:,i);
            pd_mod=model(ke).pd;%(:,:,i);
            pd=mod(-90-pd_mod,360);
            u_mod=cosd(pd).*data; v_mod=sind(pd).*data;
            cmap=cmocean('balance');
            caxis([-2 2])
          elseif strcmp(pcname,'curr')
            u_mod=model(ke).ucur;%(:,:,i);
            v_mod=model(ke).vcur;%(:,:,i);
            data=sqrt(u_mod.^2+v_mod.^2);
            cmap=cmocean('thermal');
          end

          hold on
          m_proj('lambert','long', [143.3203125-2 184.570495605469+2],'lat',[-54.453125-2 -20.8589458465576+2]);

          %colormap(ax(1),gwbreal)
          %m_pcolor(lon_ccmp,lat_ccmp,wsc_ccmp');
          %set(get(cs,'ylabel'),'string','N/m^3','fontsize',20,'fontweight','bold');
          %caxis([-2E-06 2E-06])
          colormap(ax(1),cmap) % parula, avhrr
          m_pcolor(lon_mod,lat_mod,data');
          %m_contourf(lon_mod,lat_mod,data',[0:.5:8]);
          %m_contour(lon_mod,lat_mod,data',[0:.5:8]);
          %caxis([0 .13])

          cb=colorbar;%('southoutside');
          set(get(cb,'ylabel'),'string','Sig. Wave Height (m)','fontsize',12,'fontweight','bold');
          set(cb,'fontsize',12,'fontweight','bold');

          if strcmp(pcname,'hs')
            caxis([0 14])
            %caxis([0 20])
          elseif strcmp(pcname,'wlv')
            caxis([-2 2])
          elseif strcmp(pcname,'curr')
            cmap=cmocean('thermal');
          end

          %%set(cs,'position','south');
          %shading interp
          %%caxis([nanmin(t_avhrr(:)) nanmax(t_avhrr(:))]); % 12 25

          %[cs,h]=m_contour(lon_bath,lat_bath,bath',[-1000 -1000],'color',[1 1 1],'linewidth',2);
          %clabel(cs,h,'color',[1 1 1],'fontsize',14,'fontweight','bold','LabelSpacing',2000)
          if ke==length(expts)
            %[cs,h]=m_contour(lon_mod,lat_mod,-depth_mod',[-2 -20 -50 -100 -200],'color',[1 1 1],'linewidth',2);
            %clabel(cs,h,'color',[0 0 0],'fontsize',10,'fontweight','bold','LabelSpacing',2000)
            [cs,h]=m_contour(lon_gebco,lat_gebco,bat_gebco',[-200 -200],'color',[146 36 40]./255,'linewidth',2);
            %clabel(cs,h,'color',[0 0 0],'fontsize',10,'fontweight','bold','LabelSpacing',2000)
          end

          %[cs,h]=m_contour(lon_aviso,lat_aviso,z_aviso_avg,'color','k');
          %[cs,h]=m_contour(lon_aviso,lat_aviso,z_aviso_avg,[.76:.04:.88],'color','k');
          %clabel(cs,h,'color',[0 0 0],'fontsize',14,'fontweight','bold','LabelSpacing',2000) % each contour line represents 0.02 increase/decrease

          if ke==1%length(expts)
            %%m_grid('fontname','helvetica','fontsize',22,'fontweight','bold');
            %%m_grid('xtick',[round(loni):2:round(lonf)],'xticklabels',[round(loni):2:round(lonf)],'fontname','helvetica','fontsize',20,'fontweight','bold');
          %elseif ke==length(expts)
            m_gshhs_i('patch',[.7 .7 .7],'edgecolor','k','linewidth',.5);
            m_grid('fontname','helvetica','fontsize',14,'fontweight','bold');
            %title([replace(file,'_',' '),' Hs and water level at ',datestr(time(i),'HH:MM dd/mmm/yyyy')],'fontsize',14,'fontweight','bold')
            title(['Sig. Wave Height at ',datestr(time(i),'HH:MM dd/mmm/yyyy'),' UTC'],'fontsize',14,'fontweight','bold')
          end
          %set(gca,'fontsize',14,'fontname','helvetica','fontweight','bold')
          %%m_contour(plon,plat,bathy',[200 2000],'w','linewidth',1);
          %%caxis([0 .5])

          % getting u and v
          if strcmp(vecname,'wlv')
          elseif strcmp(vecname,'hs')
            dat=model(ke).hs;%(:,:,i);
            pd_mod=model(ke).pd;%(:,:,i);
            pd=mod(-90-pd_mod,360);
            u_mod=cosd(pd).*dat; v_mod=sind(pd).*dat;
          elseif strcmp(vecname,'curr')
            u_mod=model(ke).ucur;%(:,:,i);
            v_mod=model(ke).vcur;%(:,:,i);
          end

          v_spa=40;
          if strcmp(vecname,'hs')
            m_vec(50,lon_modm(10:v_spa:end,10:v_spa:end),lat_modm(10:v_spa:end,10:v_spa:end),u_mod(10:v_spa:end,10:v_spa:end)',v_mod(10:v_spa:end,10:v_spa:end)','w','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5)
            if ke==1
              [hp ht]=m_vec(50,144,-30,8,0,'w','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5,'key',[num2str(8.0),' m']);
              set(ht,'fontsize',14,'fontweight','bold')
            end
          end
          if ke==1
            m_text(144,-29.5,['Sig. Wave Height'],'color','w','fontsize',14,'fontweight','bold') 
          end

          pd_mod=model(ke).pd;%(:,:,i);
          for o=1:length(stations)
            filo=files{stations(o)};
            filo=replace(filo,'_',' ');
            m_plot(lon_obss(o),lat_obss(o),'.w','markersize',30)
            m_plot(lon_obss(o),lat_obss(o),'.r','markersize',20)
            if o==1 
              %m_text(lon_obss(o)-6,lat_obss(o)+.2,filo,'color','k','fontsize',14,'fontweight','bold') 
            elseif o==2 
              %m_text(lon_obss(o)-8,lat_obss(o)+.2,filo,'color','k','fontsize',14,'fontweight','bold') 
            end
          end
          %m_text(lon_obs+.1,lat_obs+.18,['Obs dp: ',num2str(pd_obs(i),'%.2f'),'^o'],'fontsize',14,'fontweight','bold') 

          % model single vector 
          data=model(ke).hs;%(:,:,i);
          pd_mod=model(ke).pd;%(:,:,i);
          pd=mod(-90-pd_mod,360);
          u_mod=cosd(pd).*data; v_mod=sind(pd).*data;
          %m_vec(50,lon_mod(ilon),lat_mod(ilat),u_mod(ilon,ilat),v_mod(ilon,ilat),'b','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5)
          %[hp ht]=m_vec(50,lon_obs+.1,lat_obs+.12,2,0,'b','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5,'key',[num2str(2.0),' m']);
          %set(ht,'fontsize',14,'fontweight','bold')
          
          %% obs single vector 
          %%pd=mod(-90-pd_obs,360);
          %%u_obs=cosd(pd).*hs_obs; v_obs=sind(pd).*hs_obs;
          %m_vec(50,lon_obs,lat_obs,u_obsc(i),v_obsc(i),'m','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5)
          %%m_vec(10,lon_obs,lat_obs,u_obsl(i),v_obsl(i),'c','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5)


          if ke>1%==length(expts)
            %[range,ln,lt]=m_lldist([nanmax(lon_hr) nanmin(lon_hr)],[nanmin(lat_hr) nanmin(lat_hr)],140);
            %m_line(wrapTo360(ln),lt,'color','k','linewi',2)
            %[range,ln,lt]=m_lldist([nanmin(lon_hr) nanmax(lon_hr)],[nanmax(lat_hr) nanmax(lat_hr)],140);
            %m_line(wrapTo360(ln),lt,'color','k','linewi',2)
            %[range,ln,lt]=m_lldist([nanmin(lon_hr) nanmin(lon_hr)],[nanmin(lat_hr) nanmax(lat_hr)],140);
            %m_line(wrapTo360(ln),lt,'color','k','linewi',2)
            %[range,ln,lt]=m_lldist([nanmax(lon_hr) nanmax(lon_hr)],[nanmin(lat_hr) nanmax(lat_hr)],140);
            %m_line(wrapTo360(ln),lt,'color','k','linewi',2)

            m_plot([lon_modm(:,1)],[lat_modm(:,1)],'k','linewidth',2);
            m_plot([lon_modm(:,end)],[lat_modm(:,end)],'k','linewidth',2);
            m_plot([lon_modm(1,:)],[lat_modm(1,:)],'k','linewidth',2);
            m_plot([lon_modm(end,:)],[lat_modm(end,:)],'k','linewidth',2);
          end

          display(['Finished ploting map']); toc

          %if dcol(i)==6
          %  mod_rmse=nanrmse(obs(:,dcol(i)),modeli');
          %  mod_bias=nanmean(modeli')-nanmean(obs(:,dcol(i)));
          %  legh=[legh,{[expt,' rmse=',num2str(mod_rmse,'%.2f'),', bias=',num2str(mod_bias,'%.2f')]}];
          %end

            %end % for expt=expts

          %end % for i=1:length(time)

        end % for expt=expts

        path_dm=strcat(path_fig,'paper_1/');
        system(['mkdir -p ',path_dm]);
        figname=[path_dm,file,'_',pcname,'_',vecname,'_',datestr(time(i),'YYYY_mm_DD_HH_MM'),'.png'];
        display(['Plotting: ',figname]);
        if save_fig==1
          display(['Saving: ',figname]);
          export_fig(gcf,figname,'-png','-r150');
          display(['End of script']);
          return       
          %close
          %clf('reset')
        else
          %pause
          %clf('reset')
        end


        if save_video==1 & plot_downs==1;
          if t==time_lima(1)
            videoName = [path_dm,'/',pcname,'_',datestr(time_lima(1),'yyyymmdd_HHMMSS'),'_',datestr(time_lima(end),'yyyymmdd_HHMMSS'),''];
            display(['Saving video name: ',videoName])
            vidObj = VideoWriter(videoName,'Uncompressed AVI');
            set(vidObj,'FrameRate',5)
            open(vidObj);
          end
          display(['Printing frame time: ',figname]);
          pause(0.5)
          M = getframe(gcf);
          writeVideo(vidObj,M);
        end

        clf('reset')
        set(gcf,'color',[1 1 1])

      end % for t=time_lima

    end % fobs==stations(end)

  end % if plot_downs


  if plot_rose
      
    legh={'Obs'}; legt={'Obs'};
    
    ke=0;
    for expt=expts
      ke=ke+1;
    
      expt=expt{1};


      [ig,ltime]=grab_nc_sufix(expt);
    
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
      
          %display(['Loading: ',filename])
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
            if strcmp(expt,'NZWAVE-HR')
              wlv=squeeze(double(ncread(fname,'wlv',          [ilon ilat 1],[1 1 ltime-1]))); 
              %model(ke).wlv=wlv;
              %ucur=squeeze(double(ncread(fname,'ucur',        [ilon ilat 1],[1 1 ltime-1]))); 
              %model(ke).ucur=ucur;
              %vcur=squeeze(double(ncread(fname,'vcur',        [ilon ilat 1],[1 1 ltime-1]))); 
              %model(ke).vcur=vcur;
            else
              wlv=nan(size(hs));
              ucur=nan(size(hs));
              vcur=nan(size(hs));
            end
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
        
        %path_dm=[path_expt,'matlab/'];
        %system(['mkdir -p ',path_dm]);
        %display(['Saving: ',filename])
        %save(filename,'time_mod','model')
      
      end % for t=time_lima

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

      % wave rose options        
      Option = {'vWinds', [0:6],'anglenorth', 0, 'angleeast', 90, 'labels', {'N (0°)', 'E (90°)', 'S (180°)', 'W (270°)'}, 'freqlabelangle', 45, 'legendtype', 0};    
      %Option={'legendtype', 0, 'vWinds', [0:6]};

      tnames={'Peak Period',    'Peak Direction',    'Directional spread',   'Tm01','Tm02','Significant Wave Height','Qp','Tm','Water Level'};
      vnames={'Tp',             'Dp',                'Spr',                  'Tm01','Tm02','Hs'                     ,'Qp','Tm','Wlv'};
      for i=1:length(dcol)

        if portrait
          scrsz=[2    42   958   953];
        else
          scrsz=[1 41 1920 962];
        end

        h=figure(i);
        set(h,'position',scrsz,'color',[1 1 1],'visible','on')

        if length(runs)==1
          %subplot(length(dcol),1,i);
          lsub=1; csub=1;
        elseif mod(length(runs),2)==1
          lsub=round(length(runs)/2); csub=round(length(runs)/2)-1;
        elseif mod(length(runs),2)==0
          lsub=round(length(runs)/2); csub=round(length(runs)/2);
        end
         
        if portrait
          ax=subplot(lsub,csub,ke);
        else
          ax=subplot(csub,lsub,ke);
        end
        if strcmp(file,'Banks_Peninsula')
          maxf=35;
				elseif strcmp(file,'Baring_Head')
          maxf=65;
        end
        Options = [Option, {'axes', ax, 'cmap', 'invbone','TitleString', {['(a) ',replace(file,'_',' '), ' obs. ',tnames{dcol(i)}];''},...
        'LegendType', 0, 'LegendPosition', 'north','LegendOrientation', 'vertical', 'LabLegend',vnames{dcol(i)} , 'LegendVariable', vnames{dcol(i)}, 'MaxFrequency',maxf, 'nDirections', 12, 'scalefactor', 1.0}];%, 'vWinds', [0:.5 1 1.2 2 5]}];
        if ke==1
          [dif is]=nanmin(abs(time_obs-time_lima(1)));
          [dif ie]=nanmin(abs(time_obs-time_lima(end)));
          dp_obs=obs(is:ie,2)';
          data_obs=obs(is:ie,dcol(i))';
          [figure_handle, count, speeds, directions, Table] = WindRose(dp_obs,data_obs,Options);
          %title('Obs')
          display(['Obs'])
          [directions,sum(count,2)]
        end
        if portrait
          ax=subplot(lsub,csub,ke+1);
        else
          ax=subplot(csub,lsub,ke+1);
        end
        lett={'(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)','(m)','(n)','(o)','(p)','(q)','(r)','(s)','(t)','(u)','(v)','(w)','(x)','(y)','(z)'};
        ltype=0; if length(runs)==ke; ltype=2; end
        Options = [Option, {'axes', ax, 'cmap', 'invbone','TitleString', {[lett{ke},' ',expt];''},...
        'LegendType', ltype, 'LegendPosition', 'west','LegendOrientation','vertical', 'LabLegend', vnames{dcol(i)} , 'LegendVariable', vnames{dcol(i)}, 'MaxFrequency',maxf, 'nDirections', 12, 'scalefactor', 1.0}];
        [figure_handle, count, speeds, directions, Table] = WindRose(model.data(:,2)',model.data(:,dcol(i))',Options);
        display(['Expt: ',expt])
        [directions,sum(count,2)]
         
      end
    
    end
    
    path_dm=strcat(path_fig,'paper_1/');
    system(['mkdir -p ',path_dm]);
    figname=[path_dm,file,'_wave_rose_',replace(tnames{dcol(i)},' ','_'),'_',datestr(time_lima(1),'YYYY_mm_DD_HH_MM'),'_',datestr(time_lima(end),'YYYY_mm_DD_HH_MM'),'.png'];
    display(['Plotting: ',figname]);
    if save_fig==1
      display(['Correct legend box and press enter to save your figure']);
      pause
      display(['Saving: ',figname]);
      export_fig(gcf,figname,'-png','-r150');
      %close
      %clf('reset')
    end

  end % if plot_rose


end % fobs=fstations


% compute maps of annual average wave height and plot them



