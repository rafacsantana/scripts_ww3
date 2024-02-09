clear
close_all=1;
if close_all==1; close all; end
warning off

% Daily time
time_lima=datenum(1985,1,1,0,0,0):1:datenum(1997,08,31,0,0,0); % GLOBALWAVE-ERA5
time_lima=datenum(1997,06,08,0,0,0):1:datenum(1997,08,31,0,0,0); % GLOBALWAVE-ERA5
%time_lima=datenum(2008,1,1,0,0,0):1:datenum(2009,7,30,0,0,0); % 
%time_lima=datenum(2015,1,1,0,0,0):1:datenum(2016,12,31,0,0,0); % 
%time_lima=datenum(1988,6,28,0,0,0):1:datenum(1988,7,2,0,0,0); % NZWAVE-ERA5
%time_lima=datenum(1989,6,1,0,0,0):1:datenum(1989,6,3,0,0,0); % NZWAVE-ERA5
%time_lima=datenum(1985,1,1,0,0,0):1:datenum(1996,12,31,0,0,0); % 
%time_lima=time_lima+3;

%time_lima=datenum(2014,7,9,0,0,0):1:datenum(2014,9,3,0,0,0); % initial date
%time_lima=datenum(2018,1,6,0,0,0):1:datenum(2023,11,15,0,0,0); % Graham Harrington ECAN request on 18/09/2023

%time_lima=datenum(2014,8,5,0,0,0):1:datenum(2014,8,17,0,0,0); % 9-m wave event
%time_lima=datenum(2018,1,20,0,0,0):1:datenum(2020,12,1,0,0,0); % crossing between hindcast and ecoconnect
%time_lima=datenum(2021,1,9,0,0,0):1:datenum(2021,1,16,0,0,0); % Baring head wave direction
%time_lima=datenum(2021,5,20,0,0,0):1:datenum(2021,6,1,0,0,0); % Baring Head tide and swell
%time_lima=datenum(2021,10,19,0,0,0):1:datenum(2021,10,24,0,0,0); % Banks Peninsula head tide and swell
%time_lima=datenum(2021,10,1,0,0,0):1:datenum(2021,10,24,0,0,0); % Banks Peninsula data assimilation
time_lima=datenum(2021,6,26,12,0,0):1/24:datenum(2021,6,30,0,0,0); %
%time_lima=datenum(2021,6,29,6,0,0):1/24:datenum(2021,6,30,0,0,0); %
%time_lima=datenum(2021,1,5,0,0,0):1:datenum(2021,1,15,0,0,0); %

%time_lima=datenum(2029,12,1,0,0,0):1:datenum(2030,1,1,0,0,0); % Firts two days of output

stations=[1,2];
files={'Banks_Peninsula','Baring_Head','wairewa_lake_forsyth','Pelorus_Sound','Taumutu'};% 'SteepHead_SP_FullRecord_QC';

runs=[10];
%runs=[2,3,4,8];
runs=[1:5];
runs=[1,2];

expt_names={'GLOBALWAVE'          ,'NZWAVE'             ,'NZWAVE-ST6'         ,'NZWAVE-HR-NOTIDES','NZWAVE-HR',... 5
            'GLOBALWAVE-GFDL-CCAM','NZWAVE-ST4-GLOBALUM','NZWAVE-ST6-GLOBALUM','NZWAVE-GFDL-CCAM' ,'NZWAVE-ERA5',...10
            'GLOBALWAVE-ERA5'};            

plot_series=0;
plot_map   =1;
plot_rose  =0;
plot_disp  =0;
plot_eva   =0;

plot_aserie=0;
plot_amap  =0;

plot_obs   =1;
proc_obs   =0;
ck_mod_mat =1;
save_csv   =0;
set_xlim   =1;
save_fig   =0;
save_video =1;

correct_model=1;

% series dispersion
% 1 {'PeakPerid'}   2 {'PeakDirection'}   3 {'Directional spread'}   4 {'Tm01'}   5 {'Tm02'}   6 {'Hs'}   7 {'Qp'}  8 Tm
% 9 wlv
dcol=[1,4,5,9]; %[6,1,2]; % Tp, Tm01, Tm02
%dcol=[6,5,8,9]; %[6,1,2]; % Tp, Tm01, Tm02
dcol=[6,1,2]; % Hs, Tp, dir
%dcol=[6,1,8]; % Hs, Tp, dir
dcol=[6]; % Hs, Tp, dir

% maps
pcname ='hs'; % m_pcolor: wlv, hs
vecname='hs'; % m_vec: hs, curr

% altimeter stats
acol=[1:6];
acol=[7:12];
avals={'hs_mean_rmse','hs_maxi_rmse','hs_mean_bias','hs_maxi_bias','hs_mean_corr','hs_maxi_corr',...
       'hs_mean_frmse','hs_maxi_frmse','hs_mean_fbias','hs_maxi_fbias','hs_mean_fcorr','hs_maxi_fcorr'};
avalm={'hs_mean_mrmse','hs_maxi_mrmse','hs_mean_mbias','hs_maxi_mbias','hs_mean_mcorr','hs_maxi_mcorr',... 
       'hs_mean_fmrmse','hs_maxi_fmrmse','hs_mean_fmbias','hs_maxi_fmbias','hs_mean_fmcorr','hs_maxi_fmcorr'}; % 


colors={'k','b','r','m','c','y'};

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

if plot_map
  %fname=[pname,ww3pre,'_',ftime,'-utc_',gnames{ig},'.nc'];
  fname='/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/NZWAVE-HR/2021/01/01/00/ww3g_2021010100-utc_nzwave_hr+nzcsm.nc';
  depth_hr=double(ncread(fname,'depth',[1 1 1],[Inf Inf 1]));
  lon_hr=double(ncread(fname,'lon'));
  lat_hr=double(ncread(fname,'lat'));
  
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

scrsz=[1 1 1366 768];
scrsz=[1 1 1920 1080];
%scrsz=get(0,'screensize');

figure('position',scrsz,'color',[1 1 1],'visible','on')
hold on
set(gca,'fontsize',12,'fontweight','bold')

%m_proj('lambert','long', [min(lon_hr(:))-1 max(lon_hr(:))+1],'lat',[min(lat_hr(:))-1 max(lat_hr(:))+1]);


% biggest loop on stations
for fobs=stations

  file=files{fobs};
  file=replace(file,' ','_');
  
  % OBS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  file_obs=[path_obs,file];
  [obs,time_obs,lon_obs,lat_obs]=proc_wave_obs(file,proc_obs); %,plot_obs);   

  lon_obss=[lon_obss;lon_obs];
  lat_obss=[lat_obss;lat_obs];

  if fobs==stations(end)
 
    if plot_map
        
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


        if save_video==1 & plot_map==1;
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

    end % if plot_map

  end % fobs==stations(end)

end % for fobs=stations

if save_video==1 & plot_map==1;
  close(vidObj);
end

display(['End of script']);
