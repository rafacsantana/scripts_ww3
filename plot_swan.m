clear
close all
%clf('reset')
tic

%if strncmp(expt,'lugan',5)
time=datenum(2000,1,1,0,0,0):.5/24:datenum(2000,4,1,0,0,0);
%elseif strncmp(expt,'lv',2)
time=datenum(2020,4,1,1,0,0):.5/24:datenum(2020,4,11,0,0,0);
%time=datenum(2020,4,7,0,0,0):.5/24:datenum(2020,4,7,1,0,0);
%time=datenum(2020,4,3,0,0,0):.5/24:datenum(2020,4,8,0,0,0);
%end

runs=[1,5,8,9,10];
%runs=[3 4];
runs=[13,12];


plot_map   =0;
plot_serie =1;
plot_seriem=0;
plot_obs   =0; % obs=1 or analysis stations = 0
plot_obsm  =1; % obs stations in map

save_fig  =0;
save_video=0;

vname='Hsig'; % 'Hsig'; % Dir, Dspr, RTpeak, Tm01

expt_names={'lv_era5_tc_harold'             ,'lv_era5_tc_harold_corr'        ,'lv_tcwindgen_tc_harold_mwr'   ,'lv_tcwindgen_tc_harold_mwrcfac093',...  % 3=4; lv_tcwindgen_tc_harold_mwr = lv_tcwindgen_tc_harold_mwrcfac093
            'lv_tcwindgen_tc_harold_st6c100','lv_tcwindgen_tc_harold_st6c093','lv_tcwindgen_tc_harold_st6c93','lv_tcwindgen_tc_harold_xuz',...         % 5=6; lv_tcwindgen_tc_harold_st6c100 = lv_tcwindgen_tc_harold_st6c93
            'lv_tcwindgen_tc_harold_xuzst6' ,'lv_xuzst6_ssw_neg'             ,'SWAN'                         ,'SWAN-TCwindgen',...  % 9=10; lv_tcwindgen_tc_harold_xuzst6 = lv_xuzst6_ssw_neg
            'SWAN-ERA5'};
 
if plot_obs==1
  time_obs=[datenum(2020,4,4,6,3,54),datenum(2020,4,5,17,45,16),datenum(2020,4,6,6,0,51),datenum(2020,4,7,17,42,58)];
  hs_obs=[5.025,3.115,4.633,1.05];
  lon_obs=[166.24  ,168.722 ,165.598 ,167.792];
  lat_obs=[-14.1635,-15.4866,-16.3755,-15.9868];
  snames={'S1','S2','S3','S4'};
else
  % declustering
  lon_obs=[166.560946,167.447265,167.3129388,166.9881609];
  lat_obs=[-15.841595,-15.560540,-15.5405023,-15.6683607];
  snames={'D1','D2','Tutuba','Araki'};
end
stations=[1:length(lon_obs)];

% stations for xbeach comparison - same as Tutuba
%lonx=167.3129388;
%latx=-15.5405023;


% Figure variables
colors={'c','b','k','g','--b','--k','m','c','g'};

% winds
filewnd=('/scale_wlg_persistent/filesets/project/niwa03150/santanarc/swan/projects/VANUATU/lv_xuzst6_ssw_neg/Luganville_harold.nc');
timewnd=double(ncread(filewnd,'time'))./(24*3600)+datenum(2020,04,01,0,0,0); %:1/48:datenum(2000,01,18,0,0,0);
lonwnd =double(ncread(filewnd,'x'));
latwnd =double(ncread(filewnd,'y'));
[lonwndm,latwndm]=meshgrid(lonwnd,latwnd);  

% read TC track file
file_obs='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/TCwindgen_zhonghou/projects/lv_tcwindgen_tc_harold_mwr/TC_Harold_track_proc_final_mwr.csv';
display(['Processing: ',file_obs,'.txt']);
ob=importdata([file_obs]);
ob.textdata(1,:)
for i=1:length(ob.data)
  textdata=num2str(ob.data(i,1));
  if length(textdata)==8
    time_track(i)=datenum(textdata,'YYYYmmdd');
  else
    time_track(i)=datenum(textdata,'YYYYmmdd.HH');
  end
end
lon_track=ob.data(:,2);
lat_track=ob.data(:,3);
lon_track=interp1(time_track,lon_track,time);
lat_track=interp1(time_track,lat_track,time);


path_fig='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/swan/projects/VANUATU/figures/';
scrsz=[1 1 1366 768];
scrsz=[1 1 1920 1080];
scrsz=[1 1 1910 990];

%scrsz=get(0,'screensize');
figure('position',scrsz,'color',[1 1 1],'visible','on')
hold on
set(gca,'fontsize',12,'fontweight','bold')
      
datas=nan(length(time),length(stations),length(runs));
legh=[];

ke=0;
for ex=runs 
  ke=ke+1;
  expt=expt_names(ex);
  expt=expt{1};

  legh=[legh,{[replace(expt,'_','-')]}]; % ,' rmse=',num2str(mod_rmse,'%.2f'),' corr=',num2str(mod_corr,'%.2f'),', bias=',num2str(mod_bias,'%.2f')]}];


  matfile=['/scale_wlg_persistent/filesets/project/niwa03150/santanarc/swan/projects/VANUATU/',expt,'/Van.mat'];
  display(['Loading: ',matfile])
  load(matfile); 
  
  %return
  kt=0;
  for t=time
    kt=kt+1;

    if plot_seriem==1
      kee=0;
      legh=[];
      for ex=runs 
        kee=kee+1;
        expt=expt_names(ex);
        expt=expt{1};
        legh=[legh,{[replace(expt,'_','-')]}]; % ,' rmse=',num2str(mod_rmse,'%.2f'),' corr=',num2str(mod_corr,'%.2f'),', bias=',num2str(mod_bias,'%.2f')]}];
        matfile=['/scale_wlg_persistent/filesets/project/niwa03150/santanarc/swan/projects/VANUATU/',expt,'/Van.mat'];
        display(['Loading: ',matfile])
        load(matfile); 
        %hs=Hsig_20000101_033000
        display(['Loading: ',vname,'_',datestr(t,'yyyymmdd_HHMMSS')]);
        data=double(eval([vname,'_',datestr(t,'yyyymmdd_HHMMSS')]));
        dp=double(eval(['Dir_',datestr(t,'yyyymmdd_HHMMSS')]));
        tp=double(eval(['RTpeak_',datestr(t,'yyyymmdd_HHMMSS')]));
        uwnd=double(eval(['Windv_x_',datestr(t,'yyyymmdd_HHMMSS')]));
        vwnd=double(eval(['Windv_y_',datestr(t,'yyyymmdd_HHMMSS')]));
        if plot_serie==1 || plot_seriem==1;
          k=0;
          for s=stations
            k=k+1;
            dis=sqrt((Xp-lon_obs(s)).^2+(Yp-lat_obs(s)).^2);
            [dif loc]=nanmin(dis(:));
            datas(kt,k,kee)=data(loc);
            udatas(kt,k,kee)=uwnd(loc);
            vdatas(kt,k,kee)=vwnd(loc);
          end  
        end %if plot_serie==1
      end % for ex=runs
    else
       display(['Loading: ',vname,'_',datestr(t,'yyyymmdd_HHMMSS')]);
       data=double(eval([vname,'_',datestr(t,'yyyymmdd_HHMMSS')]));
       dp=double(eval(['Dir_',datestr(t,'yyyymmdd_HHMMSS')]));
       tp=double(eval(['RTpeak_',datestr(t,'yyyymmdd_HHMMSS')]));
       uwnd=double(eval(['Windv_x_',datestr(t,'yyyymmdd_HHMMSS')]));
       vwnd=double(eval(['Windv_y_',datestr(t,'yyyymmdd_HHMMSS')]));
          k=0;
          for s=stations
            k=k+1;
            dis=sqrt((Xp-lon_obs(s)).^2+(Yp-lat_obs(s)).^2);
            [dif loc]=nanmin(dis(:));
            datas(kt,k,ke)=data(loc);
            udatas(kt,k,ke)=uwnd(loc);
            vdatas(kt,k,ke)=vwnd(loc);
          end  
    end % if plot_seriem==1

    if plot_map==1
  
      dp=mod(-90-dp,360);
      u=cosd(dp).*data; v=sind(dp).*data;
  
      u=double(u); v=double(v);
  
      % winds
      [dif iw]=nanmin(abs(timewnd-t));
      pair=double(ncread(filewnd,'P',[1 1 iw],[Inf Inf 1])')./100;
      %uwnd=double(ncread(filewnd,'U',[1 1 iw],[Inf Inf 1])');
      %vwnd=double(ncread(filewnd,'V',[1 1 iw],[Inf Inf 1])');
      if plot_seriem==1
        subplot(4,2,[1 3 5 7])
      end
 
      hold on
      %geolimits([-18 -13],[164.5 169.5])
      %geobasemap streets
      m_proj('lambert','long', [164.5 169.5],'lat',[-18 -13])
  
      set(gca,'fontsize',12,'fontweight','bold')
      m_pcolor(Xp,Yp,data); 
      %pcolor(Xp,Yp,dp); shading flat; colorbar; olormap(hsv)
      %pcolor(Xp,Yp,tp); shading flat; colorbar; 
      %shading flat; 
      caxis([0 20]); 
      colormap([jet])
      cb=colorbar;%('southoutside');
      set(get(cb,'ylabel'),'string',vname,'fontsize',12,'fontweight','bold');
      set(cb,'fontsize',12,'fontweight','bold');
  
      m_contour(Xp,Yp,Botlev,[0 20 50 100 200 300 500 1000],'color',[.3 .3 .3]);
      %contour(Xp,Yp,tp,[round(min(dp(:))):1:round(max(dp(:)))],'k')
  
      % TC track
      [dif itr]=nanmin(abs(time-t));
      m_plot(lon_track(1:itr),lat_track(1:itr),'r','linewidth',2)
  
      % pressure
      [cs,h]=m_contour(lonwnd,latwnd,pair,[990 1000 1010],'k');
      clabel(cs,h,'fontsize',10,'fontweight','bold','LabelSpacing',1000)  

      %m_gshhs_f('patch',[.7 .7 .7],'edgecolor','k','linewidth',.5);
      m_grid('box','fancy','fontname','helvetica','fontsize',14,'fontweight','bold');
  
      lon_vec=167; lat_vec=-15.25;
      
      %waves
      vspw=4;
      lonq=double(Xp); latq=double(Yp);
      %quiver(lonq(1:vspw:end,1:vspw:end),latq(1:vspw:end,1:vspw:end),u(1:vspw:end,1:vspw:end),v(1:vspw:end,1:vspw:end),'w')
      m_vec(80,lonq(1:vspw:end,1:vspw:end),latq(1:vspw:end,1:vspw:end),u(1:vspw:end,1:vspw:end),v(1:vspw:end,1:vspw:end),'w','shaftwidth',.25,'headwidth',3,'headangle',60,'headlength',3.5)
      %[hp ht]=m_vec(10,lon_vec+.1,lat_vec+.14,2,0,'k','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5,'key',[num2str(2.0),' m']);
      %set(ht,'fontsize',14,'fontweight','bold')
      %wind
  
      vspa=20;
      %lonq=lonwnd; latq=latwnd;
      lonq=double(Xp); latq=double(Yp);
      %quiver(lonq(1:vspa:end,1:vspa:end),latq(1:vspa:end,1:vspa:end),uwnd(1:vspa:end,1:vspa:end),vwnd(1:vspa:end,1:vspa:end),'k')
      m_vec(120,lonq(1:vspa:end,1:vspa:end),latq(1:vspa:end,1:vspa:end),uwnd(1:vspa:end,1:vspa:end),vwnd(1:vspa:end,1:vspa:end),'k','shaftwidth',.75,'headwidth',3,'headangle',60,'headlength',3.5)
      %axis equal
  
      if plot_obsm==1
        m_plot(lon_obs,lat_obs,'.k','markersize',24)
        m_plot(lon_obs,lat_obs,'.m','markersize',18)
        m_text(lon_obs+0.04,lat_obs+0.04,snames,'color','black','fontsize',12,'fontweight','bold')
        %m_text(lon_obs+0.04,lat_obs+0.04,snames,'color','magenta','fontsize',12,'fontweight','bold')
      end

      title({expt;[vname,' (max=',num2str(max(data(:)),'%.2f'),') and winds (max=',num2str(3.6*max(sqrt(uwnd(:).^2+vwnd(:).^2)),'%.2f'),'km/h) at ',datestr(t,'HH:MM dd/mm/yyyy')]})

			if plot_seriem==1
        for kee=1:length(runs)
          k=0; kk=0; 
          for s=stations
            k=k+1; kk=kk+2;
            subplot(length(stations),2,kk)
            set(gca,'fontsize',12,'fontweight','bold')
            hold on
            title(['Station: ',snames{s}])
            if plot_obs==1 & kee==1
              plot(time_obs(s),hs_obs(s),'.','color',[1 0 1],'markersize',30)
            end
            plot(time,datas(:,k,kee),colors{kee},'linewidth',2)
            datetick('x','dd/mm')
            ylabel(vname)
            ylim([0 15])
            grid('on')
          end 
        end 
        if plot_obs==1
          legend([{'Altimeter Obs.'},legh])
        else
          legend(legh)
        end
      end  
 
      path_dm=[path_fig,expt,'/'];
      system(['mkdir -p ',path_dm]);
  
      figname=[path_dm,'/map_',vname,'_',datestr(t,'yyyymmdd_HHMMSS'),'.png'];
      if save_fig==1
        display(['Saving: ',figname]);
        export_fig(gcf,figname,'-png','-r150' );
      elseif save_video==1
      else
        pause
      end
  
      if save_video==1 & plot_map==1;
        if t==time(1)
          videoName = [path_dm,'/',vname,'_',datestr(time(1),'yyyymmdd_HHMMSS'),'_',datestr(time(end),'yyyymmdd_HHMMSS'),''];
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
  
    end % if plot_map==1
  
  end % for t=time
  
  if save_video==1 & plot_map==1;
    close(vidObj);
  end


  if plot_serie==1
    k=0; 
    for s=stations
      k=k+1;
      subplot(length(stations),1,k)
      set(gca,'fontsize',12,'fontweight','bold')
      hold on
      title(['Station: ',snames{s}])

      if plot_obs==1 & ke==1
        plot(time_obs(s),hs_obs(s),'.','color',[1 0 1],'markersize',30)
      end
    
      plot(time,datas(:,k,ke),colors{ke},'linewidth',2)
      datetick('x','dd/mm')
      ylabel(vname)
      ylim([0 15])
      grid('on')

    end 
 
  end %if plot_serie==1
  if plot_seriem==1
    toc
    return
  end
end % for ex=1:length(runs) 


%close all

if plot_serie==1
  k=0;
  for s=stations
    k=k+1;
    subplot(length(stations),1,k)
    if plot_obs==1
      legend([{'Altimeter Obs.'},legh])
    else
      legend(legh)
    end
  end

  path_dm=[path_fig,expt,'/'];
  system(['mkdir -p ',path_dm]);
  figname=[path_dm,'/serie_',vname,'_',datestr(time(1),'yyyymmdd_HHMMSS'),'_',datestr(time(end),'yyyymmdd_HHMMSS'),'.png'];
  if save_fig==1
    display(['Saving: ',figname]);
    export_fig(gcf,figname,'-png','-r150' );
  end

end

save_wnd_sta=1;
if save_wnd_sta

  ke=0;
  for ex=runs 
    ke=ke+1;
    expt=expt_names(ex);
    expt=expt{1};
    matfile=['/scale_wlg_persistent/filesets/project/niwa03150/santanarc/swan/projects/VANUATU/',expt,'/',expt(6:end),'_winds_near_araki.mat'];
    display(['Saving: ',matfile])
    uwnd=squeeze(udatas(:,end,ke));
    vwnd=squeeze(vdatas(:,end,ke));
    lon_sta=lon_obs(end); lat_sta=lat_obs(end);
    save(matfile,'uwnd','vwnd','time','lon_sta','lat_sta')
  end
  
end


toc
display('End of script')

