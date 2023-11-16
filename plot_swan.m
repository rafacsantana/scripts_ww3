clear
%close all
clf('reset')


expt='luganville_3'; % 3';
%expt='lv_era5_tc_harold_corr'; % 3';
expt='lv_tcwindgen_tc_harold'; % 3';

if strncmp(expt,'lugan',5)
  time=datenum(2000,1,1,0,0,0):.5/24:datenum(2000,4,11,0,0,0);
elseif strncmp(expt,'lv',2)
  time=datenum(2020,4,6,0,0,0):.5/24:datenum(2020,4,11,0,0,0);
end


vname='Hsig'; % Dir, Dspr, RTpeak, Tm01

vspw=4;
vspa=40;

path_fig='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/swan/projects/VANUATU/figures/';

load(['/scale_wlg_persistent/filesets/project/niwa03150/santanarc/swan/projects/VANUATU/',expt,'/Van.mat']);

% winds
filewnd=('/scale_wlg_persistent/filesets/project/niwa03150/santanarc/swan/projects/VANUATU/data/TCwindgen/test_LV.nc');
lonwnd =double(ncread(filewnd,'x'));
latwnd =double(ncread(filewnd,'y'));
timewnd=double(ncread(filewnd,'time'))./(24*3600)+datenum(2000,01,01,0,0,0); %:1/48:datenum(2000,01,18,0,0,0);

scrsz=[1 1 1920 1080];
scrsz=get(0,'screensize');
%figure('position',scrsz,'color',[1 1 1],'visible','on');  hold on;

for t=time
  
  %hs=Hsig_20000101_033000
  display(['Loading: ',vname,'_',datestr(t,'yyyymmdd_HHMMSS')]);
  data=double(eval(['Hsig_',datestr(t,'yyyymmdd_HHMMSS')]));
  dp=double(eval(['Dir_',datestr(t,'yyyymmdd_HHMMSS')]));
  tp=double(eval(['RTpeak_',datestr(t,'yyyymmdd_HHMMSS')]));
  uwnd=double(eval(['Windv_x_',datestr(t,'yyyymmdd_HHMMSS')]));
  vwnd=double(eval(['Windv_y_',datestr(t,'yyyymmdd_HHMMSS')]));

  dp=mod(-90-dp,360);
  u=cosd(dp).*data; v=sind(dp).*data;

  % winds
  iw=find(timewnd==t);
  %pair=double(ncread(filewnd,'P',[1 1 iw],[Inf Inf 1])')./100;
  %uwnd=double(ncread(filewnd,'U',[1 1 iw],[Inf Inf 1])');
  %vwnd=double(ncread(filewnd,'V',[1 1 iw],[Inf Inf 1])');

  hold on
  set(gca,'fontsize',12,'fontweight','bold')
  pcolor(Xp,Yp,data); shading flat; colorbar; caxis([0 20]); colormap([jet])
  %pcolor(Xp,Yp,dp); shading flat; colorbar; caxis([0 360]); colormap(hsv)
  %pcolor(Xp,Yp,tp); shading flat; colorbar; 
  contour(Xp,Yp,Botlev,[50 100 200 300 500 1000],'color',[.5 0 0]); shading flat; colorbar; 
  
  %contour(Xp,Yp,tp,[round(min(dp(:))):1:round(max(dp(:)))],'k')

  % pressure
  %[cs,h]=contour(lonwnd,latwnd,pair,'k'); % [round(min(dp(:))):1:round(max(dp(:)))],'k');
  %clabel(cs,h,'fontsize',10,'fontweight','bold','LabelSpacing',1000)  

  %waves
  lonq=Xp; latq=Yp;
  quiver(lonq(1:vspw:end,1:vspw:end),latq(1:vspw:end,1:vspw:end),u(1:vspw:end,1:vspw:end),v(1:vspw:end,1:vspw:end),'w')

  %wind
  vspa=40;
  lonq=lonwnd; latq=latwnd;
  lonq=Xp; latq=Yp;
  quiver(lonq(1:vspa:end,1:vspa:end),latq(1:vspa:end,1:vspa:end),uwnd(1:vspa:end,1:vspa:end),vwnd(1:vspa:end,1:vspa:end),'k')
  axis equal

  title([vname,' and winds (max. = ',num2str(3.6*max(sqrt(uwnd(:).^2+vwnd(:).^2)),'%.2f'),' km/h) at ',datestr(t,'HH:MM dd/mm/yyyy')])

  xlim([164.5 169.5])
  ylim([-18 -13])

  save_fig=0;
  if save_fig==1
    path_dm=[path_fig,expt,'/'];
    system(['mkdir -p ',path_dm]);
    figname=[path_dm,'/map_',vname,'_',datestr(t,'yyyymmdd_HHMMSS'),'.png'];
    display(['Saving: ',figname]);
    export_fig(gcf,figname,'-png','-r150' );
  else
    pause
  end
  clf('reset')
  set(gcf,'color',[1 1 1])

end


