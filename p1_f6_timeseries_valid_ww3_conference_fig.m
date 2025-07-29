clear
close_all=1;
if close_all==1; close all; end
warning off
tic
  
% Daily time
time_lima=datenum(1983,1,1,0,0,0):1:datenum(1983,3,31,0,0,0);    % hindcast
%time_lima=datenum(1983,1,1,0,0,0):1:datenum(2002,12,31,0,0,0);    % hindcast
time_lima=datenum(2014,1,1,0,0,0):1:datenum(2022,12,31,0,0,0);    % hindcast
%time_lima=datenum(1996,1,1,0,0,0):1:datenum(2002,12,31,0,0,0);    % hindcast
time_lima=datenum(2003,1,1,0,0,0):1:datenum(2005,12,31,0,0,0); %

%time_lima=datenum(2020,12,1,0,0,0):1:datenum(2020,12,15,0,0,0);  % 1-year experiments
time_lima=datenum(2021,1,1,0,0,0):1:datenum(2021,12,31,0,0,0);  % 1-year experiments

%time_lima=datenum(2014,7,9,0,0,0):1:datenum(2014,9,3,0,0,0);     % initial date NZLAM hindcast
%time_lima=datenum(2018,1,6,0,0,0):1:datenum(2023,11,15,0,0,0);   % Graham Harrington ECAN request on 18/09/2023

stations=2;

% 1 {'PeakPerid'}   2 {'PeakDirection'}   3 {'Directional spread'}   4 {'Tm01'}   5 {'Tm02'}   6 {'Hs'}   7 {'Qp'}  8 Tm
% 9 wlv 10 ucurr 11 vcurr
if stations==1 % Banks_Peninsula

dcol=[6,4,11]; % Hs, v
%dcol=[9,2,11]; % 
time_lima=datenum(2021,10,20,0,0,0):1:datenum(2021,10,23,0,0,0); % Banks Peninsula head tide and direction
%time_lima=datenum(2021,10,21,0,0,0):1:datenum(2021,10,22,0,0,0); % Banks Peninsula head tide and direction

elseif stations==2 % Baring_Head

dcol=[6,2,11]; % Hs, wlv
%dcol=[9,2,11]; % wlv, direction and v
time_lima=datenum(2021,1,10,0,0,0):1:datenum(2021,1,15,0,0,0);    % Baring head wave direction, time series
%time_lima=datenum(2021,1,13,0,0,0):1:datenum(2021,1,13,0,0,0);   % Baring head wave direction and currents, maps
%time_lima=datenum(2021,5,20,0,0,0):1:datenum(2021,5,27,0,0,0);   % Baring Head tide and swell

end

%time_lima=datenum(2029,12,1,0,0,0):1:datenum(2030,1,1,0,0,0); % climate projections

files={'Banks_Peninsula','Baring_Head','wairewa_lake_forsyth','Pelorus_Sound','Taumutu'};% 'SteepHead_SP_FullRecord_QC';

runs=[10];
%runs=[3,13];
runs=[4,5];

expt_names={'GLOBALWAVE'          ,'NZWAVE'             ,'NZWAVE-ST6'         ,'NZWAVE-HR-NOTIDES','NZWAVE-HR',... 5
            'GLOBALWAVE-GFDL-CCAM','NZWAVE-ST4-GLOBALUM','NZWAVE-ST6-GLOBALUM','NZWAVE-GFDL-CCAM' ,'NZWAVE-ERA5',...10
            'GLOBALWAVE-ERA5'     ,'NZWAVE-ERA5-2021'   ,'NZWAVE-ST6-NZLAM'};            

plot_series=1; % stations=0, plot_series=1, plot_coastm=1, plot_obs==0 to save the coastal points
plot_rose  =0;
plot_disp  =0;
plot_eva   =0;

plot_coastm=0; % statistical maps for coastal points
plot_map   =0; % daily station maps
plot_downs =0; % daily downscaled maps
plot_mmap  =0; % mean model maps
plot_mserie=0; % mean model series
plot_amap  =0; % altimeter maps
plot_aserie=0; % altimeter series

plot_obs   =1;
proc_obs   =0;
ck_mod_mat =1;
save_csv   =0;
save_fig   =0;
save_video =0;
save_gif   =0;

set_xlim   =1;
obs_first  =1;
add_stats  =0;
portrait   =1;

correct_model=0; % -1=scatter_kde, 0=binscatter, 1=nothing, >1 correct it
time_limab=time_lima;

% series, rose, dispersion, eva


% maps
pcname ='hs'; % m_pcolor: wlv, hs, curr
vecname='hs'; % m_vec: hs, curr

% mmaps
mcol=[1,3]; % selected variables from mvalm
mvalm={'hs_m_mean','tp_m_mean','we_m_mean'}; 
mvals={'hs_mean','tp_mean','we_mean'}; % 1=hs, 2=tp, 3=we

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

end

%fname=[pname,ww3pre,'_',ftime,'-utc_',gnames{ig},'.nc'];
fname='/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/NZWAVE/2021/01/01/00/ww3g_2021010100-utc_nzwave+nzlam.nc';
depth_mod=double(ncread(fname,'depth',[1 1 1],[Inf Inf 1]));
lon_mod=double(ncread(fname,'lon'));
lat_mod=double(ncread(fname,'lat'));

fname='/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/NZWAVE-HR/2021/01/01/00/ww3g_2021010100-utc_nzwave_hr+nzcsm.nc';
depth_hr=double(ncread(fname,'depth',[1 1 1],[Inf Inf 1]));
lon_hr=double(ncread(fname,'lon'));
lat_hr=double(ncread(fname,'lat'));

if plot_downs %|| plot_coastm
%if plot_map || plot_mmap || plot_downs || plot_coastm

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

end


if plot_coastm || stations==0

  % finding the closest points to the coast
  np=1; % number of the closest points 
  xl=[390,605]; yl=[180,520];
  xl=[390,700]; yl=[180,520]; % including the Chatam Islands
  dep=depth_mod(xl(1):xl(2),yl(1):yl(2));

  % finding y points in the x direction
  k=0;
  for i=1:size(dep,1)
    [ngap,igap,fgap,sizegap] = findgap(dep(i,:)');
    if ngap>0
      for j=1:ngap; k=k+1;
        points=[igap(j)-np:igap(j)-1,fgap(j)+1:fgap(j)+np];
        yp(k).p=points;
        xp(k).p(1:length(points))=i;
        end; end; end
  % finding x points in the y direction
  k=0;
  for i=1:size(dep,2)
    [ngap,igap,fgap,sizegap] = findgap(dep(:,i));
    if ngap>0
      for j=1:ngap; k=k+1;
        points=[igap(j)-np:igap(j)-1,fgap(j)+1:fgap(j)+np];
        xpp(k).p=points;
        ypp(k).p(1:length(points))=i;
        end; end; end
  % concatenating x and y points
  xpoints=[]; ypoints=[];
  for i=1:length(xp)
    xpoints=[xpoints,xp(i).p]; 
		ypoints=[ypoints,yp(i).p];
  end
  for i=1:length(xpp)
    xpoints=[xpoints,xpp(i).p]; 
		ypoints=[ypoints,ypp(i).p];
  end

  % cleaning repeated pair of points
	points=[xpoints',ypoints'];
  points=unique(points,'rows');
  points(:,1)=points(:,1)+xl(1)-1;
  points(:,2)=points(:,2)+yl(1)-1;

  % selecting from spectral points
  spec_points=0;
  if spec_points
    fname='/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/NZWAVE-ERA5/2021/01/01/00/ww3p_2021010100-utc_nzwave_era5.nc';
    lon_spe=double(ncread(fname,'lon'));
    lat_spe=double(ncread(fname,'lat'));
    depth_spe=double(ncread(fname,'depth'));

    a=find(lon_spe(:,1)<=178.7 & lon_spe(:,1)>=166.2 & lat_spe(:,1)<=-34 & lat_spe(:,1)>=-48);
    b=find(depth_spe(a,1)<=1000);
    c=a(b);
    a=find(lon_spe(c,1)<=170.6 & lon_spe(c,1)>=168 & lat_spe(c,1)<=-37.5 & lat_spe(c,1)>=-42.5);
    c(a)=nan; return
    c(isnan(c))=[];
    figure; hold on; plot(lon_spe(c,1),lat_spe(c,1),'.'); xlim([165 180]); ylim([-50 -30])
	end


end

lon_obss=[]; lat_obss=[];

if plot_coastm==1 & plot_series==0
	stations=1:length(points);
end


% biggest loop on stations
for fobs=stations
  
  % OBS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if plot_coastm==1 & plot_series==0
    file=['i',num2str(points(fobs,1)),'_j',num2str(points(fobs,2)),'_grid_point_'];
  elseif plot_coastm==1 & plot_series==1
    file=['bla'];
    %file=['i',num2str(points(fobs,1)),'_j',num2str(points(fobs,2)),'_grid_point_'];
  else
    file=files{fobs};
    file=replace(file,' ','_');
    file_obs=[path_obs,file];
    [obs,time_obs,lon_obs,lat_obs]=proc_wave_obs(file,proc_obs); %,plot_obs);   
  end
  
  if plot_series
  
   
    if portrait
      scrsz=[2    42   958   953];
    else
      scrsz=[1 41 1920 962];
      %scrsz=[1 1 1366 768];
      %scrsz=get(0,'screensize');
    end


 
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
        %tic
        if plot_coastm==1 & plot_series==0 % plot_series=0 needed to save the coastal points
          file=['i',num2str(points(fobs,1)),'_j',num2str(points(fobs,2)),'_grid_point_'];
        elseif plot_coastm==1 & plot_series==1
          file=['bla'];
          %file=['i',num2str(points(fobs,1)),'_j',num2str(points(fobs,2)),'_grid_point_'];
        elseif ke==1 & t==time_lima(1)
          %file=files{fobs};
          %file=replace(file,' ','_');
          %file_obs=[path_obs,file];
          %[obs,time_obs,lon_obs,lat_obs]=proc_wave_obs(file,proc_obs); %,plot_obs);   
        end
        filename=[path_dm,file,'_',datestr(t,'YYYYmmDDHH'),'.mat'];
    
        if ck_mod_mat==1 && exist(filename)==2
          display(['Loading: ',filename])
          load(filename)%,'time_mod','model')
          if exist('ucur','var')~=1
            wlv=nan(size(hs)); ucur=nan(size(hs)); vcur=nan(size(hs));
            we=nan(size(hs)); ucur=nan(size(hs)); vcur=nan(size(hs));
          end 
 
        elseif fobs==0 %plot_coastm==1

          ptime=datestr(t,'YYYY/mm/DD/HH/');
          ftime=datestr(t,'YYYYmmDDHH');
          pname=[path_expt,ptime];
          fname=[pname,ww3pre,'_',ftime,'-utc_',gnames{ig},'.nc'];
          display(['Loading: ',fname]);

          if checkin==0; % t==time_lima(1)
            checkin=1; 
            lon_mod=double(ncread(fname,'lon'));
            lat_mod=double(ncread(fname,'lat'));
          end

         	stations=1:length(points);
          for fobss=stations
            %tic; 

            file=['i',num2str(points(fobss,1)),'_j',num2str(points(fobss,2)),'_grid_point_'];
            filename=[path_dm,file,'_',datestr(t,'YYYYmmDDHH'),'.mat'];
            display(['Checking: ',filename])

            if ck_mod_mat==0 || exist(filename)~=2

              lon_obs=lon_mod(points(fobss,1));
						  lat_obs=lat_mod(points(fobss,2));
              [dif ilon]=nanmin(abs(lon_mod-lon_obs));
              [dif ilat]=nanmin(abs(lat_mod-lat_obs));
              dis=sqrt((lon_mod(ilon)-lon_obs).^2 + (lat_mod(ilat)-lat_obs).^2)*111;
              %display(['Distance between obs ',file,' and grid point is: ',num2str(dis),' km'])
          
              time=squeeze(double(ncread(fname,'time')))./24+t; time=time(1:ltime-1);
              if strncmp(ww3pre,'ww3p',4)
              elseif strncmp(ww3pre,'ww3g',4)
                hs=squeeze(double(ncread(fname,'hsig',          [ilon ilat 1],[1 1 ltime]))); hs=hs(1:end-1);
                tp=squeeze(double(ncread(fname,'tpeak',         [ilon ilat 1],[1 1 ltime]))); tp=tp(1:end-1);
                %tm01=squeeze(double(ncread(fname,'tmean01',     [ilon ilat 1],[1 1 ltime]))); tm01=tm01(1:end-1);
                %tm02=squeeze(double(ncread(fname,'tmean02',     [ilon ilat 1],[1 1 ltime]))); tm02=tm02(1:end-1);
                %tm=squeeze(double(ncread(fname,'tmean',         [ilon ilat 1],[1 1 ltime]))); tm=tm(1:end-1);
                pd=squeeze(double(ncread(fname,'peak_direction',[ilon ilat 1],[1 1 ltime]))); pd=pd(1:end-1);
                %ds=squeeze(double(ncread(fname,'directional_spread',[ilon ilat 1],[1 1 ltime]))); ds=ds(1:end-1);
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
              %wave energy (we) = (1/16)*rho*g*(hs.^2); % https://www.coastalwiki.org/wiki/Shallow-water_wave_theory#Introduction
              we=(1/16)*1025*9.8*(hs.^2);    

              display(['Saving:   ',filename])
              save(filename,'time','tp','hs','pd')%,'ds','tm01','tm02','pd','tm','wlv','ucur','vcur','we')
            end %if ck_mod_mat==0 || exist(filename)~=2
            %toc
          end % for fobss=stations
      
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
            % coastal points
            if plot_coastm==1
              lon_obs=lon_mod(points(fobs,1));
							lat_obs=lat_mod(points(fobs,2));
            end
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
          %wave energy (we) = (1/16)*rho*g*(hs.^2); % https://www.coastalwiki.org/wiki/Shallow-water_wave_theory#Introduction
          we=(1/16)*1025*9.8*(hs.^2);    

          display(['Saving: ',filename])
          save(filename,'time','tp','pd','ds','tm01','tm02','hs','tm','wlv','ucur','vcur','we')
  
        end % if ck_mod_mat==1 && exist(filename)==2
        
        if plot_coastm==0 
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

        end
 
      end % for t=time_lima

      if plot_coastm==1
        display(['Finished processing coastal stations from: ',datestr(time_lima(1)),' to ',datestr(time_lima(end))]) 
        toc
        return
      end

    
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
      tnames={'Peak Period',     'Peak Direction',    'Directional spread',   'Tm01','Mean Period (Tm02)','Significant Wave Height (Hs)','Qp','Tm','Water Level','U current','V current','Wave Energy (J/m2)'};
      lnames={'Peak Period9 (s)','Peak Direction (^o)',    'Directional spread',   'Tm01 (s)','Tm02 (s)','Hs (m)','Qp','Tm (s)','Water Level (m)','U current (m/s)','V current (m/s)','Wave Energy (J/m2)'};
      for i=1:length(dcol)
        
        subplot(length(dcol),1,i);
        set(gca,'fontsize',12,'fontweight','bold')
        hold on

        lett={'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)','(m)','(n)','(o)','(p)','(q)','(r)','(s)','(t)','(u)','(v)','(w)','(x)','(y)','(z)'};
        if dcol(i)~=11
          if ke==1 & plot_obs==1
  	        if obs_first==1;
				      if dcol(i)<=size(obs,2) 
                plot(time_obs,obs(:,dcol(i))','.','markersize',12,'color',[0 .7 0],'linewidth',2)
	   		  		end
            end
            for eix=1:length(expts); plot(nan,nan,'.','markersize',12,'color',colors{eix},'linewidth',2); end
            plot(nan,nan,'.','markersize',12,'color',[0 .7 0],'linewidth',2)
          end
          plot(time_mod,model.data(:,dcol(i))','.','markersize',12,'color',colors{ke},'linewidth',2)
        end
 
        if dcol(i)==11 & strcmp(expt,'NZWAVE-HR')
          plot(time_mod,model.data(:,dcol(i))','.','markersize',12,'color',colors{ke},'linewidth',2)
        end

        if dcol(i)~=11
          if obs_first==0; 
          	if ke==length(expts) & plot_obs==1
	    	  		if dcol(i)<=size(obs,2) 
              	plot(time_obs,obs(:,dcol(i))','.','markersize',12,'color',[0 .7 0],'linewidth',2)
	   		  	  end
          	end
          end
        end

        if set_xlim==1;
          xlim([time_lima(1) time_lima(end)])
        end
      
        %if i==3
          datetick('x','dd/mmm/yyyy','keeplimits')
        %end
        title([lett{i},' ',replace(file,'_',' '),' ',tnames{dcol(i)}])
        ylabel([lnames{dcol(i)}])
        grid('on')
        if plot_obs==1 
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
				    ylim([120 220])
          end
        else
          legh=[legh,{[expt]}];%,' rmse=',num2str(mod_rmse,'%.2f'),' corr=',num2str(mod_corr,'%.2f'),', bias=',num2str(mod_bias,'%.2f')]}];
        end

        % read and plot TELEMAC output
        if strcmp(lett{i},'(c)') & strcmp(expt,'NZWAVE-HR') & stations==2
          path_t='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/data/NZTIDE-TELEMAC/';
          file_t='modelled_data_UMGlobal';
          file_t=[path_t,file_t];
          if exist([replace(file_t,' ','_'),'.mat'])==0
            te=importdata([file_t,'.csv']);
            txt=te.textdata;
            time_t=nan(length(te.textdata)-1,1);
            for j=2:length(te.textdata)
                tx=txt{j,1};
                time_t(j-1)=datenum([tx(1:19)],'YYYY-mm-ddTHH:MM:SS');
            end
            telemac=te.data;
            save([replace(file_t,' ','_'),'.mat'],'time_t','telemac')
          else
            load([replace(file_t,' ','_'),'.mat']) % ,'time_obs','wlv')
          end
          %plot(time_t,telemac(:,end)','.','markersize',12,'color',[0 0 0],'linewidth',2)
        end
        
        % plotting a red rectangle around 21st and 22nd midnight of October, if fbos==1, i.e., Banks Peninsula
				% getting the y limits and plotting the rectangle based on that
				had=3;
        if fobs==1 & ke==length(runs)
					amin=min(get(gca,'ylim')); amax=max(get(gca,'ylim'));
					plot([datenum('2021-10-22 00:00:00')-had/24 datenum('2021-10-22 00:00:00')-had/24],[amin amax],'r','linewidth',2)
					plot([datenum('2021-10-22 00:00:00')+had/24 datenum('2021-10-22 00:00:00')+had/24],[amin amax],'r','linewidth',2)
					plot([datenum('2021-10-22 00:00:00')-had/24 datenum('2021-10-22 00:00:00')+had/24],[amin amin],'r','linewidth',2)
					plot([datenum('2021-10-22 00:00:00')-had/24 datenum('2021-10-22 00:00:00')+had/24],[amax amax],'r','linewidth',2)

					plot([datenum('2021-10-21 12:00:00')-had/24 datenum('2021-10-21 12:00:00')-had/24],[amin amax],'r','linewidth',2)
					plot([datenum('2021-10-21 12:00:00')+had/24 datenum('2021-10-21 12:00:00')+had/24],[amin amax],'r','linewidth',2)
					plot([datenum('2021-10-21 12:00:00')-had/24 datenum('2021-10-21 12:00:00')+had/24],[amax amax],'r','linewidth',2)
          plot([datenum('2021-10-21 12:00:00')-had/24 datenum('2021-10-21 12:00:00')+had/24],[amin amin],'r','linewidth',2)
        end
        % plotting a red rectangle around 1 am and 12pm on the 13th of January, if fbos==1, i.e., Banks Peninsula
				had=2;
        if fobs==2 & ke==length(runs)
					amin=min(get(gca,'ylim')); amax=max(get(gca,'ylim'));
					plot([datenum('2021-01-13 02:00:00')-had/24 datenum('2021-01-13 02:00:00')-had/24],[amin amax],'r','linewidth',2)
					plot([datenum('2021-01-13 02:00:00')+had/24 datenum('2021-01-13 02:00:00')+had/24],[amin amax],'r','linewidth',2)
					plot([datenum('2021-01-13 02:00:00')-had/24 datenum('2021-01-13 02:00:00')+had/24],[amin amin],'r','linewidth',2)
					plot([datenum('2021-01-13 02:00:00')-had/24 datenum('2021-01-13 02:00:00')+had/24],[amax amax],'r','linewidth',2)
					plot([datenum('2021-01-13 13:00:00')-had/24 datenum('2021-01-13 13:00:00')-had/24],[amin amax],'r','linewidth',2)
					plot([datenum('2021-01-13 13:00:00')+had/24 datenum('2021-01-13 13:00:00')+had/24],[amin amax],'r','linewidth',2)
					plot([datenum('2021-01-13 13:00:00')-had/24 datenum('2021-01-13 13:00:00')+had/24],[amin amin],'r','linewidth',2)
					plot([datenum('2021-01-13 13:00:00')-had/24 datenum('2021-01-13 13:00:00')+had/24],[amax amax],'r','linewidth',2)
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
      if stations==2
        subplot(length(dcol),1,3);
        legend('RiCOM','TELEMAC')
      else
        subplot(length(dcol),1,3);
        legend('RiCOM')
      end
    
    if save_fig==1
      %export_fig(gcf,'steephead_period','-png','-r150');
    end

  end % if plot_series

  path_dm=strcat(path_fig,'paper_1/');
  system(['mkdir -p ',path_dm]);
  figname=[path_dm,file,'_tides_',replace(tnames{dcol(i)},' ','_'),'_',datestr(time_lima(1),'YYYY_mm_DD_HH_MM'),'_',datestr(time_lima(end),'YYYY_mm_DD_HH_MM'),'_no_Telemac.png'];
  display(['Plotting: ',figname]);
  if save_fig==1

return

    display(['Saving: ',figname]);
    export_fig(gcf,figname,'-png','-r150');
    %close
    %clf('reset')
  end




end % fobs=fstations





% compute maps of annual average wave height and plot them

toc
