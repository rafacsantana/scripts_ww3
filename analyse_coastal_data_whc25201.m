%read tidal stations csv files
clear
close all

% sea level
ll=3;
stations=[ll]; % always use 2 for Nelson Port at Main Wharf_dataset1
lname={'hatea','port','harbour'};
file={'level_hatea_at_town_basin_07_11_24 09_19_28','level_Port_Whangarei_1987-2007','level_whangarei_harbour_at_marsden_point_07_11_24 09_46_58'};

rr=[1]; % river
river_name={'hatea','raumanga','waiarohia'};
river_stations=[rr];
file_river={'flow_hatea_at_whareora_rd_07_11_24 14_31_58','flow_raumanga_at_bernard_st_07_11_24 14_29_44','flow_waiarohia_at_lovers_lane_07_11_24 14_35_58'};
%rain_river={'hatea','otaika','waiarohia'};

ra=3; % always 3 (town centre) for Whangarei
rain_stations=[ra]; % always 3 (town centre) for Whangarei
file_rainfall={'rainfall_hatea_at_glenbervie_forest_hq_07_11_24 11_19_21','rainfall_taika_at_cemetery_road_07_11_24 11_17_58','rainfall_waiarohia_at_nrc_water_st_07_11_24 11_14_29'};

ocean_source=2;
%file_ocean={'EcoConnect_waves_hydro_Blenheim_2018010600_2024082000.mat','Hindcast_waves_Blenheim_1984010100_2023123100.mat'};

% file_obs=[path,'JPdata_dailymax_Riverflow_skewsurge.csv']; there is code to analyse this file at the end of this script

%% water level data
proc_wlv   =1; % process water level data
read_csv   =0;
apply_qaqc =2; % 0 load existing non qaqced .mat, 1 apply qaqc, 2 load qaqced .mat 
apply_harm =1; % 1 to apply harmonic analysis
read_harm  =1;
plot_data  =1; % 1
manual_qaqc=0; % 1 to run ginput
calc_tides =0; 
plot_tides =0;
use_anomaly=0; % tides analysis is done with MSL=0m, if 1, aep includes the mean sea level

skew_aep   =0;
skew_am=1;
skew_pot=1;
calc_aep   =1; % gev
calc_pot   =1; % gpd 

%% river flow data
proc_river=1; % process/plot river flow data
read_river=0; % 1=csv, 0=matlab
plot_river=0; % script stops here

%% rainfall data
proc_rainfall=1; % process/plot rainfall data
read_rainfall=0; % 1=csv, 0=matlab
plot_rainfall=0; % script stops here
daily_rain   =0; % 1=plot daily rainfall

%% wave
proc_wave=0; % process wave data
plot_wave=0; % plot wave data

crop_river_level=0;
plot_daily_max  =0;
write_jp_file   =0;

plot_all=0; % subplots of different variables

predict_river_wlv=0;
          
detrend_data=1;
save_fig=0;

nzvd_off=[1751,1751,1765]; % NZVD offset for each station

path_proj='\\niwa.local\projects\hamilton\WHC25201\Working\JPanalysis\';
path=[path_proj,'RawData\'];
path_flow=path; 
path_fig=[path_proj,'matlab_figures\'];

path_jp=[path_proj,'Input\'];

scrsz=[1 41 1920 962];

if proc_wlv==1

  for fobs=stations % 1:length(file)
  
      disp(['Processing: ',file{fobs}]);
  
      file_obs=[path,file{fobs}];
      if read_csv==1
        ob=importdata([file_obs,'.csv']);
        txt=ob.textdata;
        time_obs=nan(length(ob.textdata)-1,1);

	if strcmp(file{fobs},'level_Port_Whangarei_1987-2007')
          for j=1:length(ob.textdata)
              time_obs(j)=datenum([txt{j}],'dd-mmm-yyyy HH:MM:SS,');
          end
          wlv=ob.data.*1000;
          % processing observations
          % interpolating to a 5 minute interval
          time_int=[time_obs(1):1/24/(60/5):time_obs(end)]'; % 5 minutes
          %wlv=wlv-nanmean(wlv(:));
          wlv_int=interp1(time_obs,wlv,time_int,'spline');
          wlv_int=qaqc_time_gap(time_obs,time_int,wlv_int);

	else
          for j=2:length(ob.textdata)
              time_obs(j-1)=datenum([txt{j,2}],'yyyy-mm-dd HH:MM:SS');
          end
          wlv=ob.data;
          % processing observations
          % interpolating to a 5 minute interval
          time_int=[time_obs(1):1/24/(60/5):time_obs(end)]'; % 5 minutes
          %wlv=wlv-nanmean(wlv(:));
          wlv_int=interp1(time_obs,wlv,time_int,'spline');
          wlv_int=qaqc_time_gap(time_obs,time_int,wlv_int);
	end
  
        save([replace(file_obs,' ','_'),'.mat'],'time_obs','wlv','time_int','wlv_int') 
        return
  
      else
        %load([replace(file_obs,' ','_'),'_int.mat']) % data with gaps as in the original dataset (it takes a while to do that) 
        %wlv_int=wlv; 
        load([replace(file_obs,' ','_'),'.mat']) % ,'time_obs','wlv') 
      end	

      if apply_qaqc==2
        if exist([replace(file_obs,' ','_'),'_qaqced.mat'])==2
          disp(['loading: ',replace(file_obs,' ','_'),'_qaqced.mat'])
          load([replace(file_obs,' ','_'),'_qaqced.mat'])
	end
      end

      % removing NZVD offset
      wlv_int=wlv_int-nzvd_off(fobs);
      wlv=wlv-nzvd_off(fobs);

      if strcmp(file{fobs},'level_hatea_at_town_basin_07_11_24 09_19_28') 
  
        %applying qaqc
        if apply_qaqc==1
          data=wlv_int;
  
          gfilename='level_hatea_ginput_qaqc_outputs.txt';
  	  if exist(gfilename,'file')==2;
  	    display(['REMEMBER TO DELETE ',gfilename,', if you re using this file from a previous study.']);	
            data=qaqc_read_ginput_data(time_int,data,gfilename); close;
  	  end
  
  	% removing large gaps
          x1=728424.545195489; x2=733703.496539174; [dif loc1]=nanmin(abs(time_int-x1)); [dif loc2]=nanmin(abs(time_int-x2)); 
	  data(1:loc2)=nan;
	  %time_int(1:loc2)=[];
  
          mean_data=nanmean(wlv(:));
          mean_data=nanmean(data(:));
  
          % harmonic analysis
          if read_harm==0
            loc1=1; 
  	    [dif loc1]=nanmin(abs(time_int-733703.496539174)); 
            [tidestruc,data_oharm]=t_tide(data(loc1:end),'interval',(time_int(2)-time_int(1))*24,'start time',time_int(loc1),'latitude',-35.766627,'synthesis',1); %,'output','none'););
  	    data_oharm=t_predic(time_int,tidestruc,'synthesis',1);
            data_oharm=data_oharm+nanmean(data);
            save([replace(file_obs,' ','_'),'_harm.mat'],'data_oharm') 
          else 
            disp(['loading: ',replace(file_obs,' ','_'),'_harm.mat']) % ,'time_obs','wlv') 
            load([replace(file_obs,' ','_'),'_harm.mat']) % ,'time_obs','wlv') 
          end
          resid=data-data_oharm;
          resido=resid;
  
          data(isnan(resid))=nan; 	
          wlv_int=data;
          save([replace(file_obs,' ','_'),'_qaqced.mat'],'time_int','wlv_int')
  	
        elseif apply_qaqc==0
          load([replace(file_obs,' ','_'),'.mat']) % ,'time_obs','wlv')
  
        elseif apply_qaqc==2
          load([replace(file_obs,' ','_'),'_qaqced.mat']) % ,'time_obs','wlv')
          load([replace(file_obs,' ','_'),'_harm.mat']) % ,'time_obs','wlv') 
          resid=wlv_int-data_oharm;
  	
        end

      elseif strcmp(file{fobs},'level_Port_Whangarei_1987-2007') 
  
        %applying qaqc
        if apply_qaqc==1
          data=wlv_int;
  
          gfilename='level_port_whangarei_ginput_qaqc_outputs.txt';
  	  if exist(gfilename,'file')==2;
  	    display(['REMEMBER TO DELETE ',gfilename,', if you re using this file from a previous study.']);	
            data=qaqc_read_ginput_data(time_int,data,gfilename); close;
  	  end
  
          % finding gaps in the original data and inserting nans on the interpolated data
          %[dif iloc]=find(diff(time_obs)>1/24/(60/10)); % 10 minutes
          %for i=1:length(iloc)
          %  [dif iini]=nanmin(abs((time_obs(iloc(i))-time_int)));
          %  [dif ifin]=nanmin(abs((time_obs(iloc(i)+1)-time_int)));
          %  data(iini:ifin)=nan;
          %end
   
          %manual clean up
      	  %bad begining
          %[dif loc]=nanmin(abs(time_int-727858.80943138)); data(1:loc)=nan;
          %[dif loc]=nanmin(abs(time_int-datenum(2005,1,1))); data(1:loc)=nan;
  	  % removing large gaps
          %x1=728424.545195489; x2=733703.496539174; [dif loc1]=nanmin(abs(time_int-x1)); [dif loc2]=nanmin(abs(time_int-x2)); data(1:loc2)=nan;
          %x1=737094.839857142; x2=737111.736473575; [dif loc1]=nanmin(abs(time_int-x1)); [dif loc2]=nanmin(abs(time_int-x2)); data(loc1:loc2)=nan;
          %x1=737867.751552838; x2=737869.350375492; [dif loc1]=nanmin(abs(time_int-x1)); [dif loc2]=nanmin(abs(time_int-x2)); data(loc1:loc2)=nan;
  
          % detrending
          mean_data=nanmean(wlv(:));
          %wlv=detrend(wlv,1,'omitnan')+mean_data;
  
          % detrending
          mean_data=nanmean(data(:));
          %data=detrend(data,1,'omitnan')+mean_data;
          %wlv=data; time_obs=time_int;
  
          %% removing extremes
          %data=qaqc_extremes(data,-600,4000);
          %%% removing large local differences
          %%%if strcmp(file{fobs},'Nelson Port at Fairway Beacon') % proc_obs==1;     
          %data=qaqc_fwd_diff(data,300,0); % 300
          %data=qaqc_back_diff(data,300,0); % 300
          %data=qaqc_remove_lin_interp(data,1.000002);
          data=qaqc_time_gap(time_obs,time_int,data);
          %%end
  
          % harmonic analysis
	  if apply_harm==1
            if read_harm==0
              loc1=1; 
  	      [dif loc1]=nanmin(abs(time_int-730355.378757781)); 
              [tidestruc,data_oharm]=t_tide(data(loc1:end),'interval',(time_int(2)-time_int(1))*24,'start time',time_int(loc1),'latitude',-35.766627,'synthesis',1); %,'output','none'););
  	      data_oharm=t_predic(time_int,tidestruc,'synthesis',1);
              data_oharm=data_oharm+nanmean(data);
              save([replace(file_obs,' ','_'),'_harm.mat'],'data_oharm') 
            else 
              disp(['loading: ',replace(file_obs,' ','_'),'_harm.mat']) % ,'time_obs','wlv') 
              load([replace(file_obs,' ','_'),'_harm.mat']) % ,'time_obs','wlv') 
            end
            resid=data-data_oharm;
            resido=resid;
          end 
          % removing obs with resid surge above 1.5m and below -1.5m
          resid=qaqc_extremes(resid,-400,1500); 
          resid=qaqc_fwd_diff(resid,60,0); resid=qaqc_fwd_diff(resid,60,0); resid=qaqc_fwd_diff(resid,60,0); resid=qaqc_fwd_diff(resid,60,0); resid=qaqc_back_diff(resid,60,0); % 300
          resid=qaqc_fwd_diff(resid,60,0); resid=qaqc_fwd_diff(resid,60,0); resid=qaqc_fwd_diff(resid,60,0); resid=qaqc_fwd_diff(resid,60,0); resid=qaqc_back_diff(resid,60,0); % 300
          resid=qaqc_fwd_diff(resid,60,0); resid=qaqc_fwd_diff(resid,60,0); resid=qaqc_fwd_diff(resid,60,0); resid=qaqc_fwd_diff(resid,60,0); resid=qaqc_back_diff(resid,60,0); % 300
          resid=qaqc_run_cent_std(resid,20,3);
  
	  if apply_harm==1
            data(isnan(resid))=nan; 	
	  end
          wlv_int=data;
          save([replace(file_obs,' ','_'),'_qaqced.mat'],'time_int','wlv_int')
  	
        elseif apply_qaqc==0
          load([replace(file_obs,' ','_'),'.mat']) % ,'time_obs','wlv')
  
        elseif apply_qaqc==2
          load([replace(file_obs,' ','_'),'_qaqced.mat']) % ,'time_obs','wlv')
          load([replace(file_obs,' ','_'),'_harm.mat']) % ,'time_obs','wlv') 
          resid=wlv_int-data_oharm;
  	
        end
  
      elseif strcmp(file{fobs},'level_whangarei_harbour_at_marsden_point_07_11_24 09_46_58') 
  
        %applying qaqc
        if apply_qaqc==1
          data=wlv_int;
  
          gfilename='level_whangarei_harbour_ginput_qaqc_outputs.txt';
  	  if exist(gfilename,'file')==2;
  	    display(['REMEMBER TO DELETE ',gfilename,', if you re using this file from a previous study.']);	
            data=qaqc_read_ginput_data(time_int,data,gfilename); close;
  	  end
  
          % finding gaps in the original data and inserting nans on the interpolated data
          %[dif iloc]=find(diff(time_obs)>1/24/(60/10)); % 10 minutes
          %for i=1:length(iloc)
          %  [dif iini]=nanmin(abs((time_obs(iloc(i))-time_int)));
          %  [dif ifin]=nanmin(abs((time_obs(iloc(i)+1)-time_int)));
          %  data(iini:ifin)=nan;
          %end
   
          %manual clean up
      	  %bad begining
          %[dif loc]=nanmin(abs(time_int-727858.80943138)); data(1:loc)=nan;
          %[dif loc]=nanmin(abs(time_int-datenum(2005,1,1))); data(1:loc)=nan;
  	  % removing large gaps
          %x1=728424.545195489; x2=733703.496539174; [dif loc1]=nanmin(abs(time_int-x1)); [dif loc2]=nanmin(abs(time_int-x2)); data(1:loc2)=nan;
          %x1=737094.839857142; x2=737111.736473575; [dif loc1]=nanmin(abs(time_int-x1)); [dif loc2]=nanmin(abs(time_int-x2)); data(loc1:loc2)=nan;
          %x1=737867.751552838; x2=737869.350375492; [dif loc1]=nanmin(abs(time_int-x1)); [dif loc2]=nanmin(abs(time_int-x2)); data(loc1:loc2)=nan;
  
          mean_data=nanmean(wlv(:));
          %data=detrend(data,1,'omitnan')+mean_data;
          %wlv=data; time_obs=time_int;

          %% removing extremes
          %data=qaqc_extremes(data,-600,4000);
          %%% removing large local differences
          %%%if strcmp(file{fobs},'Nelson Port at Fairway Beacon') % proc_obs==1;     
          %data=qaqc_fwd_diff(data,300,0); % 300
          %data=qaqc_back_diff(data,300,0); % 300
          %data=qaqc_remove_lin_interp(data,1.000002);
          %%end
  
          % harmonic analysis
	  if apply_harm==1
            if read_harm==0
              loc1=1; 
  	      %[dif loc1]=nanmin(abs(time_int-730355.378757781)); 
              [tidestruc,data_oharm]=t_tide(data(loc1:end),'interval',(time_int(2)-time_int(1))*24,'start time',time_int(loc1),'latitude',-35.766627,'synthesis',1); %,'output','none'););
  	      %data_oharm=t_predic(time_int,tidestruc,'synthesis',1);
              data_oharm=data_oharm+nanmean(data);
              save([replace(file_obs,' ','_'),'_harm.mat'],'data_oharm') 
            else 
              disp(['loading: ',replace(file_obs,' ','_'),'_harm.mat']) % ,'time_obs','wlv') 
              load([replace(file_obs,' ','_'),'_harm.mat']) % ,'time_obs','wlv') 
            end
            resid=data-data_oharm;
            resido=resid;
          end 
          % removing obs with resid surge above 1.5m and below -1.5m
          %resid=qaqc_extremes(resid,-400,1500); 
          %resid=qaqc_fwd_diff(resid,60,0); resid=qaqc_fwd_diff(resid,60,0); resid=qaqc_fwd_diff(resid,60,0); resid=qaqc_fwd_diff(resid,60,0); resid=qaqc_back_diff(resid,60,0); % 300
          %resid=qaqc_fwd_diff(resid,60,0); resid=qaqc_fwd_diff(resid,60,0); resid=qaqc_fwd_diff(resid,60,0); resid=qaqc_fwd_diff(resid,60,0); resid=qaqc_back_diff(resid,60,0); % 300
          %resid=qaqc_fwd_diff(resid,60,0); resid=qaqc_fwd_diff(resid,60,0); resid=qaqc_fwd_diff(resid,60,0); resid=qaqc_fwd_diff(resid,60,0); resid=qaqc_back_diff(resid,60,0); % 300
          %resid=qaqc_run_cent_std(resid,20,3);
  
	  if apply_harm==1
            data(isnan(resid))=nan; 	
	  end
          wlv_int=data;
          save([replace(file_obs,' ','_'),'_qaqced.mat'],'time_int','wlv_int')
  	
        elseif apply_qaqc==0
          load([replace(file_obs,' ','_'),'.mat']) % ,'time_obs','wlv')
  
        elseif apply_qaqc==2
          load([replace(file_obs,' ','_'),'_qaqced.mat']) % ,'time_obs','wlv')
          load([replace(file_obs,' ','_'),'_harm.mat']) % ,'time_obs','wlv') 
          resid=wlv_int-data_oharm;
  	
        end
  
      end % if strcmp(file{fobs},'Nelson Port at Fairway Beacon') % proc_obs==1;
  
  
      % filing original data with harmonic reconstruction when cyclone fehi was present
      fix_fehi=0;    
      if fix_fehi & stations==2
        [dif i1]=nanmin(abs(time_int-737092.385580058)); [dif i2]=nanmin(abs(time_int-737092.468628188));
        wlv_int(i1:i2)=data_oharm(i1:i2)+4736.30167552436-4130.63004584364;
      end

 
      
      % PLOTTING
  
      if plot_data==1

        rmean=movmean(wlv_int,30*24*(60/5));

        % plotting data
        figure('position',[scrsz],'color',[1 1 1],'visible','on');
        hold on; 
        legh={};
        plot(time_obs,wlv,'.r','markersize',10);  legh=[legh,{['obs']}];
        plot(time_int,wlv_int,'.k'); legh=[legh,{['int']}];
        if apply_qaqc==1 & apply_harm==1
          plot(time_int,data_oharm,'.b'); legh=[legh,{['harm']}];
          plot(time_int,resido,'.','color',[1 0 1]); legh=[legh,{['resido']}];
          plot(time_int,resid,'.','color',[0 .7 0]); legh=[legh,{['resid']}];
        elseif apply_qaqc==2 & apply_harm==1
          plot(time_int,data_oharm,'.b'); legh=[legh,{['harm']}];
          resid=wlv_int-data_oharm;
          plot(time_int,resid,'.','color',[0 .7 0]); legh=[legh,{['resid']}];
        end
        plot(time_int,rmean,'.c'); legh=[legh,{['moving mean']}];
        legend(legh)
        %xlim([datenum(2009,4,21) datenum(2009,4,22)]) 
        datetick('x','dd-mmm-yy','keeplimits')
        grid on
        title(replace(file{fobs},'_',' '))
  
        if manual_qaqc>=1
  	data=wlv_int;
  
          inp=input('Zoom in to the selected area. Do you want to delete data using ginput? (y=1):');
          if inp==1
            while inp==1
              data=qaqc_write_ginput_data(time_int,data,gfilename);
              inp=input('Zoom in to the selected area. Do you want to delete data using ginput? (y=1):');
            end
          end
  
          inp=input('Do you want to apply manual corrections to original? (y=1):');
          if inp==1
            %wlv_int(isnan(data))=nan;
            data=qaqc_read_ginput_data(time_int,data,gfilename);
          end
  
          % saving data
          inp=input('Do you want to save the qaqced data? (y=1):');
          if inp==1
            wlv_int=data;
            save([replace(file_obs,' ','_'),'_qaqced.mat'],'time_int','wlv_int') 
          end
        end
  
      end
  
      % calculating tidal design values
      if calc_tides==1
        if use_anomaly==1
          data=data_oharm-nanmean(data_oharm);
        else
          data=data_oharm;
        end
        data=data./1000; % converting to meters
        ddata=diff(data);
  
        % finding mean high water springs (MHWS) (MHWS-10, MHWS-6), Cadastral mean high-water springs (MHWS-C), mean higher high-water (MHHW) and minimum high water (MinHW).
        % Mean High Water Springs (MHWS) is the average height of two successive high tides over a 24-hour period during a semi-lunation, which is about every 14 days.
        % Mean Higher High Water (MHHW) is the average height of the daily high water 
        % Minimum High Water (MinHW) is the average height of the lowest high tide recorded at a tide station each day during a 19-year period.
  
        k=0; kl=0;
        for i=2:length(data)-1
  	  if data(i-1)<data(i) & data(i+1)<=data(i)
              k=k+1; ht(k)=i;
  	  end
  	  % finding low tide
  	  if data(i-1)>data(i) & data(i+1)>=data(i)
  	    kl=kl+1; lt(kl)=i;
  	  end
        end
  
        % loop to find maximum tidal range every 14.8 days
        thigh_tide=time_int(ht); tlow_tide=time_int(lt); 
        high_tide=data(ht); low_tide=data(lt); 
        %tidal_range=high_tide-low_tide;
        %figure; hold on
        %plot(thigh_tide,tidal_range,'.b')

        % finding MHHW
	%creating a thigh_tide of round numbers 
	thigh_tided=floor(thigh_tide);
	%finding the daily max high tide
	k=0;
	for i=thigh_tided(1):thigh_tided(end)
	  k=k+1;
          loc=find(thigh_tided==i);
	  [dif loc2]=nanmax(high_tide(loc));
	  htd(k)=loc(loc2);
	end

        time_fourtn_days=time_int(1):12:time_int(end);
        [dif ismall]=nanmin(high_tide(1:24));
        [dif ini]=nanmin(abs(data-high_tide(ismall)));
        time_fourtn_days=time_int(ini):14.8:time_int(end);
        k=0;
        for i=1:length(time_fourtn_days)-1
  	  [dif iini]=nanmin(abs(thigh_tide-time_fourtn_days(i)));
  	  [dif ifin]=nanmin(abs(thigh_tide-time_fourtn_days(i+1)));
  	  [dif imax]=nanmax(high_tide(iini:ifin));
  	  k=k+1;
  	  imax_tide(k)=imax+iini-1; 
  	  % finding the second max tide
  	  %[datas locs]=sort(high_tide(imax+iini-7:imax+iini-5));
  	  [datas locs]=nanmax(high_tide([imax_tide(k)-1,imax_tide(k)+1]));
  	  if locs==1; ilmax=-1; else ilmax=1; end
  	  k=k+1;
  	  imax_tide(k)=imax_tide(k-1)+ilmax;
        end
        %imax_tide=imax_tide-1;
  
        % 29.6-day maxx tide
        time_fourtn_days=time_int(ini):29.6:time_int(end);
        for i=1:length(time_fourtn_days)-1
  	[dif iini]=nanmin(abs(time_int-time_fourtn_days(i)));
  	[dif ifin]=nanmin(abs(time_int-time_fourtn_days(i+1)));
  	[dif imax]=nanmax(data(iini:ifin));
  	[dif imin]=nanmin(data(iini:ifin));
  	imaxx_tide(i)=imax+iini; imin_tide(i)=imin+iini;
        end
        imaxx_tide=imaxx_tide-1;
  
        % smallest high tide
        [dif imhw]=nanmin(data(ht));
  
        if plot_tides
  
          legh={};
          figure('position',[scrsz],'color',[1 1 1],'visible','on');
          hold on
          plot(time_int,data,'.b'); legh=[legh,{['harm']}];
          plot(time_int(ht),data(ht),'.r','markersize',20); legh=[legh,{['high tides. MHW=',num2str(nanmean(data(ht)),'%.2f')]}];
          plot(thigh_tide(htd),high_tide(htd),'.m','markersize',20); legh=[legh,{['high tides. MHHW=',num2str(nanmean(high_tide(htd)),'%.2f')]}];
          plot(time_int(ini),data(ini),'.c','markersize',20); legh=[legh,{['first low high tide']}];
          plot(thigh_tide(imax_tide),high_tide(imax_tide),'.g','markersize',30); legh=[legh,{['MHWS=',num2str(nanmean(high_tide(imax_tide)),'%.2f')]}];
          plot(time_int(imaxx_tide),data(imaxx_tide),'.k','markersize',20); legh=[legh,{['MHWS-29.6-d=',num2str(nanmean(data(imaxx_tide)),'%.2f')]}];
          plot(time_int(ht(imhw)),data(ht(imhw)),'.c','markersize',40); legh=[legh,{['MinHW=',num2str(nanmean(data(ht(imhw))),'%.2f')]}];
          plot(time_int(lt),data(lt),'.k','markersize',20); legh=[legh,{['low tide']}];
          grid
          legend(legh)
          datetick('x','dd-mmm-yy','keeplimits')
          title(replace(file{fobs},'_',' '))

  
          % finding MHWS-10, MHWS-6, MHWS-C, MHHW and MinHW
          % plotting culative distribution function
          pctl = [94 90 50];
  	% computing hws using two high tides
  	%k=0;
  	%for i=1:2:length(imax_tide)-1; k=k+1; mhws(k)=nanmean(high_tide(imax_tide([i i+1]))); end
          %pctlv = prctile(mhws,pctl); 
          %[f,x] = ecdf(mhws);
          pctlv = prctile(data(ht),pctl); 
          [f,x] = ecdf(data(ht));
          figure
          plot(x,f,'color',[1 0 0],'LineWidth',2)
          hold on
          for k = 1:numel(pctl)
              %plot([1;1]*pctlv(k), [0;1]*pctl(k)/100, '--k')
              %plot([0;1]*pctlv(k), [1;1]*pctl(k)/100, '--k')
              text([1]*pctlv(k),[1]*pctl(k)/100,['MHWS-',num2str(100-pctl(k)),' = ',num2str(pctlv(k),'%.2f')],'fontsize',12,'fontweight','bold')
          end
          title([replace(file{fobs},'_',' '),' - HWS cumulative distribution function'])
          grid on
  
        end
  
        %rmean=data; rmean(isnan(rmean))=nanmean(rmean); 
        %yrmean=movmean(rmean,365*24*(60/5));
        %rmean=movmean(rmean,30*24*(60/5));
        %figure('position',[scrsz],'color',[1 1 1],'visible','on');
        %hold on
        %plot(time_int,rmean-nanmean(rmean))
        %plot(time_int,yrmean-nanmean(rmean),'k')
        %datetick('x','dd-mmm-yy','keeplimits')
        %legend('harmonic moving mean')
        %title(file{fobs})
  
  
      end % if calc_tides==1
  
      % computing annual exceedance probability
      if calc_aep || calc_pot % plot_data==1
  
        %aep=eprob(wlv_int); % ,365*24*(60/5));
        disp(['GEV fitted to annual maxima']);
        disp(['DO NOT click on figures while plots are being made']);
        %figure; clf
        %semilogx([0 0],[0 0],'r',[0 0],[0 0],'b'); hold on;
        ARI = 100./[63 10 5 2 1 0.5 0.2 0.1]; % This specifies the Average Recurrence Intervals (ARI) for output. Change these to what you like.
        fign=98;
        ax=figure(fign); % subplot(lsub,csub,ke);
        %set(gca,'fontsize',12,'fontweight','bold')
        %hold on
  
        % Compare to annual maxima method. Don't need this but sometimes useful to
        data = wlv_int(:)./1000; % full water level
        mean_data=nanmean(data(:));
        %data = wlv_int(:)-data_oharm(:); % residual 
        if detrend_data==1
          data=detrend(data,1,'omitnan'); % +mean_data;
          %mmean=movmean(data,30*24*(60/5),'omitnan'); % +mean_data;
          %data=data-mmean; % +mean_data;
        end
  	% reintroducing the mean
        data=data+mean_data;

        % computing median of the annual maximum
        time_intv=datevec(time_int);
        k=0;
        for i=time_intv(1,1):time_intv(end,1)
          iday=find(time_intv(:,1)==i);
          k=k+1;
          data_max(k)=nanmax(data(iday));
        end
        %data_max=(data_max(:)-mean_data)./1000; % converting to meters
  
        %data=data./1000;
        %mean_data=mean_data./1000;
        % removing the mean
        %data=data-mean_data;
  
        time = time_int(:);
        figure; plot(time,data); datetick('x','dd-mmm-yy','keeplimits')

        if calc_aep==1
  
          AMy = AMcalc(time,data,1); 
          inan=find(isnan(AMy(:,3))); AMy(inan,:)=[]; % removing nans
          AM = AMy(1:end,3);
          [gev,gevci] = gevfit(AM);
          GEV = gevinv(1-ARI2AEP(ARI),gev(1),gev(2),gev(3));
          GEV(2,:) = gevinv(1-ARI2AEP(ARI),gevci(1,1),gevci(1,2),gevci(1,3));
          GEV(3,:) = gevinv(1-ARI2AEP(ARI),gevci(2,1),gevci(2,2),gevci(2,3));
  
          %AM=AM-mean_data;
          %GEV=GEV-mean_data;
  
          figure;[AMari,AMsorted] = gringorten(AM,1,1);  hold on;
          close
  
          figure(fign);
          semilogx(AMari,AMsorted,'xk'); 
          hold on;
          semilogx(ARI,GEV(1,:),'k','LineWidth',2)
          semilogx(ARI,GEV(2,:),'--k',ARI,GEV(3,:),'--k')
          %hold off
          % loop in ARI to plot the values
          for i=1:length(ARI)
            text(ARI(i),GEV(1,i),num2str(GEV(1,i),'%.2f'),'color','r','fontsize',12,'fontweight','bold')
          end
          title(['GEV-AM (red) and GPD-POT (pink) - total water level. Median of annual maxima = ',num2str(nanmedian(data_max),'%.2f'),' m'])
  
        end % if calc_aep==1

        % WATER LEVEL
        if calc_pot==1
          % POT method
          % ERRORS: Confidence intervals and standard errors can not be computed reliably.
          figure
          cdfplot(data) % PF cdf = cumulative density function. find suitable thresholds to constrain GPthres.m. This is just a rough guide. looks like something between 2 and 6 metres would work (between 0.5 and 1)
     	  %return
          tolerance =0.0001; % sometimes gpfit doesn't like having zero values, so add a tiny (inconsequential) amount
          % GPD
          tthr = 3; % peaks over threshold (POT) separated by at least 3 days. Most weather systems in NZ are separated by 4-7 days
          %ythr = GPthres(time,data+tolerance,3.700-(mean_data),4.700-(mean_data),.100,tthr,1); % Automate threshold selection
          if strcmp(file{fobs},'level_hatea_at_town_basin_07_11_24 09_19_28')
	     thr1=3.20-nzvd_off(fobs)./1000; 
	     thr2=3.400-nzvd_off(fobs)./1000; 
	     dthr=.01;
	  
	  elseif strcmp(file{fobs},'level_whangarei_harbour_at_marsden_point_07_11_24 09_46_58')
	     thr1=2-nzvd_off(fobs)./1000; 
	     thr2=3.500-nzvd_off(fobs)./1000; 
	     dthr=.05;
          end  	     

          ythr = GPthres(time,data+tolerance,thr1,thr2,dthr,tthr,1); % Automate threshold selection
          disp(['GPD threshold selected = ' num2str(ythr) ' m']);
          figure;
          [GPD,gpp,gpci,Sy,ARI_POT,POT_sorted] = GPDfit(data+tolerance,time,tthr,ythr,ARI); % PF
          %close
  
          %POT_sorted=POT_sorted-mean_data;
          %GPD=GPD-mean_data;
  
          figure(ax)
          semilogx(ARI_POT,POT_sorted,'o','color','b')
          hold on
          semilogx(ARI,GPD(1,:),'-', 'color','b','LineWidth',2)
          semilogx(ARI,GPD(2,:),'--','color','b')
          semilogx(ARI,GPD(3,:),'--','color','b')
          % loop in ARI to plot the values
          for i=1:length(ARI)
            text(ARI(i),GPD(1,i),['',num2str(GPD(1,i),'%.2f')],'color','m','fontsize',12,'fontweight','bold')
          end
        end
  
        xlim([1 ARI(end)]) 
        set(gca,'fontsize',12,'fontweight','bold')
        xlabel('Average recurrence interval  (years)')
        ylabel('Sea level (m)')
        %ylabel('H_s')
        grid('on')
        set(gcf,'position',[2 42 958 954],'color',[1 1 1],'visible','on');
        title(['EVA of total water level. Median of annual maxima = ',num2str(nanmedian(data_max),'%.2f'),' m'])
        title({['EVA. Median of annual maxima = ',num2str(nanmedian(data_max),'%.2f'),' m'],['GPD threshold selected = ',num2str(ythr)]})


      end % if calc_aep || calc_pot


      if skew_aep
        disp(['GEV fitted to annual maxima']);
        disp(['DO NOT click on figures while plots are being made']);
        %figure; clf
        %semilogx([0 0],[0 0],'r',[0 0],[0 0],'b'); hold on;
        ARI = 100./[63 10 5 2 1 0.5 0.2 0.1]; % This specifies the Average Recurrence Intervals (ARI) for output. Change these to what you like.
        %wlv_int=interp1(time_int,wlv_int,time_int,'linear');
        detrend_resid=1;
        if detrend_resid==1
          %resid=detrend(resid,1,'omitnan'); % +mean_data;
          mean_data=nanmean(wlv_int(:));
          wlv_int=detrend(wlv_int,1,'omitnan')+mean_data;
        end
        %inan=find(wlv_int==nan); wlv_int(inan)=data_oharm(inan);
        wlv_int=wlv_int./1000;
        data_oharm=data_oharm./1000 + nanmean(wlv_int)-nanmean(data_oharm./1000);

        rmean=movmean(wlv_int,30*24*(60/5)); 
        [skew,TW,HW,HT_SS,HT,time_skew,iok,SSHF,tSSHF,HT_SSHF,iHF,SSHF_TW] = skewsurge(time_int,wlv_int,data_oharm,rmean);
        time_skew=time_int(iok);
        % Compare to annual maxima method. Don't need this but sometimes useful to
        data = skew(:); % full water level
        mean_data=nanmean(data(:));
        %data = wlv_int(:)-data_oharm(:); % residual 
        if detrend_data==1
          data=detrend(data,1,'omitnan'); % +mean_data;
          %mmean=movmean(data,30*24*(60/5),'omitnan'); % +mean_data;
          %data=data-mmean; % +mean_data;
        end
        % reintroducing the mean
        data=data+mean_data;
        
        fign=100;
        ax=figure(fign); % subplot(lsub,csub,ke);
        %set(gca,'fontsize',12,'fontweight','bold')
        %hold on
  
        % computing median of the annual maximum
        time_intv=datevec(time_skew);
        k=0;
        for i=time_intv(1,1):time_intv(end,1)
          iday=find(time_intv(:,1)==i);
          k=k+1;
          data_max(k)=nanmax(data(iday));
        end
        data_max=(data_max(:)-mean_data); % converting to meters
  
        % removing the mean
        %data=data-mean_data;
  
        time = time_skew(:);
        figure; plot(time,data); datetick('x','dd-mmm-yy','keeplimits')

        if skew_am==1
  
          AMy = AMcalc(time,data,1); 
          inan=find(isnan(AMy(:,3))); AMy(inan,:)=[]; % removing nans
          AM = AMy(1:end,3);
          [gev,gevci] = gevfit(AM);
          GEV = gevinv(1-ARI2AEP(ARI),gev(1),gev(2),gev(3));
          GEV(2,:) = gevinv(1-ARI2AEP(ARI),gevci(1,1),gevci(1,2),gevci(1,3));
          GEV(3,:) = gevinv(1-ARI2AEP(ARI),gevci(2,1),gevci(2,2),gevci(2,3));
  
          %AM=AM-mean_data;
          %GEV=GEV-mean_data;
  
          figure;[AMari,AMsorted] = gringorten(AM,1,1);  hold on; % figure with the annual maxima sorted
          close
  
          figure(fign);
          semilogx(AMari,AMsorted,'xk'); 
          hold on;
          semilogx(ARI,GEV(1,:),'k','LineWidth',2)
          semilogx(ARI,GEV(2,:),'--k',ARI,GEV(3,:),'--k')
          %hold off
          % loop in ARI to plot the values
          for i=1:length(ARI)
            text(ARI(i),GEV(1,i),num2str(GEV(1,i),'%.2f'),'color','r','fontsize',12,'fontweight','bold')
          end
          title(['GEV-AM (red) and GPD-POT (pink) - total water level. Median of annual maxima = ',num2str(nanmedian(data_max),'%.2f'),' m'])
          title({['EVA. Median of annual maxima = ',num2str(nanmedian(data_max),'%.2f'),' m'],['GPD threshold selected = ',num2str(ythr)]})
  
        end % if calc_aep==1

        % WATER LEVEL
        if skew_pot==1
          % POT method
          % ERRORS: Confidence intervals and standard errors can not be computed reliably.
          figure
          cdfplot(data) % PF cdf = cumulative density function. find suitable thresholds to constrain GPthres.m. This is just a rough guide. looks like something between 2 and 6 metres would work (between 0.5 and 1)
          %return
          tolerance =0.0001; % sometimes gpfit doesn't like having zero values, so add a tiny (inconsequential) amount
          % GPD
          tthr = 3; % peaks over threshold (POT) separated by at least 3 days. Most weather systems in NZ are separated by 4-7 days
          %ythr = GPthres(time,data+tolerance,3.700-(mean_data),4.700-(mean_data),.100,tthr,1); % Automate threshold selection
          if strcmp(file{fobs},'level_hatea_at_town_basin_07_11_24 09_19_28')
             thr1=0.02; thr2=.400; dthr=.01;

	  elseif strcmp(file{fobs},'level_whangarei_harbour_at_marsden_point_07_11_24 09_46_58')
             thr1=0.02; thr2=.400; dthr=.005;

          else

          end  	     

          ythr = GPthres(time,data+tolerance,thr1,thr2,dthr,tthr,1); % Automate threshold selection
          disp(['GPD threshold selected = ' num2str(ythr) ' m']);
          figure;
          [GPD,gpp,gpci,Sy,ARI_POT,POT_sorted] = GPDfit(data+tolerance,time,tthr,ythr,ARI); % PF
          %close
  
          %POT_sorted=POT_sorted-mean_data;
          %GPD=GPD-mean_data;
  
          figure(ax)
          semilogx(ARI_POT,POT_sorted,'o','color','b')
          hold on
          semilogx(ARI,GPD(1,:),'-', 'color','b','LineWidth',2)
          semilogx(ARI,GPD(2,:),'--','color','b')
          semilogx(ARI,GPD(3,:),'--','color','b')
          % loop in ARI to plot the values
          for i=1:length(ARI)
            text(ARI(i),GPD(1,i),['',num2str(GPD(1,i),'%.2f')],'color','m','fontsize',12,'fontweight','bold')
          end
        end
  
        xlim([1 100]) 
        xlim([1 ARI(end)]) 
        set(gca,'fontsize',12,'fontweight','bold')
        xlabel('Average recurrence interval  (years)')
        ylabel('Sea level (m)')
        %ylabel('H_s')
        grid('on')
        set(gcf,'position',[2 42 958 954],'color',[1 1 1],'visible','on');
        title({['EVA of total skew surge. Median of annual maxima = ',num2str(nanmedian(data_max),'%.2f'),' m'],['GPD threshold selected = ',num2str(ythr)]})

      end % if skew_aep
  
  end % for fobs=stations
  
  % display initial and times of water level
  disp(['Initial time of water level: ',datestr(time_int(1))]);
  disp(['Final time of water level: ',datestr(time_int(end))]);

end % if proc_wlv==1

if proc_wave==1

  % loading wave data
  w=open([path,file_ocean{ocean_source}]); 
  w.time_mod=w.time_mod+0.5;
  time_mod=w.time_mod;
  wlv_mod=w.wlv_mod.*1000+nanmean(data_oharm);
  wlv_mod=qaqc_fwd_diff(wlv_mod,200,0); % 300
  wlv_mod=qaqc_back_diff(wlv_mod,200,0); % 300
  
  wlv_modf=wlv_mod; inan=find(isnan(wlv_modf)); wlv_modf(inan)=nanmean(wlv_mod);
  wlv_modf=fourfilt(wlv_modf,length(wlv_mod),w.time_mod(2)-w.time_mod(1),length(w.time_mod)/2,3); % filtering amount in days
  wlv_modf(inan)=nan;
  if ocean_source==1 
    [tidestruc,data_mharm]=t_tide(wlv_mod,'interval',(w.time_mod(2)-w.time_mod(1))*24,'start time',w.time_mod(1),'latitude',-41.504489,'synthesis',1); %,'output','none'););
  else  
    data_mharm=wlv_mod;
  end
  data_mharm=data_mharm+nanmean(wlv_mod);
  resid_mod=wlv_mod-data_mharm;
  
  hs_mod=w.hs_mod;
  tm_mod=w.tm_mod;
  pd_mod=w.pd_mod;

  if plot_wave==1
	  
    % Annual Excedence Probability
    %aep=eprob(wlv_int); % ,365*24*(60/5));
    disp(['DO NOT click on figures while plots are being made']);
    %figure; clf
    %semilogx([0 0],[0 0],'r',[0 0],[0 0],'b'); hold on;
    if ocean_source==1 
      ARI = 100./[63 20 10 5 2 1]; % This specifies the Average Recurrence Intervals (ARI) for output. Change these to what you like.
      min_riv=1; max_riv=5; int_riv=.1; % GPD
      %min_riv=5; max_riv=15; int_riv=.5; % GPD
      tolerance =0.001; % sometimes gpfit doesn't like having zero values, so add a tiny (inconsequential) amount
    else
      ARI = 100./[63 20 10 5 2 1 0.5 0.25]; % This specifies the Average Recurrence Intervals (ARI) for output. Change these to what you like.
      min_riv=.8; max_riv=4; int_riv=.1; % GPD
      %min_riv=5; max_riv=15; int_riv=.5; % GPD
      tolerance =0.001; % sometimes gpfit doesn't like having zero values, so add a tiny (inconsequential) amount
    end

    fign=99;
    ax=figure(fign); % subplot(lsub,csub,ke);
    %set(gca,'fontsize',12,'fontweight','bold')
    %hold on
    
    % Compare to annual maxima method. Don't need this but sometimes useful to
    data = hs_mod; %
    time=time_mod(:);
    
    % computing median of the annual maximum
    time_intv=datevec(time_mod);
    k=0;
    for i=time_intv(:,1):time_intv(end,1)
      iday=find(time_intv(:,1)==i);
      k=k+1;
      if ~isempty(data(iday))		
        data_max(k)=nanmax(data(iday));
      else
	data_max(k)=nan;
      end
    end
    %data_max=(data_max(:)-mean_data)./1000; % converting to meters
    
    AMy = AMcalc(time,data,1); AM = AMy(:,3);
    AM(isnan(AM))=[]; % removing nans
    [gev,gevci] = gevfit(AM);
    GEV = gevinv(1-ARI2AEP(ARI),gev(1),gev(2),gev(3));
    GEV(2,:) = gevinv(1-ARI2AEP(ARI),gevci(1,1),gevci(1,2),gevci(1,3));
    GEV(3,:) = gevinv(1-ARI2AEP(ARI),gevci(2,1),gevci(2,2),gevci(2,3));
    if calc_aep==1
      disp(['GEV fitted to annual maxima']);
      figure;[AMari,AMsorted] = gringorten(AM,1,1);  hold on;
      %close
      figure(fign);
      semilogx(AMari,AMsorted,'x')
      hold on
      semilogx(ARI,GEV(1,:),'k','LineWidth',2)
      semilogx(ARI,GEV(2,:),'--k',ARI,GEV(3,:),'--k')
      % loop in ARI to plot the values
      for i=1:length(ARI)
        text(ARI(i),GEV(1,i),['GEV ',num2str(GEV(1,i),'%.2f')],'color','r','fontsize',12,'fontweight','bold')
      end
      title(['Wave GEV-AM. Median of annual maxima = ',num2str(nanmedian(data_max),'%.2f'),' m^3/s'])
      xlabel('Average recurrence interval  (years)')
      %ylabel('H_s')
      grid('on')
      set(gcf,'position',[2 42 958 954],'color',[1 1 1],'visible','on');
    end
    if calc_pot==1
	    
      figure
      cdfplot(data) % PF cdf = cumulative density function. find suitable thresholds to constrain GPthres.m. This is just a rough guide. looks like something between 2 and 6 metres would work (between 0.5 and 1)
      % GPD
      tthr = 3; % peaks over threshold (POT) separated by at least 3 days. Most weather systems in NZ are separated by 4-7 days
    
      ythr = GPthres(time,data+tolerance,min_riv,max_riv,int_riv,tthr,1); % Automate threshold selection 10,50 worked well
      disp(['GPD threshold selected = ' num2str(ythr) ' ']);
      figure;
      [GPD,gpp,gpci,Sy,ARI_POT,POT_sorted] = GPDfit(data+tolerance,time,tthr,ythr,ARI); % PF
      %close
      figure(fign)
      semilogx(ARI_POT,POT_sorted,'o','color','b'),
      hold on
      semilogx(ARI,GPD(1,:),'-', 'color','b','LineWidth',2)
      semilogx(ARI,GPD(2,:),'--','color','b')
      semilogx(ARI,GPD(3,:),'--','color','b')
      % loop in ARI to plot the values
      for i=1:length(ARI)
        %text(ARI(i),GPD(1,i),num2str(GPD(1,i),'%.2f'),'color','m','fontsize',12,'fontweight','bold')
        text(ARI(i),GPD(1,i),['',num2str(GPD(1,i),'%.2f')],'color','m','fontsize',12,'fontweight','bold')
      end
      title(['Wave EVA. Median of annual maxima = ',num2str(nanmedian(data_max),'%.2f'),' m^3/s'])
      xlabel('Average recurrence interval  (years)')
      %ylabel('H_s')
      grid('on')
      set(gcf,'position',[2 42 958 954],'color',[1 1 1],'visible','on');
    end
    
    %xlim([1 100]) 
    set(gca,'fontsize',12,'fontweight','bold')
    return

  end % if plot_wave==1
	
end % if proc_wave==1



if proc_river

  % processing river flow data
  for fobs=river_stations % 1:length(file)
  
      file_robs=[path,file_river{fobs}]; tic;
      if read_river==1
        disp(['Processing csv: ',file_robs]);
        ob=importdata([file_robs,'.csv']);
        txt=ob.textdata;
        time_river=nan(length(ob.textdata)-2,1);
        for j=2:length(ob.textdata)
            time_river(j-1)=datenum([txt{j,2}],'yyyy-mm-dd HH:MM:SS');
        end
        river=ob.data(:,1); % m3/s
        save([replace(file_robs,' ','_'),'.mat'],'time_river','river') 

      else
        disp(['Loading matlab: ',file_robs]);
        load([replace(file_robs,' ','_'),'.mat']) % 

        if strcmp(file_river{fobs},'Wairau at Tuamarina 5 minute instantaneous flows')

          time_rivern=time_river; rivern=river;

          load([path_flow,'Wairau_at_Barnetts_Bank_5_minute_instantaneous_flows.mat'])
	  river=river*0.001; % converting to m3/s
	  [dif loc]=nanmin(abs(time_river-datenum(1999,7,30)));
	  river=river(loc:end); time_river=time_river(loc:end);
	  % cocatenating the two datasets
	  river=[rivern;river]; time_river=[time_rivern;time_river];

        end

      end	
      toc

      % processing observations
      % interpolating to 1 minute interval
      %time_int=[time_obs(1):1/24/(60/5):time_obs(end)]'; % 5 minutes
      %wlv=wlv-nanmean(wlv(:));
      %wlv_int=interp1(time_obs,wlv,time_int);

      % display initial and times of river flow
      disp(['Initial time of river flow: ',datestr(time_river(1))]);
      disp(['Final time of river flow: ',datestr(time_river(end))]);

  
      if plot_river==1;
  
        legh={};
        figure;
        figure('position',[scrsz],'color',[1 1 1],'visible','on');
        set(gca,'fontsize',12,'fontweight','bold')
        hold on
        plot(time_river,river,'.b'); legh=[legh,{['raw']}];
  
        grid
        legend(legh)
        datetick('x','dd-mmm-yy','keeplimits')
        title(replace(file_river{fobs},'_',' '))
        xlabel('Flow (m^3/s)')

	% Annual Excedence Probability
        %aep=eprob(wlv_int); % ,365*24*(60/5));
        disp(['DO NOT click on figures while plots are being made']);
        %figure; clf
        %semilogx([0 0],[0 0],'r',[0 0],[0 0],'b'); hold on;
        ARI = 100./[63 10 5 2 1 0.5 0.2 0.1]; % This specifies the Average Recurrence Intervals (ARI) for output. Change these to what you like.
        fign=99;
        ax=figure(fign); % subplot(lsub,csub,ke);
        %set(gca,'fontsize',12,'fontweight','bold')
        %hold on

        % Compare to annual maxima method. Don't need this but sometimes useful to
        data = river; % full water level
        %data = wlv_int(:)-data_oharm(:); % storm surge 
        %detrend_data=1;
        if detrend_data==1
          %data=detrend(data,1,'omitnan'); % +mean_data;
        end
        time=time_river(:);

        % computing median of the annual maximum
        time_intv=datevec(time_river);
        k=0;
        for i=time_intv(1,1):1:time_intv(end,1)
          iday=find(time_intv(:,1)==i);
          k=k+1;
          if ~isempty(data(iday))		
            data_max(k)=nanmax(data(iday));
          else
            data_max(k)=nan;
          end
        end
        %data_max=(data_max(:)-mean_data)./1000; % converting to meters

        if calc_aep==1
          disp(['GEV fitted to annual maxima']);
          AMy = AMcalc(time,data,1); AM = AMy(:,3);
	  AM(isnan(AM))=[]; % removing nans
          [gev,gevci] = gevfit(AM);
          GEV = gevinv(1-ARI2AEP(ARI),gev(1),gev(2),gev(3));
          GEV(2,:) = gevinv(1-ARI2AEP(ARI),gevci(1,1),gevci(1,2),gevci(1,3));
          GEV(3,:) = gevinv(1-ARI2AEP(ARI),gevci(2,1),gevci(2,2),gevci(2,3));

          figure;[AMari,AMsorted] = gringorten(AM,1,1);  hold on;
          %close
          figure(fign);
          semilogx(AMari,AMsorted,'x')
          hold on
          semilogx(ARI,GEV(1,:),'k','LineWidth',2)
          semilogx(ARI,GEV(2,:),'--k',ARI,GEV(3,:),'--k')
          % loop in ARI to plot the values
          for i=1:length(ARI)
            text(ARI(i),GEV(1,i),['GEV ',num2str(GEV(1,i),'%.2f')],'color','r','fontsize',12,'fontweight','bold')
          end
          title([replace(file_river{fobs},'_',' '),' GEV-AM. Median of annual maxima = ',num2str(nanmedian(data_max),'%.2f'),' m^3/s'])
          xlabel('Average recurrence interval  (years)')
          %ylabel('H_s')
          grid('on')
          set(gcf,'position',[2 42 958 954],'color',[1 1 1],'visible','on');
        end

        % RIVER
        if calc_pot==1
	  data=river;

          figure
          cdfplot(data) % PF cdf = cumulative density function. find suitable thresholds to constrain GPthres.m. This is just a rough guide. looks like something between 2 and 6 metres would work (between 0.5 and 1)
          tolerance =0.1; % sometimes gpfit doesn't like having zero values, so add a tiny (inconsequential) amount
          % GPD
          tthr = 3; % peaks over threshold (POT) separated by at least 3 days. Most weather systems in NZ are separated by 4-7 days

	  if strcmp(file_river{fobs},'flow_hatea_at_whareora_rd_07_11_24 14_31_58')
            min_riv=9; max_riv=60; int_riv=1;
            min_riv=9; max_riv=30; int_riv=1;
	  elseif strcmp(file_river{fobs},'flow_raumanga_at_bernard_st_07_11_24 14_29_44')
            min_riv=2; max_riv=20; int_riv=.5;
	  elseif strcmp(file_river{fobs},'flow_waiarohia_at_lovers_lane_07_11_24 14_35_58')
            min_riv=2; max_riv=20; int_riv=1;
	  end

          ythr = GPthres(time,data+tolerance,min_riv,max_riv,int_riv,tthr,1); % Automate threshold selection 10,50 worked well
          disp(['GPD threshold selected = ' num2str(ythr) ' ']);
          figure;
          [GPD,gpp,gpci,Sy,ARI_POT,POT_sorted] = GPDfit(data+tolerance,time,tthr,ythr,ARI); % PF
          %close
          figure(fign)
          semilogx(ARI_POT,POT_sorted,'o','color','b'),
          hold on
          semilogx(ARI,GPD(1,:),'-', 'color','b','LineWidth',2)
          semilogx(ARI,GPD(2,:),'--','color','b')
          semilogx(ARI,GPD(3,:),'--','color','b')
          % loop in ARI to plot the values
          for i=1:length(ARI)
            %text(ARI(i),GPD(1,i),num2str(GPD(1,i),'%.2f'),'color','m','fontsize',12,'fontweight','bold')
            text(ARI(i),GPD(1,i),['',num2str(GPD(1,i),'%.2f')],'color','m','fontsize',12,'fontweight','bold')
          end
          title([replace(file_river{fobs},'_',' '),' EVA. Median of annual maxima = ',num2str(nanmedian(data_max),'%.2f'),' m^3/s'])
          title({['EVA. Median of annual maxima = ',num2str(nanmedian(data_max),'%.2f'),' m'],['GPD threshold selected = ',num2str(ythr)]})
          xlabel('Average recurrence interval  (years)')
          %ylabel('H_s')
          grid('on')
          set(gcf,'position',[2 42 958 954],'color',[1 1 1],'visible','on');
        end

%        return
        %xlim([1 100]) 
        xlim([1 ARI(end)]) 
        set(gca,'fontsize',12,'fontweight','bold')
      end
      
  end

end



if proc_rainfall

  % processing river flow data
  for fobs=rain_stations % 1:length(file)
  
      file_robs=[path,file_rainfall{fobs}]; tic;
      if read_rainfall==1
        disp(['Processing csv: ',file_robs]);
        ob=importdata([file_robs,'.csv']);
        txt=ob.textdata;
        time_rain=nan(length(ob.textdata)-2,1);
        for j=2:length(ob.textdata)
            time_rain(j-1)=datenum([txt{j,2}],'yyyy-mm-dd HH:MM:SS');
        end
        rain=ob.data(:,1); % m3/s

	% computing accumulated daily rainfall in mm
        disp('computing accumulated daily rainfall in mm');
	time_raind=round(time_rain(1)):1:(time_rain(end));
	raind=nan(length(time_raind),1);
	for i=1:length(time_raind)-1
          [dif loc1]=nanmin(abs(time_raind(i)-time_rain));
	  [dif loc2]=nanmin(abs(time_raind(i+1)-time_rain));
          loc2=loc2-1;
	  if abs(time_raind(i)-time_rain(loc1))<6/24 & abs(time_rain(loc2)-time_raind(i+1))<6/24
	    raind(i)=nansum(rain(loc1:loc2));
	  else
            raind(i)=nan;
	  end
	end
	%time_rain=time_raind; rain=raind;

        save([replace(file_robs,' ','_'),'.mat'],'time_rain','rain','time_raind','raind')
        return

      else

        disp(['Loading matlab: ',file_robs]);
        load([replace(file_robs,' ','_'),'.mat']) % 

      end	
      toc

      % processing observations
      % interpolating to 1 minute interval
      %time_int=[time_obs(1):1/24/(60/5):time_obs(end)]'; % 5 minutes
      %wlv=wlv-nanmean(wlv(:));
      %wlv_int=interp1(time_obs,wlv,time_int);

      % display initial and times of river flow
      disp(['Initial time of rain fall: ',datestr(time_rain(1))]);
      disp(['Final time of rain fall: ',datestr(time_rain(end))]);
  
      if daily_rain==1;
        disp('using accumulated daily rainfall in mm');
        time_rain=time_raind; rain=raind;
      end

      %QAQC
      if strcmp(file_rainfall{fobs},'rainfall_hatea_at_glenbervie_forest_hq_07_11_24 11_19_21')
	% inserting nan before 15th of December 1987      
	[dif loc]=nanmin(abs(time_rain-datenum(1987,12,15)));      
	rain(1:loc)=nan;
	[dif loc]=nanmin(abs(time_raind-datenum(1987,12,15)));
	raind(1:loc)=nan;
      end

      if plot_rainfall==1;
  
        legh={};
        figure;
        figure('position',[scrsz],'color',[1 1 1],'visible','on');
        set(gca,'fontsize',12,'fontweight','bold')
        hold on
        plot(time_rain,rain,'.b'); legh=[legh,{['raw']}];
  
        grid
        legend(legh)
        datetick('x','dd-mmm-yy','keeplimits')
        title(replace(file_rainfall{fobs},'_',' '))
        xlabel('Rainfall (mm)')

	% Annual Excedence Probability
        %aep=eprob(wlv_int); % ,365*24*(60/5));
        disp(['DO NOT click on figures while plots are being made']);
        %figure; clf
        %semilogx([0 0],[0 0],'r',[0 0],[0 0],'b'); hold on;
        ARI = 100./[63 10 5 2 1 0.5 0.2 0.1]; % This specifies the Average Recurrence Intervals (ARI) for output. Change these to what you like.
        ax=figure; %(fign); % subplot(lsub,csub,ke);
        %set(gca,'fontsize',12,'fontweight','bold')
        %hold on

        % Compare to annual maxima method. Don't need this but sometimes useful to
        data = rain; % full water level
        %data = wlv_int(:)-data_oharm(:); % storm surge 
        %detrend_data=1;
        if detrend_data==1
          %data=detrend(data,1,'omitnan'); % +mean_data;
        end
        time=time_rain(:);

        % computing median of the annual maximum
        time_intv=datevec(time_rain);
        k=0;
        for i=time_intv(1,1):1:time_intv(end,1)
          iday=find(time_intv(:,1)==i);
          k=k+1;
          if ~isempty(data(iday))		
            data_max(k)=nanmax(data(iday));
          else
            data_max(k)=nan;
          end
        end
        %data_max=(data_max(:)-mean_data)./1000; % converting to meters

        if calc_aep==1
          disp(['GEV fitted to annual maxima']);
          AMy = AMcalc(time,data,1); AM = AMy(:,3);
	  AM(isnan(AM))=[]; % removing nans
          [gev,gevci] = gevfit(AM);
          GEV = gevinv(1-ARI2AEP(ARI),gev(1),gev(2),gev(3));
          GEV(2,:) = gevinv(1-ARI2AEP(ARI),gevci(1,1),gevci(1,2),gevci(1,3));
          GEV(3,:) = gevinv(1-ARI2AEP(ARI),gevci(2,1),gevci(2,2),gevci(2,3));

          figure;[AMari,AMsorted] = gringorten(AM,1,1);  hold on;
          %close
          figure(ax);
          semilogx(AMari,AMsorted,'xk')
          hold on
          semilogx(ARI,GEV(1,:),'k','LineWidth',2)
          semilogx(ARI,GEV(2,:),'--k',ARI,GEV(3,:),'--k')
          % loop in ARI to plot the values
          for i=1:length(ARI)
            text(ARI(i),GEV(1,i),['GEV ',num2str(GEV(1,i),'%.2f')],'color','r','fontsize',12,'fontweight','bold')
          end
          title([replace(file_rainfall{fobs},'_',' '),' GEV-AM. Median of annual maxima = ',num2str(nanmedian(data_max),'%.2f'),' mm/day'])
          xlabel('Average recurrence interval  (years)')
          %ylabel('H_s')
          grid('on')
          set(gcf,'position',[2 42 958 954],'color',[1 1 1],'visible','on');
        end

        if calc_pot==1
	  data=rain;

          figure
          cdfplot(data) % PF cdf = cumulative density function. find suitable thresholds to constrain GPthres.m. This is just a rough guide. looks like something between 2 and 6 metres would work (between 0.5 and 1)
          tolerance =0.1; % sometimes gpfit doesn't like having zero values, so add a tiny (inconsequential) amount
          % GPD
          tthr = 3; % peaks over threshold (POT) separated by at least 3 days. Most weather systems in NZ are separated by 4-7 days

	  if strcmp(file_rainfall{fobs},'rainfall_hatea_at_glenbervie_forest_hq_07_11_24 11_19_21')
            min_riv=45; max_riv=100; int_riv=5;
	  elseif strcmp(file_rainfall{fobs},'rainfall_taika_at_cemetery_road_07_11_24 11_17_58')
            min_riv=20; max_riv=60; int_riv=5;
	  elseif strcmp(file_rainfall{fobs},'rainfall_waiarohia_at_nrc_water_st_07_11_24 11_14_29')
            min_riv=45; max_riv=100; int_riv=1;
	  end

          ythr = GPthres(time,data+tolerance,min_riv,max_riv,int_riv,tthr,1); % Automate threshold selection 10,50 worked well
          disp(['GPD threshold selected = ' num2str(ythr) ' ']);
          figure;
          [GPD,gpp,gpci,Sy,ARI_POT,POT_sorted] = GPDfit(data+tolerance,time,tthr,ythr,ARI); % PF
          %close
          figure(ax)
          semilogx(ARI_POT,POT_sorted,'o','color','b'),
          hold on
          semilogx(ARI,GPD(1,:),'-', 'color','b','LineWidth',2)
          semilogx(ARI,GPD(2,:),'--','color','b')
          semilogx(ARI,GPD(3,:),'--','color','b')
          % loop in ARI to plot the values
          for i=1:length(ARI)
            %text(ARI(i),GPD(1,i),num2str(GPD(1,i),'%.2f'),'color','m','fontsize',12,'fontweight','bold')
            text(ARI(i),GPD(1,i),['',num2str(GPD(1,i),'%.2f')],'color','m','fontsize',12,'fontweight','bold')
          end
          title([replace(file_rainfall{fobs},'_',' '),' EVA. Median of annual maxima = ',num2str(nanmedian(data_max),'%.2f'),' mm/day'])
          title({['EVA. Median of annual maxima = ',num2str(nanmedian(data_max),'%.2f'),' m'],['GPD threshold selected = ',num2str(ythr)]})
          xlabel('Average recurrence interval  (years)')
          %ylabel('H_s')
          grid('on')
          set(gcf,'position',[2 42 958 954],'color',[1 1 1],'visible','on');
        end

%        return
        %xlim([1 100]) 
        xlim([1 ARI(end)]) 
        set(gca,'fontsize',12,'fontweight','bold')
      end


  end

end


if crop_river_level==1;

  %wlv_int=interp1(time_int,wlv_int,time_int,'linear');
  detrend_resid=1;
  if detrend_resid==1
    %resid=detrend(resid,1,'omitnan'); % +mean_data;
    mean_data=nanmean(wlv_int(:));
    wlv_int=detrend(wlv_int,1,'omitnan')+mean_data;
  end
  wlv_int=wlv_int./1000;
  data_oharm=data_oharm./1000 + nanmean(wlv_int)-nanmean(data_oharm./1000);
  %inan=find(wlv_int==nan); wlv_int(inan)=data_oharm(inan);

  %wlv_mod=wlv_mod./1000;
  %data_mharm=data_mharm./1000 - nanmean(data_mharm./1000);
  %stormsurge=wlv_mod-data_mharm;
  %if ocean_source==1
  rmean=movmean(wlv_int,30*24*(60/5)); 
  %else
  %  rmean=movmean(wlv_mod,30*24*(60/60)); % 30-day filter (Hindcast data is 60 minutes)
  %end
  [skew,TW,HW,HT_SS,HT,time_skew,iok,SSHF,tSSHF,HT_SSHF,iHF,SSHF_TW] = skewsurge(time_int,wlv_int,data_oharm,rmean);
  %wlv_int=data_oharm; % using harmonic instead of total water level
  time_skew=time_int(iok);
  %skew(skew<-.5)=nan; skew(skew>.7)=nan;
  HT=HT_SS; % high tide

  % cropping series to the same time period
  disp('cropping series to maintain only the overlapping period'); tic
  time_crop=[min([time_river(1),time_int(1),time_rain(1)]):max([time_river(end),time_int(end),time_rain(end)])]';
  disp(['time_crop covers from ',datestr(time_crop(1)),' to ',datestr(time_crop(end))]);

  [dif i1]=nanmin(abs(time_crop(1)-time_rain));
  [dif i2]=nanmin(abs(time_crop(end)-time_rain));
  rain=rain(i1:i2);
  time_rain=time_rain(i1:i2);

  [dif i1]=nanmin(abs(time_crop(1)-time_river));
  [dif i2]=nanmin(abs(time_crop(end)-time_river));
  river=river(i1:i2);
  time_river=time_river(i1:i2);

  %[dif ii1]=nanmin(abs(time_crop(1)-time_int));
  %[dif ii2]=nanmin(abs(time_crop(end)-time_int));
  %wlv_int=data_oharm(ii1:ii2);
  %time_int=time_int(ii1:ii2);

  % HT, skew, and time_skew 
  wlv_int=HT; time_int=time_skew; % using high tide instead of total water level

  [dif ii1]=nanmin(abs(time_crop(1)-time_skew));
  [dif ii2]=nanmin(abs(time_crop(end)-time_skew));
  skew=skew(ii1:ii2);
  time_skew=time_skew(ii1:ii2);
  wlv_int=wlv_int(ii1:ii2); % using high tide instead of total water level
  time_int=time_int(ii1:ii2);

  toc
  

  disp('finding daily maxima and plotting data'); tic
  if plot_daily_max==1
          
    % integer time
    time_cp=time_crop; % unique(floor(time_river));
    % reducing daily time to remove first and last day
    time_cp=time_cp(2:end-1);
    t_cvec=datevec(time_cp);
    time_rdvec=datevec(time_river);
    time_intvec=datevec(time_rain); % rain
    time_skewv=datevec(time_skew);
  
    riverd=nan(length(time_cp),1); skewd=nan(length(time_cp),1); wlvd=nan(length(time_cp),1); raind=nan(length(time_cp),1); 
    for i=1:length(time_cp)
  
      iday=find(t_cvec(i,1)==time_rdvec(:,1) & t_cvec(i,2)==time_rdvec(:,2) & t_cvec(i,3)==time_rdvec(:,3));
      if ~isempty(iday)
        riverd(i)=nanmax(river(iday)); 
      else
        riverd(i)=nan;
      end
  
      iday=find(t_cvec(i,1)==time_intvec(:,1) & t_cvec(i,2)==time_intvec(:,2) & t_cvec(i,3)==time_intvec(:,3)); 
      if ~isempty(iday)
        raind(i)=nanmax(rain(iday));
      else
        raind(i)=nan;
      end

      iday=find(t_cvec(i,1)==time_skewv(:,1) & t_cvec(i,2)==time_skewv(:,2) & t_cvec(i,3)==time_skewv(:,3)); 
      if ~isempty(iday)
        wlvd(i)=nanmax(wlv_int(iday));
        skewd(i)=nanmax(skew(iday));
      else
        wlvd(i)=nan;
        skewd(i)=nan;
      end
  
    end
    toc
  
    % cocatenating river, skew and high tide
    data=[riverd(:),skewd(:),wlvd(:),raind(:)]; % 
  
    % plot daily maxima
    figure('position',[scrsz],'color',[1 1 1],'visible','on');
    % in three subplots
    for i=1:length(data(1,:))
      subplot(length(data(1,:)),1,i)
      plot(time_cp,data(:,i))
      if i==2 || i==3; title(file{stations}); 
      elseif i==1; title(file_river{rr});
      elseif i==4; title(file_rainfall{ra});
      end
      xlim([time_cp(1) time_cp(end)])
      datetick('x','dd-mmm-yy','keeplimits')
      grid on
    end
  end
  
  
  % Writing data to csv files for JP analysis
  if write_jp_file
  
    file_obs=[path_jp,lname{ll},'_',river_name{rr},'.csv'];
    disp(['Processing: ',file_obs]);
  
    fda=file_obs; % [filename(1:end-4),'.csv'];
    display(['Saving: ',fda]);
    if exist(fda)==2
       system(['del ',fda]);
    end
  
    % columns in the generated file must be: year, month, day, river flow, skew surge, hight tide = ("yr", "mo", "dy", "flow", "SS","HT") 
    % DAILY MAX of river flow, skew surge, and high tide
    % there shouldn't be a header in the file
  
    %pdd=deg2rad(pdd); % converting to radians
    %if ocean_source==1
    data=[t_cvec(:,1:3),data]; % subtracting the mean water level to reach NZVD2016
    %else
    %  data=[t_cvec(:,1:3),riverd(:),hsd(:),tmd(:),pdd(:)]; % subtracting the mean water level to reach NZVD2016
    %end
  
    csvdata=data;
    csvdata=num2cell(csvdata);
    %for i=1:length(time_mod); 
    %  csvtime{i,1}=datestr(time_mod(i),'yyyy-mm-ddTHH:MM:SS'); 
    %end
    csvhead={'year','month','day','flow','ss','ht','rain'};
  
    % concatenating str matrices
    csvall=[csvdata];
    %csvall=[csvtime,csvdata];
    csvall=[csvhead;csvall];
    %csvwrite(fda,csvall);
    T = cell2table(csvall); % ,'VariableNames',csvhead);
    writetable(T,fda,'WriteVariableNames',0)
    fclose('all');
  
  end


end




if plot_all

  for nvar=1
	  
    if nvar==1
      data=wlv_int;
      time_data=time_int;
    elseif nvar==2
      data=rain;
      time_data=time_rain;
    elseif nvar==3
      data=river;
      time_data=time_river;
    end
    if ll==3
      % find data during cyclone gabrielle
      [dif loc]=nanmin(abs(time_data-datenum(2023,2,13)));
    else
      % finding max data
      [dif loc]=nanmax(data);
    end

    % plotting sea level, rain fall and river flow with focus on maximum sea level centered on a 4-day window
    figure('position',[scrsz],'color',[1 1 1],'visible','on');
    hold on
    for subp=1:3
      subplot(3,1,subp); 
      set(gca,'fontsize',16,'fontweight','bold')
      hold on
      if subp==1
        plot(time_int,data_oharm./1000,'b','linewidth',2)
        plot(time_int,wlv_int./1000,'k','linewidth',2)
        %plot(time_int(loc),wlv_int(loc)./1000,'.r','markersize',20)
        ylabel('Sea level (m)')
        legend('Harmonic','Observed')%,'Maximum')
        tname='Sea level (m) at Hatea River near Town Basin';
        tname='Sea level (m) at Whangarei Harbour at Marsden Point';
      elseif subp==2
        plot(time_rain,rain,'.-k','linewidth',2)
        ylabel('Rainfall (mm/hour)')
	yyaxis right
        plot(time_raind+.5,raind,'.b','linewidth',2,'markersize',20)
        ylabel('Rainfall (mm/day)')
	set(gca,'ycolor','b')
	yyaxis left
	legend('mm/hour','mm/day')
        if ra==1
          tname='Rainfall at Hatea at Glenbervie Forest HQ';
        elseif ra==2
          tname='Rainfall at Otaika at Cemetery Rd';
        elseif ra==3
          tname='Rainfall at NRC Water St';
        end
      elseif subp==3
        plot(time_river,river,'.-k','linewidth',2)
        ylabel('River flow (m^3/s)')
 	if river_stations==1
          tname='Hatea River flow (m^3/s) at Whareora Rd';
	elseif river_stations==2
          tname='Raumanga Stream flow (m^3/s) at Bernard St';
	elseif river_stations==3
          tname='Waiarohia Stream flow (m^3/s) at Lovers Lane';
	end
      end
      title(tname)

      xlim([time_data(loc)-5,time_data(loc)+5])
      datetick('x','dd-mmm-yy','keeplimits')
      grid on
    end


    if save_fig==1
      fname=[path_fig,'timeseries_1_rr_',num2str(rr),'_river_',river_name{river_stations},'_max_var_',num2str(nvar),'_',file{stations},'.png'];
      disp(['Saving: ',fname])
      export_fig(gcf,fname,'-png','-r300');
      %print('-dpng','-r300',[fname])
    end
    
  end

end

return





if proc_river==1 & plot_river_ocean==1 & crop_river_level==0

  [dif loc]=nanmin(abs(time_river-time_int(1)));
  river=river(loc:end); time_river=time_river(loc:end);
  % removing repeated times from river data
  [time_river,ia,ic]=unique(time_river);
  river=river(ia);
  riveri=interp1(time_river,river,time_int);

  % plotting scatter plot of resid and riveri
  resid=wlv_int-data_oharm;

  % filtering the residual in 73h using fourfilt
  inan=find(isnan(resid)); residf=resid; residf(inan)=nanmean(residf);
  residf=fourfilt(residf,length(resid),time_int(2)-time_int(1),length(time_int)/2,1);
  residf(inan)=nan;
  data_oharm(inan)=nan; 

  if predict_river_wlv==1

    % plotting scatter plot of resid and riveri
    figure('position',[scrsz],'color',[1 1 1],'visible','on');
    hold on
    scatter(riveri,resid)
    xlabel('River flow (m^3/s)')
    ylabel('Residual (m)')
    title(replace(file{fobs},'_',' '))
    % computing correlation coefficient and linear regression model and texting on the figure
    [r,p]=corrcoef(riveri,resid,'rows','complete');
    text(0.1,0.9,['r = ',num2str(r(1,2),'%.2f')],'Units','normalized','fontsize',12,'fontweight','bold')
    text(0.1,0.8,['p = ',num2str(p(1,2),'%.2f')],'Units','normalized','fontsize',12,'fontweight','bold')
    % computing linear regression model
    mdl = fitlm(riveri,resid);
    % plotting the linear regression model
    plot(mdl)
    % displaying the linear regression model on the figure
    text(0.1,0.7,['y = ',num2str(mdl.Coefficients.Estimate(1),'%.2f'),' + ',num2str(mdl.Coefficients.Estimate(2),'%.2f'),'x'],'Units','normalized','fontsize',12,'fontweight','bold')
    % predicting water level from river flow
    wlvr=mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)*riveri;

    %plotting scatter plot of residf and riveri
    figure('position',[scrsz],'color',[1 1 1],'visible','on');
    hold on
    scatter(riveri,residf)
    xlabel('River flow (m^3/s)')
    ylabel('Residual (m)')
    title(replace(file{fobs},'_',' '))
    % computing correlation coefficient and linear regression model and texting on the figure
    [r,p]=corrcoef(riveri,residf,'rows','complete');
    text(0.1,0.9,['r = ',num2str(r(1,2),'%.2f')],'Units','normalized','fontsize',12,'fontweight','bold')
    text(0.1,0.8,['p = ',num2str(p(1,2),'%.2f')],'Units','normalized','fontsize',12,'fontweight','bold')
    % computing linear regression model
    mdl = fitlm(riveri,residf);
    % plotting the linear regression model
    plot(mdl)
    % displaying the linear regression model on the figure
    text(0.1,0.7,['y = ',num2str(mdl.Coefficients.Estimate(1),'%.2f'),' + ',num2str(mdl.Coefficients.Estimate(2),'%.2f'),'x'],'Units','normalized','fontsize',12,'fontweight','bold')
    % predicting water level from river flow
    wlvrf=mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)*riveri;

    % log fit
    figure('position',[scrsz],'color',[1 1 1],'visible','on');
    hold on
    scatter(riveri,residf)
    xlabel('River flow (m^3/s)')
    ylabel('Residual (m)')
    title(replace(file{fobs},'_',' '))
    % remove nans from the data
    inan=find(isnan(residf)); residfn=residf; residfn(inan)=[]; riverin=riveri; riverin(inan)=[];
    inan=find(isnan(riverin)); residfn(inan)=[]; riverin(inan)=[];
    % Fit a logarithmic model
    fit_residfl = fit(riverin, residfn, 'a*log(x) + b');
    plot(fit_residfl, 'r-'); % Fitted curve
    wlvrfl = feval(fit_residfl, riveri);

    % plotting the residuals
    %resid=resid-predict(mdl,riveri);
    %plot(riveri,resid,'.r')
  
  end


  close all

  % Plotting data
  figure('position',[scrsz],'color',[1 1 1],'visible','on');
  hold on; 
  legh={};

  yyaxis left
  %ylabel('River flow (m^3/s)')


  % making the axis color the same as the line
  ax = gca;
  ax.YAxis(1).Color = 'k';
  ax.YAxis(2).Color = 'k';

  ylabel('Water level (mm) and River flow (m^3/s)')
  %plot(time_obs,wlv,'.r');  legh=[legh,{['obs']}];
  plot(time_int,wlv_int,'.k'); legh=[legh,{['int']}];
  %if apply_qaqc==2
  plot(time_int,data_oharm,'.b'); legh=[legh,{['harm']}];
  %plot(time_int,resido,'.','color',[1 0 1]); legh=[legh,{['resido']}];
  %elseif apply_qaqc==1
  %  plot(time_int,data_oharm,'.b'); legh=[legh,{['harm']}];
 
  if predict_river_wlv==1
    plot(time_int,wlvr,'.','color',[.7 .7 0]); legh=[legh,{['wlv river']}];
    plot(time_int,wlvrf,'.','color',[.5 .5 .5]); legh=[legh,{['wlv river (residf)']}];
    %plot(time_int,wlvrfl,'.','color',[.7 .7 .7]); legh=[legh,{['wlv river (residf)']}];
    plot(time_int,residf,'.','color',[1 0 1]); legh=[legh,{['resid filtered']}];
    plot(time_int,residf-wlvr,'.','color',[0 1 1]); legh=[legh,{['ocean (residf-wlvr)']}];
    plot(time_int,residf-wlvrf,'.','color',[1 0 0]); legh=[legh,{['ocean (residf-wlvrf)']}];
    %plot(time_int,residf-wlvrfl,'.','color','b'); legh=[legh,{['ocean (residf-wlvrfl)']}];

  end

  %plot(time_int,rmean,'.c'); legh=[legh,{['moving mean']}];

  plot(w.time_mod,data_mharm,'-','color','r','linewidth',1); legh=[legh,{['harm mod']}];
  plot(w.time_mod,wlv_mod,'-','color',[1 0 1],'linewidth',2); legh=[legh,{['wlv mod']}];

  plot(time_int,resid,'.','color',[0 .7 0]); legh=[legh,{['resid']}];
  %plot(w.time_mod,wlv_modf,'-','color',[0 1 1],'linewidth',2); legh=[legh,{['wlv modf']}];
  plot(w.time_mod,resid_mod,'-','color','y','linewidth',2); legh=[legh,{['resid mod']}];

  plot(time_river,river,'.','color',[.3 .3 .3]); legh=[legh,{['river']}];

  % plotting on the right-hand axis
  yyaxis right
  plot(w.time_mod,w.hs_mod,'.','color',[146 36 40]./255); legh=[legh,{['hs mod']}];
  %plot(w.time_mod,w.hs_mod,'.','color','r'); legh=[legh,{['hs mod']}];
  %plot(w.time_mod,w.tp_mod,'.','color',[0 0 1]); legh=[legh,{['tp mod']}];
  ylabel('Wave height (m) and peak period (s)')

  yyaxis left
  %xlim([datenum(2009,4,21) datenum(2009,4,22)]) 
  xlim(w.time_mod([1 end]))	
  datetick('x','dd-mmm-yy','keeplimits')
  grid on
  legend(legh)
  title(replace(file{fobs},'_',' '))


end


return

resid_mod=resid_mod;
%oceanf=w.tm_mod;
%oceanf=interp1(time_int,residf-wlvrf,w.time_mod);
oceanf=interp1(time_int,resid,w.time_mod);
% plot scatter plot of ocean and resid_mod
figure('position',[scrsz],'color',[1 1 1],'visible','on');
hold on
scatter(resid_mod,oceanf)
% computing correlation coefficient and linear regression model and texting on the figure
[r,p]=corrcoef(resid_mod,oceanf,'rows','complete');
text(0.1,0.9,['r = ',num2str(r(1,2),'%.2f')],'Units','normalized','fontsize',12,'fontweight','bold')
text(0.1,0.8,['p = ',num2str(p(1,2),'%.2f')],'Units','normalized','fontsize',12,'fontweight','bold')
% computing linear regression model
mdl = fitlm(resid_mod,oceanf);
% plotting the linear regression model
plot(mdl)
% displaying the linear regression model on the figure
text(0.1,0.7,['y = ',num2str(mdl.Coefficients.Estimate(1),'%.2f'),' + ',num2str(mdl.Coefficients.Estimate(2),'%.2f'),'x'],'Units','normalized','fontsize',12,'fontweight','bold')
xlabel('Mod resid wlv (m)')
ylabel('Obs resod wlv (m)')
title(replace(file{fobs},'_',' '))


%return

figure('position',[scrsz],'color',[1 1 1],'visible','on');
hold on
plot(time_mod,oceanf-resid_mod,'.k')
datetick('x','dd-mmm-yy','keeplimits')
% compute the rate of increase of the oceanf-resid_mod


%figure('position',[scrsz],'color',[1 1 1],'visible','on');
%hold on
mdlo = fitlm(time_mod,oceanf-resid_mod);
plot(mdlo)
text(0.1,0.9,['r = ',num2str(mdlo.Rsquared.Ordinary,'%.2f')],'Units','normalized','fontsize',12,'fontweight','bold')
text(0.1,0.8,['p = ',num2str(mdlo.Coefficients.pValue(2),'%.2f')],'Units','normalized','fontsize',12,'fontweight','bold')
text(0.1,0.7,['y = ',num2str(mdlo.Coefficients.Estimate(1),'%.2f'),' + ',num2str(mdlo.Coefficients.Estimate(2),'%.2f'),'x'],'Units','normalized','fontsize',12,'fontweight','bold')
xlabel('Time')
ylabel('Ocean - resid mod (m)')
title(replace(file{fobs},'_',' '))




