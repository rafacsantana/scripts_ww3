clear
close all


path='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/';
path_mod=[path,'model/'];
path=[path,'data/obs/'];


proc_data=1;
save_fig =0;


file='SteepHead_SP_FullRecord_QC';
%    {'date'}    {'PeakPerid'}    {'PeakDirection'}    {'Directional spread'}    {'Tm01'}    {'Tm02'}    {'Hs'}    {'Qp'}

if proc_data==1
  
  display(['Processing: ',path,file,'.txt']);
  ob=importdata([path,file,'.txt']);
  ob.textdata(1,:)
  for i=2:length(ob.textdata)
    time_obs(i-1)=datenum(ob.textdata{i,1},'YYYY-mm-ddTHH:MM:SS');
  end
  obs=ob.data;
  lat_obs=-43+(45/60); lon_obs=173+(20/60);
  save([path,file,'.mat'],'time_obs','obs','lat_obs','lon_obs')

else

  display(['Loading: ',path,file,'.mat']);
  load([path,file,'.mat'])%,'time','obs')
  lat_obs=-43+(45/60); lon_obs=173+(20/60);

end

scrsz=[1 1 1366 768];
%scrsz=get(0,'screensize');
figure('position',scrsz,'color',[1 1 1],'visible','on')
hold on
set(gca,'fontsize',12,'fontweight','bold')

dcol=[1,4,5]; %[6,1,2]; % Hs, Tp, dir
tnames={'PeakPerid',    'PeakDirection',    'Directional spread',   'Tm01','Tm02','Hs','Qp'};
for i=1:3

  subplot(3,1,i);
  set(gca,'fontsize',12,'fontweight','bold')
  hold on
  plot(time_obs,obs(:,dcol(i))','.k','linewidth',2)
  %if i==3
    datetick('x')
  %end
  title(tnames{dcol(i)})

end


if save_fig==1
  export_fig(gcf,'steephead_period','-png','-r150');
end


