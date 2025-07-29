% this script loads different time slices and combines them into a single including nan between periods with data available
clear
close all

path='/scale_wlg_nobackup/filesets/nobackup/niwa03150/WAVE/hindcast/NZWAVE-GFDL-CCAM/matlab/';

files={'Tairua_2005010100_2014123000.mat', 'Tairua_2030010100_2039123100.mat',  'Tairua_2050010100_2069123100.mat',  'Tairua_2080010100_2099123100.mat'};


% loading each file and concatenating the dataset

time_modc=[]; %hs_modc=[]; tp_modc=[]; tm01_modc=[]; tm02_modc=[]; tm_modc=[]; pd_modc=[]; md_modc=[]; ds_modc=[]; wlv_modc=[];
%ucur_modc=[]; %vcur_modc=[]; we_modc=[];
display('Loading files')
for f=1:length(files)
	filename=[path files{f}];
	display(['Loading ' filename])
	load(filename)

  time_modc=[time_modc;time_mod];

end


% creating a hourly time vector from beginning to end of the dataset
% and inserting nans between periods with data available
time_i=time_modc(1):1/24:time_modc(end);
model=nan(length(time_i),12);

disp(['Time vector created from ' datestr(time_i(1)) ' to ' datestr(time_i(end))])
for f=1:length(files)
	filename=[path files{f}];
	display(['Loading ' filename])
	load(filename)

  % finding the beginning and end of the time period and inserting into model
	iloc=find(time_i==time_mod(1));
	floc=find(time_i==time_mod(end));
	model(iloc:floc,1)=tp_mod;
	model(iloc:floc,2)=pd_mod;
	model(iloc:floc,6)=hs_mod;
	model(iloc:floc,4)=tm01_mod;
	model(iloc:floc,5)=tm02_mod;
	model(iloc:floc,3)=ds_mod;
	model(iloc:floc,7)=nan; % md_mod;
	model(iloc:floc,8)=tm_mod;
	model(iloc:floc,9)=wlv_mod; % we_mod
	model(iloc:floc,10)=ucur_mod;
	model(iloc:floc,11)=vcur_mod;
	model(iloc:floc,12)=nan; % we_mod;

end

% plotting the data
figure
for i=1:3
    subplot(3,1,i)
    if i==1; ii=6; elseif i==2; ii=1; elseif i==3; ii=2; end

    plot(datetime(datevec(time_i)),model(:,ii),'.k')
    title(['Column ' num2str(ii)])
end


tp_mod=       model(:,1); 
pd_mod=       model(:,2);
hs_mod=       model(:,6);
tm01_mod=     model(:,4);
tm02_mod=     model(:,5);
ds_mod=       model(:,3);
md_mod=       model(:,7);
tm_mod=       model(:,8);
wlv_mod=      model(:,9);% we_mod
ucur_mod=    model(:,10);
vcur_mod=    model(:,11);


time_mod=time_i;

fda=[filename(1:end-25),datestr(time_i(1),'YYYYmmDDHH'),'_',datestr(time_i(end),'YYYYmmDD'),'00.mat'];
display(['Saving: ',fda]);
save(fda,'time_mod','tp_mod','pd_mod','md_mod','ds_mod','tm01_mod','tm02_mod','hs_mod','tm_mod','wlv_mod','ucur_mod','vcur_mod')% ,'we',model)












