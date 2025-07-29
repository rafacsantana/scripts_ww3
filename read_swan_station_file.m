clear;
close all

path='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/swan/projects/VANUATU/lv_bnd_4_xbeach_2/';
file=[path,'Locs.out'];

bnd=importdata(file);
time_str=num2str(bnd.data(:,1));

for i=1:size(time_str,1)

  time_stri=replace(time_str(i,:),' ','');

  if length(time_stri)==8
    time(i)=datenum(time_stri,'YYYYmmdd');
  elseif length(time_stri)==10
    time(i)=datenum(time_stri,'YYYYmmdd.HH');
  elseif length(time_stri)==11
    time(i)=datenum(time_stri,'YYYYmmdd.HH');
  elseif length(time_stri)==12
    time(i)=datenum(time_str(i,:),'YYYYmmdd.HHMM');
  end

end

plot(time,time)
datetick('x')

