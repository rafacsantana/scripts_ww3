clear

close all

path='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/data/bathy/nzbathy_2016/';

file=[path,'nzbathy_2016.tif'];

t=Tiff(file);

imageData = read(t);

imshow(imageData);

title('NZ 2016 bathhy (m)')

colorbar




