clear 
close all


filenc='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/swan/projects/VANUATU/data/gebco_2023_n-12.0_s-19.0_w163.0_e170.0.nc';
el=double(ncread(filenc,'elevation')');
X=ncread(filenc,'lon');
Y=ncread(filenc,'lat');

% rectangle that contains the three islands not shown in gebco data
loni=167.5;lonf=168.3;lati=-16;latf=-14.6;
ilon=find(X>loni & X<lonf);
ilat=find(Y>lati & Y<latf);
elev=el(ilat,ilon);
elev(elev>-700)=10;
el(ilat,ilon)=elev;

%masking the southernmost island
lati=-16.013835;loni=168.137064; latf=-15.450488; lonf=168.224314;
ilon=find(X>loni & X<lonf);
ilat=find(Y>lati & Y<latf);
elev=el(ilat,ilon);
elev=10;
el(ilat,ilon)=elev;

figure; hold on
pcolor(X,Y,el); 
shading flat; colorbar
contour(X,Y,el,[0 0],'k'); 
contour(X,Y,el,[-700 0],'r'); 

%return

[X,Y]=meshgrid(X,Y); X=X(:); Y=Y(:);
el=el(:);

filein='/scale_wlg_persistent/filesets/project/niwa03150/santanarc/swan/projects/VANUATU/data/bathy_islands.mat';
save(filein,'X','Y','el')


fileout='~/waves/swan/projects/VANUATU/lv_tcwindgen_tc_harold_cfac/depth_islands.dat';
if exist(fileout)==2
   system(['rm -rf ',fileout]);
end


buildgrid(164.504,-17.9958,169.504,-12.9958,5,5,0.04,0,filein,0,fileout,0);

%
% X0 = X origin (bottom left)
% Y0 = Y origin (botton left)
% X1 = X coord (top right) - to assess if coordinates assending (or not)
% Y1 = Y coord (top right)- to assess if coordinates are assending (or not)
% xlength = length of grid in x direction (m or degrees)
% ylength = length of grid in y direction (m or degrees)
% delta = grid resolution (m or degrees) - always positive
% theta = grid rotation around XO,YO in degrees. +angle is counter clockwise, -angle is clockwise
% scatter = xyz bathymetry grid data 'filename'
% lnd = Land value - anything greater than lnd is set to EXCEPTION -999
%       this should be a positive number
% fileout = 'filename' of output file (IDLA3 format)
% dat = datum shift
% Created 31/05/2022 © Climatise Ltd

 
