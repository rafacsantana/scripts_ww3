function buildgrid(X0,Y0,X1,Y1,xlength,ylength,delta,theta,scatter,lnd,fileout,dat)
%
% Script to build mesh for SWAN Format IDIA3 
% **Assumes negative values for depth***
%
% use = buildgrid_swan(X0,Y0,X1,Y1,xlength,ylength,delta,theta scatter,lnd,fileout,dat)
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
%
%% Import scatter data
%
load(scatter);
x=X;
y=Y;
z=el;
z=z*-1;
%
%% Check if coordinates are assending (or not)
%
x_chk=X1-X0; % + = assending, - = not
y_chk=Y1-Y0; % + = not, - = assending
%
% Define the x and y range (and flip if needed)
%
if x_chk<0
    x_rge=X0+xlength:(delta*-1):X0;
    x_rge=flip(x_rge,2); 
else
    x_rge=X0:delta:X0+xlength;
end
%
if y_chk>0
    y_rge=Y0+ylength:(delta*-1):Y0;
    y_rge=flip(y_rge,2); 
else
    y_rge=Y0:delta:Y0+ylength;
end
%
%% Build mesh
%
[xq,yq] = meshgrid(x_rge,y_rge);
%
%% Rotate grid
%
% Define a theta degree counter-clockwise rotation matrix
%
theta=deg2rad(theta);
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
%
% Rotate data row by row
% Define the x- and y-data for the original line we would like to rotate
ref = size(xq);
for i=1:ref(1)
% create a matrix of these points
v = [xq((i),:); yq((i),:)];
% Define origin for basis ofrotation
x_center = X0;
y_center = Y0;
% create a matrix which will be used later in calculations
center = repmat([x_center; y_center], 1, ref(2));
% Do the rotation...
s = v - center;     % shift points in the plane so that the center of rotation is at the origin
so = R*s;           % apply the rotation about the origin
vo = so + center;   % shift again so the origin goes back to the desired center of rotation
%
% Extract rotated grids
%
xq((i),:) = vo((1),:);
yq((i),:) = vo((2),:);
end
%
%% Interpolate to grid (note this extrapolates as well and care is needed if
% you dont have enough data extendign beyound the mesh boundary.  
%
F = scatteredInterpolant(x,y,z,'linear','nearest');
zq = F(xq,yq);
%
% Plot bathymetry
%
surf(xq,yq,zq)
view(2)
axis equal
%%
% Set land points to EXCEPTION
%
ind=find(zq<lnd);
zq(ind)=-999;
%
% Plot grid
%
figure
mesh(xq,yq,zq)
view(2)
axis equal
%%
% Export data to file
%
format='%8.4f ';
filename=fileout;
fid=fopen(filename,'w');
%
% Write mxinp and myinp header
%
ref=size(xq);
fprintf(fid,'%% mxinp');
fprintf(fid, '%8.0f ',ref(2)-1);
fprintf(fid,'\n');
fprintf(fid,'%% myinp');
fprintf(fid, '%8.0f ',ref(1)-1);
fprintf(fid,'\n');
%
% Write bathy grid data
%
for i=1:size(zq,1)
   fprintf(fid,format,zq(i,:));
   fprintf(fid,'\n');
end
fclose(fid);
return
