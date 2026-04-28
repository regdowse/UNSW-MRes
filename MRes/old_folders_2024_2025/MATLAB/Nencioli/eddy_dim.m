function  [lonlat warn land box large]=eddy_dim(nc_u,nc_v, ... 
    day,C_I,C_J,lon,lat,Lonmin,Lonmax,Latmin,Latmax, ... 
    mask,fac,rad,a)  
% eddy_dim.m  (call function)  
% eddy_dim(nc_u,nc_v,day,C_I,C_J,lon,lat,Lonmin,Lonmax,Latmin,Latmax,mask, fac,rad,a) computes the 
shape of the eddy defined by the center C_I,C_J.  
% 
% nc_u and nc_v are the time series of the 2D u and v velocity field;  
% day is the day of the 2D velocity field to analyze ; 
% C_I and C_J are the indices of the eddy center  
% lon and lat are l ongitudes and latitudes of points in the domain;  
% Lonmin, Lonmax, Latmin and Latmax are the extremes of the domain;  
% mask is the matrix that defines sea (1) and land points (0);  
% fac is the factor to enlarge the area where the shape is computed;  
% rad is used to define the area where the shape is computed;  
% a is one of the parameters from the eddy detection; is used here to  define eddy radius in case no closed 
contour of PSI is found;  
% 
% OUTPUT:  
% lonlat  is the array containing longitudes (first row) and latitudes(second row) of the vertices defining the 
eddy shape;  
% warn, land, box and large are flags that give some informations on the procedure used to compute the 
eddy shape. See 'mod_eddy_shapes.m' fo r further details.  
% 
% 'compute_psi' is used to compute the streamfunction field integrating u -and v -component f velocity. 
Check the documentation in 'compute_psi.m' for further details.  
% 
% 'max_curve' is used to compute eddy shape (defined as the largest  closed contour of PSI around the eddy 
center across which velocity magnitude increases). Check the documentation in 'max_curve.m' for further 
details.  
% 
%-------------------------  
%   Ver. 1.2 May.2010  
%   Ver. 1.1 Dec.2009  
%   Authors: Francesco Nencioli , francesco.nencioli@opl.ucsb.edu,  
%            Charles Dong, cdong@atmos.ucla.edu  
%-------------------------  
% increase the dimensions of the area where eddy shape is computed  
rad=rad*fac;  
% initialized the flags  
land=0;  
warn=0;  
box=0;  
% find day index of  the surface field to analyze -----------  
nc=netcdf(nc_u);  
time=nc{ 'day'}(:); 
close(nc)  
iday=find(time==day);  
% load 2D velocity field -----------------------------------  
nc=netcdf(nc_u);  
u=nc{ 'ssu'}(iday,:,:);  
close(nc)  
nc=netcdf(nc_v);  
v=nc{ 'ssv'}(iday,: ,:); 
close(nc)  
%-----------------------------------------------------------  
% mask velocity data  
u(mask==0)=NaN;  
v(mask==0)=NaN;  
% Added option for periodic East -West boundaries  
% (first and last column of the domain are identical)  
global  periodic  
if periodic==1  
    % size of old/smaller domain  
    dI=size(u(:,1:end -1),2); 
    % u and v are expanded by adding the the whole domain before and after  
    % (they are now 3x their original size: in case you wanna reduce this  
    % make sure to change also lat and lon expansion in mod_eddy_shapes.m)  
    u=[u(:,1:end -1),u,u(:,2:end)];  
    v=[v(:,1:end -1),v,v(:,2:end)];  
    % adjust C_I to larger domain  
    C_I=C_I+dI;  
end 
  
% center coordinates (C_J and C_I are center indices in the whole domain)  
c_lat=lat(C _J,C_I);  
c_lon=lon(C_J,C_I);  
  
% resize coordinate and velocity matrix  
% (making sure not to go outside the domain)  
lat=lat(max(C_J -rad,1):min(C_J+rad,size(lat,1)), ... 
    max(C_I -rad,1):min(C_I+rad,size(lat,2)));  
lon=lon(max(C_J -rad,1):min(C_J+rad,size( lon,1)), ... 
    max(C_I -rad,1):min(C_I+rad,size(lon,2)));  
v=v(max(C_J -rad,1):min(C_J+rad,size(v,1)), ... 
    max(C_I -rad,1):min(C_I+rad,size(v,2)));  
u=u(max(C_J -rad,1):min(C_J+rad,size(u,1)), ... 
    max(C_I -rad,1):min(C_I+rad,size(u,2)));  
  
% indices of the eddy center in the smaller area  
[c_j c_i]=find(lat==c_lat & lon==c_lon);  
  
% inspect if there are land -points in the area:  
% if land -points are present, the area is resized to have only ocean -points  
% (rad becomes the shortes distance from the center t o land)  
[yl xl]=find(isnan(u));  
if ~isempty(xl)  
    land=1; % update flag  
    % compute new 'rad' --------------------------------------  
    % distance in grid points from the center to the land points  
    dpts=sqrt((yl -c_j).^2+(xl -c_i).^2);  
    % problem if the land point is a vertex of the region  
    if floor(min(dpts)/sqrt(2))==min(dpts)/sqrt(2)  
        rad2=floor(min(dpts)/sqrt(2)) -1; 
    else 
        rad2=floor(min(dpts)/sqrt(2));  
    end 
    % rad2 cannot be less than a -1 (bc of fourth constraint)  
    if rad2<a -1 
        rad2=a -1; 
    end 
    % update 'rad' (rad2 has to be bigger than rad at previous fac)  
    if rad2>rad/fac*(fac -1) 
        rad=rad2;  
    else 
        rad=rad/fac*(fac -1); 
    end 
    %---------------------------------------------------------  
    % resize coordinates and velocities  
    lat=lat(max(c_j -rad,1):min(c_j+rad,size(lat,1)), ... 
        max(c_i -rad,1):min(c_i+rad,size(lat,2)));  
    lon=lon(max(c_j -rad,1):min(c_j+rad,size (lon,1)), ... 
        max(c_i -rad,1):min(c_i+rad,size(lon,2)));     
    v=v(max(c_j -rad,1):min(c_j+rad,size(v,1)), ... 
        max(c_i -rad,1):min(c_i+rad,size(v,2)));  
    u=u(max(c_j -rad,1):min(c_j+rad,size(u,1)), ... 
        max(c_i -rad,1):min(c_i+rad,siz e(u,2)));  
    % indices of the eddy center in the new smaller area  
    [c_j c_i]=find(lat==c_lat & lon==c_lon);  
end 
% compute velocity magnitude  
vel=sqrt(u.^2+v.^2);  
  
% convert lon and lat into km distance matrices  
% (needed to compute the streamfunction  field) 
for i=1:size(lon,1)  
    km_i(i,:)=[0; cumsum(m_lldist(lon(i,:),lat(i,:)))];  
end 
for i=1:size(lat,2)  
    km_j(:,i)=[0; cumsum(m_lldist(lon(:,i),lat(:,i)))];  
end 
  
% compute grid km spacing -----------  
% (used for psi computation and  
%  for conversi on of eddy boundary from km to lat lon)  
% 
% conversion factor from km to deg lat (equal lat spacing)  
km_dlat=m_lldist([lon(1,1) lon(1,1)],[lat(1,1) lat(2,1)]);  
km2lat=(lat(2,1) -lat(1,1))/km_dlat;  
% conversion factor from km to deg lon (varying lon  spacing with lat)  
for i=1:size(lon,1)  
    km_dlon(i,1)=m_lldist([lon(i,1) lon(i,2)],[lat(i,1) lat(i,2)]);  
end 
dlon=diff(lon,1,2);  
dlon=dlon(:,1);  
km2lon=dlon./km_dlon;  
%------------------------------------  
  
% compute streamfunction field  
psi=compute_psi( u,v,km_dlon,km_dlat);  
% eddy center position in the small area (kilometric coordinates)  
km_cj=km_j(c_j,c_i);  
km_ci=km_i(c_j,c_i);  
% compute eddy shape  
[eddy_lim large]=max_curve(km_i,km_j,psi,km_ci,km_cj,vel);  
  
% in case there is no streamfunction curve c losed around the center the  
% eddy dimension is set to a -1 points radius  
if isempty(eddy_lim)  
    warn=1; % update flag  
    R=min([km_i(1,1+(a -1))-km_i(1,1) km_i(end,1+(a -1))-km_i(end,1) ... 
        km_j(1+(a -1),1) -km_j(1,1) km_j(end,1) -km_j(end -(a-1),1)]) ; 
    theta = 0:pi/20:2*pi;  
    X = (R*cos(theta)+km_ci);  
    Y = (R*sin(theta)+km_cj);  
    eddy_lim=[X;Y];  
    % in case eddy shape is too close to the area boundary, R is further  
    % reduced to 1  
    if ~isempty(find(eddy_lim<0 , 1))  
        R=min([ km_i(1,2) -km_i(1,1) km_i(end,2) -km_i(end,1) ... 
        km_j(2,1) -km_j(1,1) km_j(end,1) -km_j(end -1,1)]);  
        theta = 0:pi/20:2*pi;  
        X = (R*cos(theta)+km_ci);  
        Y = (R*sin(theta)+km_cj);  
        eddy_lim=[X;Y];  
    end 
end 
  
% parameter tha t defines if the computed eddy shape is too close to the  
% area boundaries.  
lim_bound=1.5 * ...      % distance is set in km  
    m_lldist([lon(c_j,1) lon(c_j,2)],[lat(c_j,1) lat(c_j,1)]);  
% if the eddy shape is less than 'lim_bound' from teh area boundar ies, and  
% the area boundaries does not coincides with the domain boundaries, then  
% the flag box is set to 1 so that eddy shape will be computed in a larger  
% area (see also mod_eddy_shapes.m)  
if (min(eddy_lim(1,:)) -min(min(km_i)) < lim_bound || ... 
        max(max(km_i)) -max(eddy_lim(1,:)) < lim_bound || ... 
        min(eddy_lim(2,:)) -min(min(km_j)) < lim_bound || ... 
        max(max(km_j)) -max(eddy_lim(2,:)) < lim_bound) && ... 
        (min(min(lon))~=Lonmin && ... 
        max(max(lon))~=Lonmax && ... 
        min(min(lat))~=Latmin && ... 
        max(max(lat))~=Latmax)  
    box=1;    
end 
  
% convert eddy_lim from km to geographic coordinates -----------  
% initialize lonlat  
lonlat=zeros(size(eddy_lim));  
% convert point by point  
for j=1:length(eddy_lim (1,:)) 
  
    % 1) lat (constant km spacing)  
    diff_j=eddy_lim(2,j) -km_j(:,1);  
    dkm_j=min(diff_j(diff_j>=0));  
    lat_j=find(diff_j==dkm_j);  
    lonlat(2,j)=lat(lat_j,1)+dkm_j*km2lat;  
    % 2) lon (km spacing varies with lat)  
    diff_i=eddy_lim(1,j) -km_i(lat_j,:);  
    dkm_i=min(diff_i(diff_i>=0));  
    lon_i=find(diff_i==dkm_i);  
    lonlat(1,j)=lon(lat_j,lon_i)+dkm_i*km2lon(lat_j);  
end 
 
