function  [in_out_vel]=interpolate_limits(pt_i,pt_j,km_i,km_j,dkm,vel,dir)  
% interpolate_limits.m  (call function ) 
% interpolate_limits(pt_i,pt_j,km_i,km_j,dkm,vel,dir) interpolate velocity  magnitudes across the four 
extremes of a closed conotur.  
% 
% pt_i  and pt_j are the location in kilometric coordinate of the curve  extreme under consideration;  
% km_i and km_j are kilometric coordinate matrices (coordinate origin is the south -east corner of the area);  
% dkm is the distance to and from teh extreme where v elocity magnitude  is interpolated;  
% vel is the velocity magnitude field, used to check the increase in  velocity across the closed contour.  
% dir is a string that identify the type of extreme (North 'N', East  'E', South 'S' or West 'W');  
% 
% OUTPUT:  
% in_o u_vel contains the velocity magnitude 'dkm' kilometers before and  after the curve extreme; the two 
are compared in 'max_curve.m'.  
% 
%-------------------------  
%   Ver. 1.2 May.2010  
%   Ver. 1.1 Dec.2009  
%   Authors: Francesco Nencioli, francesco.nencioli@o pl.ucsb.edu,  
%            Charles Dong, cdong@atmos.ucla.edu  
%-------------------------  
  
% include only 'n' points around that point to interpolate vel  
% (griddata is a very slow process) --------------------------------  
% find the closest grid point to t he curve extreme  
dist=sqrt((km_i -pt_i).^2+(km_j -pt_j).^2);  
[d_j d_i]=find(dist==min(min(dist)));  
n=4; % at least 4 points away to avoid a qhull precision warning!!!  
% resize coordinate and velocity matrices  
svel=vel(max(d_j -n,1):min(d_j+n,size(vel,1)), ... 
    max(d_i -n,1):min(d_i+n,size(vel,2)));  
skm_i=km_i(max(d_j -n,1):min(d_j+n,size(vel,1)), ... 
    max(d_i -n,1):min(d_i+n,size(vel,2)));  
skm_j=km_j(max(d_j -n,1):min(d_j+n,size(vel,1)), ... 
    max(d_i -n,1):min(d_i+n,size(vel,2)));  
%--------------------- ----------------------------------------------  
% interpolate vel across the curve along different directions depending on  
% different extremes  
switch  (dir) 
    case 'N' 
        pts=[pt_i,pt_j -dkm;pt_i,pt_j+dkm];  
    case 'S' 
        pts=[pt_i,pt_j+dkm;pt_ i,pt_j -dkm];  
    case 'W' 
        pts=[pt_i+dkm,pt_j;pt_i -dkm,pt_j];  
    case 'E' 
        pts=[pt_i -dkm,pt_j;pt_i+dkm,pt_j];  
end 
% interpolate vel  
in_out_vel=griddata(skm_i,skm_j,svel,pts(:,1),pts(:,2));  
 
