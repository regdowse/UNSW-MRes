function  mod_eddy_centers(nc_u,nc_v,nc_dim,a,b,path_out)  
% mod_eddy_centers.m  (call function)  
% mode_eddy_centers(nc_u,nc_v,nc_dim,a,b,pat h_out) detect the eddy centers present in the domain, 
defined by nc_dim , for each time step of the time series of the 2 -D velocity field defined by nc_u and nc_v, 
using the parameters a and b.  
% For a description of the input parameters see param_eddy_tracking.m.  
% Eddy centers are saved in the structure array [path_out,'eddy _centers']:  
% centers(t).day : day when the eddy was detected  
% centers(t).type(n) : eddy type (1 => cyclonic; -1 => anticyclonic)  
% centers(t).lat(n) : eddy center latitude  
% centers(t).lon(n) : eddy center longitude  
% centers(t).i(n) : eddy center column  index 
% centers(t).j(n) : eddy center row index  
% (t is the time index; n is the number of eddies detected at t)  
% 'uv_search' is used to determine the points in the domain that satisfy all 4 constraints. Check the 
documentation in uv_search.m for further  details.  
%-------------------------  
%   Ver. 1.2 May.2010  
%   Ver. 1.1 Dec.2009  
%   Authors: Francesco Nencioli, francesco.nencioli@opl.ucsb.edu,  
%            Charles Dong, cdong@atmos.ucla.edu  
%-------------------------  
  
%-------  load time,  coordinates and mask -------  
nc=netcdf(nc_dim);  
Lat=nc{ 'lat'}(:); 
Lon=nc{ 'lon'}(:); 
mask=nc{ 'mask' }(:); 
close(nc)  
  
nc=netcdf(nc_u);  
time=nc{ 'day'}(:); 
close(nc)  
%---------------------------------------------  
% preallocate centers array  
centers = repmat( struct( 'day',{},'type' ,{},'lat',{},'lon',{}, ... 
    'j',{},'i',{}),1, 10);  
  
% cycle through time steps  
for i=1:length(time)  
    disp([ 'Searching day ' ,num2str(time(i)), ' ... ']) 
    % eddy centers for a given day  
    [eddy_fld_uv eddy_fld_c eddy_fld]=uv_ search(nc_u,nc_v,Lon,Lat,mask, ... 
        i,a,b);  
    % save eddy positions in struct array  
    centers(i).day=time(i);  
    if ~isempty(eddy_fld)  
        for ii=1:length(eddy_fld(:,1))  
            centers(i).type(ii)=eddy_fld(ii,3);  
            centers(i) .lat(ii)=eddy_fld(ii,1);  
            centers(i).lon(ii)=eddy_fld(ii,2);  
            [centers(i).j(ii) centers(i).i(ii)]=find(Lon==eddy_fld(ii,2) ... 
                & Lat==eddy_fld(ii,1));  
        end 
    end 
end 
  
save([path_out, 'eddy_centers' ],'centers' ) 
 
