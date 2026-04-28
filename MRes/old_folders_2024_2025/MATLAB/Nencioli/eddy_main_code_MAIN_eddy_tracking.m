% MAIN_eddy_tracks  (General procedure for automatic eddy  recognition ) 
% MAIN_eddy_tracking is the main fun ction of the eddy detection and  tracking . It returns position of the 
% centers, dimensions and  tracks of the eddies detected from th e time series of a 2 -D velocity  field. 
%-------------------------  
%   Ver. 1.3 May.2022  
%   Ver. 1.2 May.2010  
%   Ver. 1.1 Dec.2009  
%   Authors: Francesco Nencioli, francesco.nencioli@opl.ucsb.edu,  
%            Charles Dong, cdong@atmos.ucla.edu  
%            Yu Liu, yliu@nuist.edu.cn    update 2015 -01-24 
%            Lingxiao Liu, yx@nuist,edu.cn update 202 2-05-01 
%-------------------------  
  
clear;clc  
warning off 
% Set the eddy automatic detection algorithm related directories and parameters , users should modifi y 
% param_edd y_tracking.m according to their  settings  
param_eddy_tracking  
% Based on the start and end time of eddy detection, cyclic detection is carried out year by year  
for yy=yearmin:yearmax  
    disp([ 'start.........' ,num2str(yy)]);  
    % Set the time dimension for a whole year and select the detection time  
    tin=find(timenum>=datenum(yy,1,1) & timenum<=datenum(yy+1,1,1) -1); 
    % From the eddy detection time selected above, set the start time and end time  
    starttime=timenum(tin(1));  
    endtime=timenum(tin(end));  
    % Set a general catalog of eddy data output in years  
uv_file_in=[data_in,num2str(yy), '/']; 
% Set a eddy background field data directory under the general directory  
    data_out=[Main_dir, 'Eddy_' ,postfix, '/Data_' ,num2str( yy),'/']; 
path_in=data_out;  
% Set eddy center s, boundaries  and tracks data directory under the general directory  
path_out=[Main_dir, 'Eddy_' ,postfix, '/Tracks_' ,num2str(yy), '/']; 
% If the eddy detection time spans years, create a total tracking directory  
    path_out_all=[Main_dir, 'Eddy_' ,postfix, '/Tracks_all' ,'/']; 
% Create the above directory  
    mkdir(data_out);  
    mkdir(path_out);  
    mkdir(path_out_all);  
    % Read the abnormal data of geostrophic velocity  for quality control and interpolation , and save  it as 
background field data  
MAIN_create_surf_fields_UV_AVISO_grided(uv_file_in,data_out,starttime,endtime,latmin,latmax,lonmin,lonmax);  
% Read the background field of geostrophic velocity  anomaly and the position of latitude and 
longitude grid data  
% nc_u : full name of the netcdf file with the zonal component of  velocity (ssu) and the time index (day) 
% nc_v: full name of the netcdf file with the meridional component of  velocity (ssv) and the time 
index (day)  
% nc_dim: full name of the netcdf file with the  domain coordinates   (longitude and latitude) and the 
velocity mask (land -points=0;  ocean -points=1)  
    nc_u=[path_in, 'ssu.nc' ]; 
    nc_v=[path_in, 'ssv.nc' ]; 
    nc_dim=[path_in, 'lon_lat.nc' ]; 
    % Identify the eddy center  
    mod_eddy_centers( nc_u,nc_v,nc_dim,a,b,path_out)  
    % Calculate eddy boundaries  
    mod_eddy_shapes(nc_u,nc_v,nc_dim,rad,a,path_out)  
     
end 
% Integrate eddy centers and boundaries in the same or different years for cross -year tracking  
full_time_tracking(postfix,yearmin,yearmax,Main_dir,path_out_all,degree,r)
 
