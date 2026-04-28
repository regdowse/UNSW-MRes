% param_eddy_tracksing  (Call a function ) 
% param_e ddy_tracking sets user defined paths and parameters  
%-----------------------------------  
%Invoke the netcdf  (https://downloads.unidata.ucar.edu/netcdf/ ) and m_map ( https://www.eoas.ubc.ca/
~rich/map.html ) toolkit s 
addpath(genpath( './mexcdf/' )); 
addpath(genpath( './m_map/' ));                                               
% Create data output directory  
Main_dir= '../Eddy_detection/' ;                                              
mkdir(Main_dir);  
% DUACS data is stored in the directory by year  
data_in= './';                                            ¨ 
 
%%%%%%%%%%%%%%%% spatial and temporary %%%%%%%%%%%%%%%%%%%%  
% Detection eddy  area code  name 
postfix= 'CS';                                                             
% Year of beginning of eddy  detection  
yearmin=2019;                     
% Year of end ing of eddy detection  
yearmax=2019;  
% Date  of beginning of eddy detection  
eddy_date_min = datenum( 2019,01,01);      
% Date  of ending of eddy detection  
eddy_date_max = datenum(2020,12,31);      
% Eddy  detection time range  
timenum  = eddy_date_min:1:eddy_date_max;  
% Scope of eddy  detection area  
latmin=8;      
latmax=20;      
lonmin=270;     
lonmax= 300;     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%Sets the number of kernels used for parallel computing  
CorNum=2;  
% a: first parameter for the detection of eddy centers  (number of grid points for v and u velocity reversal;  
a-1 is used to inspe ct the rotation of the vectors around the eddy center)  
a=4; 
%b: second parameter for the detection of eddy centers   (number of grid points to define the area to detect 
velocity minimum)  
b=3; 
% numbuer of grid points to define the initial area to compute  eddy dimensions  
rad=12;  
% distance in  degree  to define the area used to derive eddy tracks  
degree=1;r=1.2;  
% flag to activate options for East -West periodic (e.g. global) domains.  
global  periodic  
periodic=0;  
% the best fit to an ellipse for the given set of points , constrain the eddy  pattern  
ellipse_threshold=0.85;  
 
