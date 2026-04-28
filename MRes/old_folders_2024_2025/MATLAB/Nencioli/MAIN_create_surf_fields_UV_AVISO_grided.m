function  MAIN_create_surf_fields_UV_AVISO_grided(uv_file_in,data_out, ... 
       starttime,endtime,latmin,latmax,lonmin,lonmax)  
% MAIN_create_surf_fields_UV_AVISO_grided (call function)  
% Read the abnormal data of geospin  speed for quality control and interpolation , and save it as background 
field data  
% uv_file_in is the directory where the geostrophic field anomaly data is located;  
% data_out is eddy detection area background field output directory;  
% starttime is eddy star t detection time  
% endtime is eddy end detection time  
% latmin,latmax,lonmin,lonmax are eddy detection area range  
 
% Detection area background field output variable naming  
out_name={'ssu','ssv','ssh'};  
% eddy detection time length  
time=starttime:1:endtime;    
% Read the abnormal data list of abnormal geostrophic velocity data  
file=dir([uv_file_in,'*.nc']);  
% Read the file time in the abnormal geostrophic velocity data  list 
for i=1:length(file)  
    fname=file(i).name(isstrprop(file(i).name,'digit'));  
    fname=fname(2:9) ; 
    yy=str2num(fname(1:4));  
    mm=str2num(fname(5:6));  
    dd=str2num(fname(7:8));  
    clear fname  
    ftime(i)=datenum(yy,mm,dd);  
    clear yy mm dd  
end 
% Select the eddy detection area from the abnormal geostrophic velocity data  
ngrid=[ uv_file_in,file(1).name];  
lat=ncread(ngrid,'latitude');  
lon=ncread(ngrid,'longitude');  
lon(lon<0)=lon(lon<0)+360;  
Nlat=find(lat<=latmax&lat>=latmin);  
ymin=min(Nlat);ymax=max(Nlat);  
Nlon=find(lon<=lonmax&lon>=lonmin);  
xmin=min(Nlon);xmax=max(Nlon);  
lat=lat (Nlat);  
lon=lon(Nlon);  
[Lon,Lat]=meshgrid(lon,lat);  
% load mask to mask values of u, v in the rho grid  
u=ncread(ngrid,'ugosa');  
u=squeeze(u(Nlon,Nlat))';  
v=ncread(ngrid,'vgosa');  
v=squeeze(v(Nlon,Nlat))';  
vel=sqrt(u.^2+v.^2);  
 
% delete([path,'/',nam e]);  
%%%%% %use etopo5 to create mask %%%%%%  
% Use etopo5 data to make study area mask file  
nc = netcdf( 'etopo5.nc' , 'nowrite' ); 
tlon=nc{ 'topo_lon' }(:); 
tlon(tlon<0)=tlon(tlon<0)+360;  
tlat=nc{ 'topo_lat' }(:); 
topo=nc{ 'topo' }(:); 
Ntlat=find(tlat<=latmax& tlat>=latmin);  
Ntlon=find(tlon<=lonmax&tlon>=lonmin);  
tlat=tlat(Ntlat);  
tlon=tlon(Ntlon);  
topo=topo(Ntlat,Ntlon);  
[TX TY]=meshgrid(tlon,tlat);  
topo=interp2(TX,TY,topo,Lon,Lat);  
topo(topo>=0)=0;  
topo(topo<0)=1;  
mask_rho=topo;  
mask_rho(vel>300)=0;   
% mask created using flagged values of u  
% interpolate lon lat and mask to a higher res grid  
[X Y]=meshgrid(1:size(Lon,2),1:size(Lat,1));  
[Xx Yy]=meshgrid(1:.67:size(Lon,2),1:.67:size(Lat,1));  
Lon=interp2(X,Y,Lon,Xx,Yy);  
Lat=interp2(X,Y,Lat,Xx,Yy);  
mask_rho=inter p2(X,Y,mask_rho,Xx,Yy);  
mask_rho(mask_rho<1)=0;  
% create netcdf file with lat lon  
nw = netcdf([data_out, 'lon_lat.nc' ],'clobber' ); 
% Create dimensions  
nw('x') = size(Lon,2);  
nw('y') = size(Lat,1);  
% Create variables and attributes  
nw{'lon'} = ncdouble( 'y', 'x'); 
nw{'lon'}.long_name = ncchar( 'longitude of RHO -points' ); 
nw{'lon'}.long_name = 'longitude of RHO -points' ; 
nw{'lon'}.units = ncchar( 'degree_east' ); 
nw{'lon'}.units = 'degree_east' ; 
  
nw{'lat'} = ncdouble( 'y', 'x'); 
nw{'lat'}.long_name = ncchar( 'latitude of RHO -points' ); 
nw{'lat'}.long_name = 'latitude of RHO -points' ; 
nw{'lat'}.units = ncchar( 'degree_north' ); 
nw{'lat'}.units = 'degree_north' ; 
  
nw{'mask' } = ncdouble( 'y', 'x'); 
nw{'mask' }.long_name = ncchar( 'mask on RHO -points' ); 
nw{'mask' }.long_name = 'mask on RHO -points' ; 
nw{'mask' }.option_0 = ncchar( 'land' ); 
nw{'mask' }.option_0 = 'land' ; 
nw{'mask' }.option_1 = ncchar( 'water' ); 
nw{'mask' }.option_1 = 'water' ; 
close(nw)  
% create netcdf file for surface fields  
for i=1:size(out_nam e,2) 
    nw = netcdf([data_out,char(out_name(i)), '.nc'], 'clobber' ); 
    %Create dimensions  
    nw('x') = size(Lon,2);  
    nw('y') = size(Lat,1);  
    nw('time' )= size(time,2);  
    %Create variables and attributes  
    nw{'day'} = ncdouble( 'time' ); 
    nw{'day'}.long_name = ncchar( 'day'); 
    nw{'day'}.long_name = 'day'; 
    nw{'day'}.units = ncchar( 'days' ); 
    nw{'day'}.units = 'days' ; 
    nw{char(out_name(i))} = ncdouble( 'time' ,'y', 'x'); 
    close(nw)  
end 
%------------------------------------------------------------  
% fill the netcdf file --------------------------------------  
nw=netcdf([data_out, 'lon_lat.nc' ],'write' ); 
nw{'lat'}(:)=Lat;  
nw{'lon'}(:)=Lon;  
nw{'mask' }(:)=mask_rho;  
close(nw)  
clear lat lon pm pn mask_rho  
% fill the files with surface field  
disp( 'Create input files!' ) 
nw1=netcdf([data_out, 'ssu.nc' ], 'write' ); 
nw2=netcdf([data_out, 'ssv.nc' ], 'write' ); 
nw3=netcdf([data_out, 'ssh.nc' ], 'write' ); 
% Read the abnormal geostrophic velocity data  in the study area and store it in netcdf format  
for i=1:length(time)  
    disp([num2str(i/length(time)*100), '%']) 
    t=time(i);  
    index_t=find(t==ftime);  
    filename=[uv_file_in,file(index_t).name];  
    clear index_t  
    nc = filename;  
     
    ssu=ncread(filename, 'ugosa' );   
    ssu=ssu(xmin:xmax,ymin:ymax)';  
    ssu(ssu<= -2e+9)=nan;  
    ssu=interp2(X,Y,ssu,Xx,Yy);  
    nw1{ 'day'}(i)=time(i);   
    nw1{ 'ssu'}(i,:,:)=ssu;  
    clear ssu 
     
    ssv=ncread(filename, 'vgosa' ); 
    ssv=ssv( xmin:xmax,ymin:ymax)';  
    ssv(ssv<= -2e+9)=nan;  
    ssv=interp2(X,Y,ssv,Xx,Yy);  
    nw2{ 'day'}(i)=time(i);   
    nw2{ 'ssv'}(i,:,:)=ssv;  
    clear ssv 
     
    ssh=ncread(filename, 'sla'); 
    ssh=ssh(xmin:xmax,ymin:ymax)';  
    ssh(ssh<= -2e+9)=nan;  
    ssh=i nterp2(X,Y,ssh,Xx,Yy);  
    nw3{ 'day'}(i)=time(i);   
    nw3{ 'ssh'}(i,:,:)=ssh;  
    clear ssh 
end 
 
close(nw1)  
close(nw2)  
close(nw3)  
  
  
  
 
 
