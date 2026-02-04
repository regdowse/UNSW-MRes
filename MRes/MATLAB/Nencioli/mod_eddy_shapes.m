function  mod_eddy_shapes(nc_u,nc_v,nc_dim,rad,a,path_out)  
% mod_eddy_shapes.m  (call function)  
% mod_eddy_shapes(nc_u,nc_v,nc_dim,rad,a,path_out) computes the shapes of  the eddies identified by the 
centers detected with mod_eddy_centers.m  and saves them in  [path_out,'eddy_shapes'];  
% nc_u and nc_v are the time series of the 2D u and v velocity field;  
% nc_dim defines the domain dimensions;  
% rad defines the area where the streamfunction (PSI) is first computed  
% a is one of the parameters from the eddy dete ction constraints  
% path_out is the path where the algorithm output is saved  
% (For a description of the input parameters see param_eddy_tracking.m)  
% Eddy shapes are saved in the structure array [path_out,'eddy_shapes']:  
% shapes(t).lonlat(n) : lonlat is a [2xm double] cell array, containing  the longitude and latitude position of 
the m vertices that define the  boundaries of the n -th eddy of the t -th day;  
% Another file which contain s information on the process to  compute eddy shapes is saved as 
[path_out,' warnings_shapes']:  
% warn_shapes(t).land(n) : 1 if land was found, and the area where the  eddy shape was computed was 
resized;  
% warn_shapes(t).no_curve(n): 1 if no closed contour of PSI was found  around the eddy center. Eddy is 
assumed  circular with radiu s "a-1"; 
% warn_shapes(t).large_curve(n): 1 if no closed contour around the  center with increasing velocity across. 
Shape  is defined by just the larges closed contour;  
% warn_shapes(t).fac(n): number of times the area where the shape is  computed was enlarg ed; 
% (n is the total number of eddy for a given day; t the total number of  days) 
% The function returns also a log file [path_out,'log_eddy_shapes.txt']  which contains additional 
information on the process to compute eddy  shapes.  
% (NOTE: if the file alre ady exist, the new log file will be append to it)  
% 'eddy_dim' is used to compute eddy shapes. Check the documentation in eddy_dim.m for further details.  
%-------------------------  
%   Ver. 1.2 May.2010  
%   Ver. 1.1 Dec.2009  
%   Authors: Francesco Nencioli , francesco.nencioli@opl.ucsb.edu,  
%            Charles Dong, cdong@atmos.ucla.edu  
%-------------------------  
  
% begin the log file  
diary([path_out, 'log_eddy_shapes.txt' ]); 
%----------------------------------------------  
% load eddy centers  
load( [path_out, 'eddy_centers' ]); 
% load time variable  
nc=netcdf(nc_u);  
time=nc{ 'day'}(:); 
close(nc)  
% load coordinates  
nc=netcdf(nc_dim);  
lon=nc{ 'lon'}(:); 
lat=nc{ 'lat'}(:); 
close(nc)  
% load mask  
nc=netcdf(nc_dim);  
mask=nc{ 'mask' }(:); 
close(nc)  
%--------------- -------------------------------  
% Added option for periodic East -West boundaries  
% (first and last column of the domain are identical)  
global  periodic  
if periodic==1  
    % check if minimum and maximum longitude are constant with latitude  
    if ~isempty (find(diff(min(lon,[],2))~=0)) | ... 
            ~isempty(find(diff(max(lon,[],2))~=0))  
        error( 'Problem with your domain: Lonmax and Lonmin vary with Lat!' ) 
    end 
    mlon=min(lon(1,:));  
    Mlon=max(lon(1,:));  
    % shifted longitudinal positions  
    Dlon=lon+(Mlon -mlon);  
    dlon=lon -(Mlon -mlon);  
    % lat and lon are expanded by adding the shifted positions to the  
    % beginning and to the end of the domain;  
    % (they are now 3x their original size: in case you wanna reduce this  
    % make sure to change also u and v expansion in eddy_dim.m)  
    lon=[dlon(:,1:end -1),lon,Dlon(:,2:end)];  
    lat=[lat(:,1:end -1),lat,lat(:,2:end)];  
end 
% coordinates extremes  
Lonmin=min(min(lon));  
Lonmax=max(max(lon));  
Latmin=min(min(lat));  
Latmax=max(max(lat));  
  
% Compute eddy shape --------------------------  
% preallocate shape and warning array  
shapes = repmat(struct( 'lonlat' ,{}),1,length(centers));  
warn_shapes = repmat(struct( 'land' ,{},'no_curve' ,{}, ... 
    'large_curve' ,{},'fac',{}),1,length(centers));  
% loop through all days of the time series  
for i=1:length(centers)  
    disp([ 'Day ' ,num2str(time(i)), ' %------------- ']) 
    % loop through all centers detected for a given day  
    for ii=1:length(centers(i).type)  
        % factor to increase the area where PSI is computed, if eddy  
        % dimesions are too big!!!  
        fac=1; % initially set to no increase  
        % lonlat is the computed shape;  
        % the others are flags output by eddy_dim, and saved in  
        %  'warnings_shapes';  
        % box is a flag that indicates if eddy shape is close to the area  
        %  boundaries (1.5 grid points; limit set in eddy_dim);  
        [lonlat warn land box large]=eddy_dim(nc_u,nc_v, ... 
            centers(i).day, centers(i).i(ii),centers(i).j(ii), ... 
            lon,lat,Lonmin,Lonmax,Latmin,Latmax, ... 
            mask,fac,rad,a);  
        % temporary save eddy_shape  
        tmp_lonlat=lonlat;  
        % log file messages -------------------------  
        if large==1  || warn==1 || land==1 || box==1  
            disp([ ' Eddy ' ,num2str(ii), ':']) 
            if land==1  
                disp( '   Found land!!!' ) 
            end 
            if warn==1  
                disp([ '   Warning: No streamlines ' , ... 
                    'closed around the center!!!!!!' ]) 
            end 
            if large==1  
                disp( '   Largest closed curve' ) 
            end 
            if large==1 || warn==1 || land==1  
                disp( ' ') 
            end 
        end 
        %----------------------------------------------------------  
        %----------------------------------------------------------  
        % area where psi is computed is enalrged while the eddy shape  
        % is close to the area boundaries (1.5  grid points; limit set in  
        % eddy_dim), and no land is found within the area  
        while  land==0 && box==1 && large==0 && warn==0  
            fac=fac+1; % this determine larger area  
            disp([ '   Big eddy: going to fac = ' ,num2str(fac)])  
            % compute eddy shape in the larger area  
            [lonlat warn land box large]=eddy_dim(nc_u,nc_v, ... 
                centers(i).day,centers(i).i(ii),centers(i).j(ii), ... 
                lon,lat,Lonmin,Lonmax,Latmin,Latmax, ... 
                mask,fac,rad,a);  
            % flags -----------------------------------  
            if land==1  
                disp( '     Found land!!!' ) 
            end 
            % if no closed curve in the larger area then final eddy shape  
            % is the on e computed in the smaller area  
            if warn==1 || large==1  
                disp([ '     No closed or largest curve at fac ' , ... 
                    num2str(fac), ' back to curve at fac ' , ... 
                    num2str(fac -1),'!!!']) 
                lonlat=tmp_lonlat;  
            % if eddy shape still close to the area boundaries then  
            % temporary save lonlat  
            elseif  land==0  
                tmp_lonlat=lonlat;  
            end 
            if warn==1 || large==1 || box==0 || land==1  
                disp( ' ') 
            end 
            %-------------------------------------------  
        end 
        %----------------------------------------------------------  
        % This is to fix the flags saved in case going to larger area  
        % resulted in no closed conotur of PSI around the center  
        if (large==1 || warn==1) && fac>1  
            fac=fac -1; 
            large=0;  
            warn=0;  
        end 
        % in case of periodic domain the eddy dimensions are shifted back  
        % to fit within the initial domain boundaries  
        if periodic==1  
            lonlat(1,lonlat(1,:)>Mlon)=lonlat(1,lonlat(1,:)>Mlon) -(Mlon -mlon);  
            lonlat(1,lonlat(1,:)< mlon)=lonlat(1,lonlat(1,:)<mlon)+(Mlon -mlon);  
        end 
        % save dim in a struct array  
        shapes(i).lonlat(ii)={lonlat};  
        % warnings from shape computation  
        warn_shapes(i).land(ii)=land;  
        warn_shapes(i).no_curve(ii)=warn;  
        warn_shapes(i).large_curve(ii)=large;  
        warn_shapes(i).fac(ii)=fac;  
    end 
end 
% output eddy_shapes and shapes warnings  
save([path_out, 'eddy_shapes' ],'shapes' ) 
% save([path_out,'warnings_shapes'],'warn_shapes')  
% close log file  
diary off 
 
