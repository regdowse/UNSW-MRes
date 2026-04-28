function  mod_eddy_tracks(nc_dim,degree,r,path_out)  
% mod_eddy_tracks.m  (call function)  
% mod_eddy_tracks(nc_dim,r,path_out) computes eddy tracks from the eddies  identified in 
'mod_edd_centers.m', and saves them in the file  [path_out,'eddy_tracks'];  
% Eddy  tracks are connected by comparing the detected eddy fields at  sucessive time steps. An eddy at 
time 't+1' is assumed to be the  new position of an eddy of the same type detected at 't', if the two centers are 
found within an area of 'r' km centered around the center position at time 't'.  
% If no eddies are found within the area at 't+1', a second search is  performed at 't+2' within an enlarged 
area. 
% If no eddies are found at 't+2' the eddy is considered dissipated and  the track closed.  
% Tracks with lengt h shorter than 'cut_off' days are not recorded.  
% ('cut_off' can be modified at the beginning of this function).  
% In case two or more eddies are found within the same area, the track is  created connecting the center at 
time 't+1' closest to the center at  time 't'.  
% nc_dim defines the domain dimensions;  
% r (in km) defines the area used to search for eddy centers at successive timesteps;  
% path_out is the path where the algorithm output is saved;  
% (For a description of the input parameters see param_eddy _tracking.m)  
% Eddy tracks are saved in the structure array [path_out,'eddy_tracks']:  
% ('n' is the total number of track recorded; 'm' is the lenght of a given  track)  
% tracks(n).day(m) : days when the eddy centers forming the track where  detected;  
% trac ks(n).type(m) : type of eddy;  
% tracks(n).lat(m) : latitude of the eddy centers forming the track;  
% tracks(n).lon(m) : longitude of the eddy centers forming the track;  
% tracks(n).j(m) : meridional index of the eddy centers forming the track;  
% tracks(n). i(m) : zonal index of the eddy centers forming the track;  
% tracks(n).shapes{m} : shapes of the eddies forming the track; (each  cell is a 2xL array containing the 
longitudinal (first row) and latitudinal (second row) positions of the 'L' vertices defining  the eddy shape for 
the 'm' -th day.  
% [path_out,'removed_tracks'] contains the number of tracks removed  because shorter than 'cut_off' days.  
% The function returns also a log file [path_out,'log_eddy_tracks.txt']  which contains warnings relative to 
the pre sence of multiple eddies within the searching area used to determine the tracks.  
% (NOTE: if the file already exist, the new log file will be append to it)  
% The same information is also saved in [path_out,'warnings_tracks'], an Rx3 array (R is the total n umber 
of warnings), where each column  represent in the order:  day the warning occurred; type of area used (1 for 
area at t+1, 2 for area at t+2); number of eddies detected within the area;  
%-------------------------  
%   Ver. 1.2 May.2010  
%   Ver. 1.1 Dec.2009  
%   Authors: Francesco Nencioli, francesco.nencioli@opl.ucsb.edu,  
%            Charles Dong, cdong@atmos.ucla.edu  
%-------------------------  
  
% define the length (in days) below which a track is not recorded  
cut_off=1;  
% begin the log file  
diary( [path_out, 'log_eddy_tracks.txt' ]); 
%---------------------------------------------  
% load eddy centers and eddy shapes  
load([path_out, 'eddy_centers' ]); 
load([path_out, 'eddy_shapes' ]); 
%load coordinates  
nc=netcdf(nc_dim);  
Lat=nc{ 'lat'}(:); 
lat=Lat(:,1);  
Lon= nc{'lon'}(:); 
close(nc)  
%---------------------------------------------  
  
% intitialize search and tracks structure arrays  
% - search contains all the open tracks;  
% - tracks contains all the closed tracks which will be saved;  
% - tracks contains all the closed tracks which will be saved;  
tracks=struct( 'day',{},'type' ,{},'lat',{},'lon',{}, ... 
    'j',{},'i',{},'shapes' ,{});        
search=struct( 'day',{},'type' ,{},'lat',{},'lon',{}, ... 
    'j',{},'i',{},'shapes' ,{});       
% initialize the warning array  
warn_tracks=[];  
% loop through all the days in which eddies were detected  
for i=1:length(centers)  
    day=centers(i).day;  
    disp([ 'Searching day ' ,num2str(day), ' ... ']) 
    % variables containing data of all the eddies detected for the  
    % current day  
    eddy_lat=centers(i).lat;  
    eddy_lon=centers(i).lon;  
    eddy_type=centers(i).type;  
    eddy_i=centers(i).i;  
    eddy_j=centers(i).j;  
    eddy_shapes=shapes(i).lonlat;  
         
    % Begin updating eddy tracks ----------------------------------  
    if ~isempty(eddy_type) % be sure eddies were present that specific day  
        % if first time -step, then all the eddies are open tracks  
        % (added to search)  
        if i==1 
            for i2=1:length(eddy_type)  
                search(i2).day=day;  
                search(i2).type=eddy_type(i2);  
                search(i2).lat=eddy_lat(i2);  
                search(i2).lon=eddy_lon(i2);  
                search(i2).j=eddy_j(i2);  
                search(i2).i=eddy_i(i2);  
                search(i2).shapes= eddy_shapes(i2);  
            end 
        % if not first day, then open tracks from previus days are updated  
        % with eddies detected for the current day  
        else 
            for i2=1:length(eddy_type)  
                % lon and lat of eddy center s and type at current time  
                ind_now(i2,1)=eddy_lon(i2);  
                ind_now(i2,2)=eddy_lat(i2);  
                ind_now(i2,3)=eddy_type(i2);  
            end 
            % first: loop all open tracks from the previous day  
            for i2=1:length(search)  
                % find lon and lat of the eddy centers in open tracks from  
                % previous day  
                if (day - search(i2).day(end))==1  
                    ind_old(1)=search(i2).lon(end);  
                    ind_old (2)=search(i2).lat(end);  
                    ind_old(3)=search(i2).type(end);  
                    % find centers at t2 that are within 'r' km from the  
                    % center at t1:  
                    % - create distance array  
                    dist=ones(size(ind_now(:,1)))*9999;  
                    % - compute distances  
                    for i3=1:size(ind_now,1)  
                       if(degree==0)  
                        dist(i3,1)=m_lldist([ind_now(i3,1) ind_old(1)], ... 
                            [ind_now(i3,2) ind_old(2)]);  
                       else 
                        dist(i3,1)=sqrt((ind_now(i3,1) -ind_old(1))^2+ ... 
                            (ind_now(i3,2) -ind_old(2))^2);  
                       end 
                    end 
                    % - find eddy of same type at less than 'r' km  
                    mov=find(dist<r & ind_now(:,3)==ind_old(3));  
                    % update open tracks with center from current day  
                    if ~isempt y(mov)  
                        % case in which more than one center is found  
                        % within the area  
                        if max(size(mov))>1  
                            % display and update warning  
                            disp([ 'Warning t+1 day: ' ,num2str(max(size(mov))), ' possibilities' , ... 
                                ', day ' ,num2str(day), ' !!!!!!' ]) 
                            warn_tracks(end+1,:)=[day 1 max(size(mov))];  
                            % track is updated with the closest center  
                            mov=find(dist==min(dist(mov)));  
                            % if still more than one eddy at same minimum  
                            % distance (very exceptional) then the first  
                            % eddy is chosen  
                            if max(size(mov))>1  
                                mov=mov(1);  
                            end 
                        end 
                        % add info of the new eddy center to the open  
                        % track array  
                        search(i2).day=[search(i2).day; day];  
                        search(i2).type=[search(i2).type; eddy_type(mov)];  
                        search(i2).lat=[search(i2).lat; eddy_lat(mov)];  
                        search(i2).lon=[search(i2).lon; eddy_lon(mov)];  
                        search(i2).j=[search(i2).j; eddy_j(mov)];  
                        search(i2).i=[search(i2).i; eddy_i(mov)];  
                        search(i2).shapes=[search(i2).shapes; eddy_shapes( mov)];  
                        % remove index of updated eddy center from array of  
                        % centers not updated yet  
                        ind_now(mov,:)=NaN;  
                    end 
                end 
            end 
            % secon d: loop all open tracks from two days before  
            % (enlarged searching radius 'r + r/2')  
            for i2=1:length(search)  
                if (day - search(i2).day(end))>1  
                    ind_old(1)=search(i2).lon(end);  
                    ind_old(2)=search(i2).lat(end);  
                    ind_old(3)=search(i2).type(end);  
                    % larger area  
                    r2=r+ r/2;  
                    % find centers at t2 that are within 'r2' km from the  
                    % center at t 1: 
                    % - create distance array  
                    dist=ones(size(ind_now(:,1)))*9999;  
                    % - compute distances  
                    for i3=1:size(ind_now,1)  
                      if(degree==0)  
                        dist(i3,1)=m_lldist([ind_now(i3,1) ind_old(1)], ... 
                            [ind_now(i3,2) ind_old(2)]);  
                       else 
                        dist(i3,1)=sqrt((ind_now(i3,1) -ind_old(1))^2+ ... 
                            (ind_now(i3,2) -ind_old (2))^2);  
                       end 
                    end 
                    % - find eddy of same type at less than 'r' km  
                    mov=find(dist<r2 & ind_now(:,3)==ind_old(3));  
                    % update open tracks with center from curre nt day  
                    if ~isempty(mov)  
                        % case in which more than one center is found  
                        % within the area  
                        if max(size(mov))>1  
                            % display and update warning  
                            disp([ 'Warning t+2 day: ' ,num2str(max(size(mov))), ' possibilities' , ... 
                                ', day ' ,num2str(day), ' !!!!!!' ]) 
                            warn_tracks(end+1,:)=[day 2 max(size(mov))];  
                            % track is updated with the closest center  
                            mov=find(dist==min(dist(mov)));  
                            % if still more than one eddy at same minimum  
                            % distance (very except ional) then the first  
                            % eddy is chosen  
                            if max(size(mov))>1  
                                mov=mov(1);  
                            end 
                        end 
                        % add info of  the new eddy center to the open  
                        % track array  
                        search(i2).day=[search(i2).day; day];  
                        search(i2).type=[search(i2).type; eddy_type(mov)];  
                        search(i2).lat =[search(i2).lat; eddy_lat(mov)];  
                        search(i2).lon=[search(i2).lon; eddy_lon(mov)];  
                        search(i2).j=[search(i2).j; eddy_j(mov)];  
                        search(i2).i=[search(i2).i; eddy_i(mov)];  
                        search(i2).shapes=[search(i2).shapes; eddy_shapes(mov)];  
                        % remove index of updated eddy center from array of  
                        % centers not updated yet  
                        ind_now(mov,:)=NaN;  
                    end 
                end 
            end 
            % remaining eddies from present day are new tracks;  
            % added to open track array  
            eddy_type=eddy_type(~isnan(ind_now(:,1)));  
            eddy_lat=eddy_lat(~isnan(ind_now(:,1)));  
            eddy_lon=eddy_lon(~isnan(ind_now(:,1)));  
            eddy_j=eddy_j(~isnan(ind_now(:,1)));  
            eddy_i=eddy_i(~isnan(ind_now(:,1)));  
            eddy_shapes=eddy_shapes(~isnan(ind_now(:,1)));  
            if ~isempty(eddy_type)  
                for i3=1:length(eddy_type)  
                    search(length(search)+1).day=day;  
                    search(length(search)).type=eddy_type(i3);  
                    search(length(search)).lat=eddy_lat(i3);  
                    search(length(search)).lon=eddy_ lon(i3);  
                    search(length(search)).j=eddy_j(i3);  
                    search(length(search)).i=eddy_i(i3);  
                    search(length(search)).shapes=eddy_shapes(i3);  
                end 
            end 
            % eddies are considered dissipated when tracks are not updated  
            % for longer than two days; tracks are then removed from open  
            % tracks array and stored in the closed track array.  
            for i2=1:length(search)  
                % find dissipated tracks  
                moved(i2)= search(i2).day(end) < day -1; 
                if moved(i2)  
                    % move tracks to closed tracks array  
                    tracks(length(tracks)+1)=search(i2);  
                end 
            end 
            % remove tracks from open track array  
            search(moved==1)=[];  
            % clear some variables (re -initialized at next loop)  
            clear ind_now  moved  
        end 
    end 
    % end updating eddy tracks for current day  
end 
% Add tracks that are still in the search array at the end of the  
% time -series  
for i=1:length(search)  
    tracks(length(tracks)+1)=search(i);  
end 
% remove tracks shorter than cut_off days  
for i=1:length(tracks)  
    short(i)=length(tracks(i).day)<cut_off;  
end 
tracks(short)=[];  
short=sum(short);  
  
for i=1:length(tracks)  
    day=tracks(i).day;  
    disp([num2str(i/length(tracks)*100), '%']); 
    for j=1:length(day)  
        lon=tracks(i).lon(j);  
        lat=tracks(i).lat(j);  
        shapes=tracks(i).shapes{j};  
    for k=1:length(shapes)  
        long(k)=m_lldist([lon,shapes(1,k)],[lat,shapes(2,k)]);  
    end 
        tracks(i).size(j)=nanmean(long);  
        clear long 
    end 
end 
  
% save variables  
save([path_out, 'removed_tracks' ],'short' ) 
save([path_out, 'eddy_tracks' ],'tracks' ,'-v6') 
save([path_out, 'warnings_tracks' ],'warn_tracks' ) 
% close log file  
diary off 
 
 
 
