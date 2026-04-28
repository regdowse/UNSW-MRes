% eddy  and sea -air-biological multi -parameter matching (here is an example of sst  data)  
clear;  
clc; 
tic; 
warning off; 
% eddy tracks file  
eddy_tracks_file= 'eddy_tracks.mat' ; 
% sst data file  
file_dir= '/Data_sst/' ; 
% output directory  
data_out_dir= '/match /'; 
mkdir(data_out_dir);  
% Open the  eddy  tracking file  
load(eddy_tracks_file);  
% Read sst file time list  
sst_files=dir(file_dir);  
sst_files(1:2,:)=[];  
for n=1:length(sst_files)  
    stime(n,:)=str2num(sst_files(n).name(5:12));  
end 
clear n 
% Cycle matching from the first eddy  
for i=1:length(tracks)  
disp([num2str(i/length(tracks)*100), '%']); 
% Screen for eddies with a lifetime of 4 weeks or more  
if length(tracks(i).day)>=28  
% Create a eddy  matching multiparameter information structure  
      match_eddy=repmat(struct( 'index' ,{},'day',{},'type' ,{},'radius' ,{}, 'shapes' ,{},'c_lat' ,{},'c_lon' ,{}, 
'sst',{}, 'sst_lon' ,{}, 'sst_lat' ,{},'ssu',{}, 'ssv',{},'uv_lon' ,{}, 'uv_lat' ,{}),1,10);  
        clear eddy_index  
        eddy_index=[ 'MS_' ,num2str(i, '%06d' )]; 
        num=1;  
   % Read eddy  daily information  
        for i_day=1:length(tracks(i).day)  
            day=tracks(i).day(i_day);  
            type=tracks(i).type(i_day);  
            c_lat=tracks(i).lat(i_day);  
            c_lon=tracks(i).lon(i_day ); 
            radius=tracks(i).size(i_day);  
            shape=tracks(i).shapes{i_day,1};  
            shape_lon=shape(1,:);  
            shape_lat=shape(2,:);  
% Complete the conversion of east -west cross -border 0 degree line  
            if c_lon<180  
                shape_lon(shape_lon>180)=shape_lon(shape_lon>180) -360; 
            else 
                shape_lon(shape_lon<180)=shape_lon(shape_lon<180)+360;  
            end 
% Set the eddy  to match the background field range (this is 2.5 times the radius of the eddy ) 
            R=max([max(abs(shape_lat -c_lat)),max(abs(shape_lon -c_lon))]);  
            lat_max=c_lat+R*2.5;  
            lat_min=c_lat -R*2.5;  
            lon_max=c_lon+R*2.5;  
            lon_min=c_lon -R*2.5;  
% Find background field data that matches th e eddy  
            timestr = datestr(day, 'yyyymmdd' ); 
            eddy_time = str2num(timestr);  
            index=find(eddy_time==stime);  
            if ~isempty(index)  
                sst_file=[file_dir,sst_files(index).name];  
                load( sst_file)  
                clear index  eddy_time  
    %Select background field data within the eddy  matching range  
                if lon_max>360  
                    lon(lon<180)=lon(lon<180)+360;  
                    lon=[lon(721:end,:);lon(1:720,:)];  
                    sst=[sst(721:end,:);sst(1:720,:)];  
                    index_lon=find(lon>=lon_min & lon<=lon_max);  
                    index_lat=find(lat>=lat_min & lat<=lat_max);  
                    sst_lon=lon(index_lon);  
                    sst_lon( sst_lon>360)=sst_lon(sst_lon>360) -360; 
                    sst_lat=lat(index_lat);  
                    sst=sst(index_lon,index_lat);  
                    clear lon lat index_lon  index_lat  
                elseif  lon_min<0  
                    lon(lon>180)=lon (lon>180) -360; 
                    lon=[lon(721:end,:);lon(1:720,:)];  
                    sst=[sst(721:end,:);sst(1:720,:)];  
                    index_lon=find(lon>=lon_min & lon<=lon_max);  
                    index_lat=find(lat>=lat_min & lat<=lat_max);  
                    sst_lon=lon(index_lon);  
                    sst_lon(sst_lon<0)=sst_lon(sst_lon<0)+360;  
                    sst_lat=lat(index_lat);  
                    sst=sst(index_lon,index_lat);  
                    clear lon lat index_lon  index_lat  
                else 
                    index_lon=find(lon>=lon_min & lon<=lon_max);  
                    index_lat=find(lat>=lat_min & lat<=lat_max);  
                    sst_lon=lon(index_lon);  
                    sst_lat=lat(index_lat);  
                    sst=sst(index_lon,index_lat);  
                    clear lon lat index_lon  index_lat  
                end 
            else 
                sst_lon=nan;  
                sst_lat=nan;  
                sst=nan;  
            end 
% Store the physical property information of the eddy  itself with the matched background 
field 
            match_eddy(1).ID=eddy_index;  
            match_eddy(1).day(num)=day;  
            match_eddy(1).type(num)=type;  
            match_eddy(1).radius(num)=r adius;  
            match_eddy(1).shapes(num)=mat2cell(shape,size(shape,1),size(shape,2));  
            match_eddy(1).c_lat(num)=c_lat;  
            match_eddy(1).c_lon(num)=c_lon;  
            match_eddy(1).sst(num)=mat2cell(sst,size(sst,1),size(sst,2));  
            match_eddy(1).sst_lat(num)=mat2cell(sst_lat,size(sst_lat,1),size(sst_lat,2));  
            match_eddy(1).sst_lon(num)=mat2cell(sst_lon,size(sst_lon,1),size(sst_lon,2));  
            clear day type radius  shape  shape_lon  shape_lat  lon_res  daystr  year timestr  c_lon  c_lat  
            clear sst_file  R lat_max  lat_min  lon_max  lon_min  sst_lat  sst_lon  sst nc ssu ssv uv_lon  
uv_lat  
            num=num+1;  
        end 
% Save output after matching based on a whole life cycle of eddie s 
        if ~isempty( match_eddy)  
            save([data_out_dir,eddy_index, '_match_eddy.mat' ],'match_eddy' ,'-v7.3' ); 
            clear match_eddy  eddy_index  
        end 
        toc; 
    end 
end 
 
