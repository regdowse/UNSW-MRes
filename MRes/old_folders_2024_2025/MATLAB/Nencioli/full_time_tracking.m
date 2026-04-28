function  full_time_tracking(postfix,yearmin,yearmax,Main_dir,path_out_all,degree,r)  
% full_time_tracking .m (call function ) 
% According to the identified eddy  center and boundary data, the data of more than 1 year is merged and 
tracked across years  
% postfix  is detection eddy area code  
% yearmin  and yearmax  are y ear of beginning  and ending  of eddy detection  
% Main_dir  is output directory  
% path_out_all  is total output directory   
% degree  is in one latitude  
% r is distance in degree to define the area u sed to derive eddy tracks  
if yearmin~=yearmax  
    for i=yearmin:yearmax  
        file_in=[Main_dir, 'Eddy_' ,postfix, '/Tracks_' ,num2str(i), '/']; 
% Read eddy  center and boundary data  
        load([file_in, 'eddy_centers.mat' ]) 
        load([file_in, 'eddy_shapes.mat' ]) 
        if i==yearmin  
            c_all=centers;  
            s_all=shapes;  
            clear centers  shapes  
        else 
            n=length(c_all)+1;  
            for ii=1:length(centers)  
                c_all(n)=centers(ii);  
                s_all(n)=shapes(ii);  
                n=n+1;  
            end 
            clear centers  shapes  
        end 
    end 
    centers=c_all;  
    shapes=s_all;  
    clear c_all  s_all  
    path_out=path_out_all;  
    clear path_out_all  
    save( [path_out, 'eddy_centers.mat' ],'centers' ,'-v7.3' ); 
    save([path_out, 'eddy_shapes.mat' ],'shapes' ,'-v7.3' );     
    mod_eddy_tracks (degree,r,path_out)  
end 
 
