function  [eddy_uv eddy_c eddy]=uv_search(nc_u,nc_v,lon,lat,mask,day,a,b)  
% uv_search.m  (mtool_quiverkey)  
% uv_search(nc_u,nc_v,lon,lat,mask,day,a,b) identifies the points in the  domain which satisfy the four 
velocity constraints.  
% nc_u and nc_v  are the time series of the 2D u and v velocity field;  
% on and lat are longitudes and latitudes of points in the domain;  
% mask is the matrix that defines sea (1) and land points (0);  
% day is the time index that defines the 2D velocity field to analize;  
% a nad b are the two parameters used to apply the four constraints;  
% OUTPUT:  
% eddy_uv and eddy_c are the positions of the points that satisfy the  first two, and the first three 
constraints, respectively. They are  output only for debugging purpose (if ne eded).  
% eddy is the vector containing the poisitions of the estimated eddy  centers and the type of eddy 
(cyclonic=1; anticyclonic= -1), which are then saved in mod_eddy_centers.m (center positions are saved as 
latitude and longitude)  
%--------------------- ---- 
%   Ver. 1.2 May.2010  
%   Ver. 1.1 Dec.2009  
%   Authors: Francesco Nencioli, francesco.nencioli@opl.ucsb.edu,  
%            Charles Dong, cdong@atmos.ucla.edu  
%-------------------------  
  
% load 2D velocity fields ---------------  
nc=netcdf(nc_u);  
u=nc{ 'ssu'}(day,:,:);  
close(nc)  
  
nc=netcdf(nc_v);  
v=nc{ 'ssv'}(day,:,:);  
close(nc)  
%----------------------------------------  
% mask velocity fields  
v(mask==0)=NaN;  
u(mask==0)=NaN;  
% this parameter prevents the constraints to be applied to points too close  
% to the domain boundaries which would result in an index error  
borders=max(a,b)+1;  
% Added option for periodic East -West boundaries  
% (first and last column of the domain are identical)  
global  periodic  
if periodic==1  
    % u, v, lat and lon matrices are expan ded by adding the ending columns  
    % to the beginning of the domain, and the beginning columns at its end.  
    % this way the four constraints are applied also to the points closer  
    % than 'borders' points to the domain boundaries  
    v=[v(:,(end -borders+1:end) -1),v,v(:,(1:borders)+1)];  
    u=[u(:,(end -borders+1:end) -1),u,u(:,(1:borders)+1)];  
    lon(:,end)=lon(:,1);  
    lat(:,end)=lat(:,1);  
    lon=[lon(:,(end -borders+1:end) -1),lon,lon(:,(1:borders)+1)];  
    lat=[lat(:,(end -borders +1:end) -1),lat,lat(:,(1:borders)+1)];  
end 
% compute velocity magnitude (used in third constraint)  
vel=sqrt(u.^2+v.^2);  
%------------------------  
% This further masking should be uncommented only if there are problems  
% with mask  
% u(vel==0)=NaN;  
% v(vel==0 )=NaN;  
% vel(vel==0)=NaN;  
%------------------------  
% initialize arrays for eddy center position  
eddy_uv=[0 0];  
eddy_c=[0 0];  
eddy=[0 0 0];  
% this parameter prevents index errors when the third and fourth  
% constraints are applied close to the domain boundaries  
bound=size(vel);  
  
% loop through all latitudinal sections of v  
% (borders defines the first and last sections where the constraints can be  
% applied within the domain)  
for i=borders:1:length(v(:,1)) -borders+1  
    wrk=v(i,:); %latitudinal section of v  
    % first constraint: v component ---------------------------  
    % criteria 1: find indices of points of zero crossing  
    s=sign(wrk);  
    indx=find(diff(s)~=0 & ~isnan(diff(s)));  
    % remove points too close to the East and West domain b oundary  
    indx(indx<borders | indx>length(wrk) -borders)=[];  
    % criteria 2: inspect increase of v velocity "a" points away from the  
    % zero crossing points  
    % (var is used to identify cyclonic (1) and anticiclonicy ( -1) eddies;  
    % if a point f ails to satisfy one of the constraints var is set to 0)  
    for ii=1:length(indx)  
        % anticyclonic  
        if wrk(indx(ii))>=0  
            if wrk(indx(ii) -a)>wrk(indx(ii)) ... 
                    && wrk(indx(ii)+1+a)<wrk(indx(ii)+1)  
                var=-1; 
            else 
                var=0;  
            end 
        % cyclonic  
        elseif  wrk(indx(ii))<0  
            if wrk(indx(ii) -a)<wrk(indx(ii)) ... 
                    && wrk(indx(ii)+1+a)>wrk(indx(ii)+1)  
                var=1;  
            else 
                var=0;  
            end 
        end 
    %------------------------------------------------------------  
    % second constraint: u component ----------------------------  
        % inspect reversal of u and its increase "a" points away for  the 
        % points that satisfy the first constraint;  
        % if the criterion is respected the point is saved in eddy_uv;  
        % anticyclonic  
        if var== -1 
            if (u(i-a,indx(ii))<=0 && u(i -a,indx(ii))<=u(i -1,indx(ii)) && ... 
                    u(i+a,indx(ii))>=0 && u(i+a,indx(ii))>=u(i+1,indx(ii))) || ... 
                    (u(i-a,indx(ii)+1)<=0 && u(i -a,indx(ii)+1)<=u(i -1,indx(ii)+1) && ... 
                    u(i+a,indx(ii)+1)>=0 && u(i+a,indx(ii)+1)>=u(i+1,indx(ii)+1))  
                var=-1; 
                eddy_uv=[eddy_uv; ... 
                    lat(i,indx(ii)) lon(i,indx(ii)); ... 
                    lat(i,indx(ii)+1) lon(i,indx(ii)+1)];  
            else 
                var=0;  
            end 
        % cyclonic  
        elseif  var==1  
            if (u(i-a,indx(ii))>=0 && u(i -a,indx(ii))>=u(i -1,indx(ii)) && ... 
                    u(i+a,indx(ii))<=0 && u(i+a,indx(ii))<=u(i+1,indx(ii))) || ... 
                    (u(i-a,indx(ii)+1)>=0 && u(i -a,indx(ii)+1)>=u(i -1,indx(ii)+1) && ... 
                    u(i+a,indx(ii)+1)<=0 && u(i+a,indx(ii)+1)<=u(i+1,indx(ii)+1))  
                var=1;  
                eddy_uv=[eddy_uv; ... 
                    lat(i,indx(ii)) lon(i,indx(ii)); ... 
                    lat(i,indx(ii)+1) lon(i,indx(ii)+1 )]; 
            else 
                var=0;  
            end 
        end 
    %------------------------------------------------------------  
    % third constraint: velocity minimum ------------------------  
        % find the velocity minimum within the searc hing area defined by  
        % "b" around the points that satisfy the first two constraints  
        if var~=0  
            % velocity magnitude, latitude and longitude within the  
            % searching area  
            srch=vel(i -b:i+b,indx(ii) -b:indx (ii)+1+b);  
            slat=lat(i -b:i+b,indx(ii) -b:indx(ii)+1+b);  
            slon=lon(i -b:i+b,indx(ii) -b:indx(ii)+1+b);  
            % position of the velocity minimum within the searching area  
            [X Y]=find(srch==min(min(srch)));  
            % second searching area centered around the velocity minimum  
            % (bound prevents this area from extending outside the domain)  
            srch2=vel(max((i -b)+(X -1)-b,1):min((i -b)+(X -1)+b,bound(1)), ... 
                max((indx(ii) -b)+(Y -1)-b,1): min((indx(ii) -b)+(Y -1)+b,bound(2)));  
            % if the two minima coincide then it is a local minima  
            if ~(min(min(srch2))==min(min(srch)))  
                var=0;  
            else 
                eddy_c=[eddy_c; slat(X(1),Y(1)) slon(X(1),Y(1) )]; 
            end 
        end 
    %------------------------------------------------------------  
    % fourth constraint: vector rotation ------------------------  
        % check the rotation of the vectors along the boundary of the area  
        % "a -1" around the points which satisfy the first three constraints  
        d=a-1; 
        if var~=0  
            % indices of the estimated center in the large domain  
            [i1 i2]=find(lat==slat(X(1),Y(1)) & lon==slon(X(1),Y(1)));  
            % velocities w ithin "a -1" points from the estimated center  
            u_small=u(max(i1 -d,1):min(i1+d,bound(1)), ... 
                max(i2 -d,1):min(i2+d,bound(2)));  
            v_small=v(max(i1 -d,1):min(i1+d,bound(1)), ... 
                max(i2 -d,1):min(i2+d,bound(2)) ); 
            % constraint is applied only if sea -points are within the area  
            if isempty(find(isnan(u_small)))  
                % boundary velocities  
                u_bound=[u_small(1,:),u_small(2:end,end)', ... 
                    u_small( end,end -1:-1:1),u_small(end -1:-1:1,1)'];  
                v_bound=[v_small(1,:),v_small(2:end,end)', ... 
                    v_small(end,end -1:-1:1),v_small(end -1:-1:1,1)'];  
                % vector defining which quadrant each boundary vector  
                % belongs to  
                quadrants=zeros(size(u_bound));  
                quadrants(u_bound>=0 & v_bound>=0)=1;  
                quadrants(u_bound<0 & v_bound>=0)=2;  
                quadrants(u_bound<0 & v_bound<0)=3;  
                quadrants(u_bound >=0 & v_bound<0)=4;  
                % used identify which is the firts fourth quadrant vector  
                spin=find(quadrants==4);  
                % apply the constraint only if complete rotation and not  
                % all vectors in the fourth quadrant  
                if ~isempty(spin) && ... 
                        size(spin,2)~=size(quadrants,2)  
                    % if vectors start in 4 quadrant, then I add 4 to all  
                    % quandrant positions from the first 1 occurrence  
                    if spin(1)==1  
                        spin=find(quadrants~=4);  
                        spin=spin(1) -1; 
                    end 
                    quadrants(spin(end)+1:end)=quadrants(spin(end)+1:end)+4;  
                    % inspect vector rotation:  
                    % - no consecutive vectors more than one quadrant away  
                    % - no backward rotation  
                    if isempty(find(diff(quadrants)>1 , 1)) && ... 
                            isempty(find(diff(quadrants)<0 , 1)) 
                        eddy=[eddy;slat(X(1),Y(1)) slon(X(1),Y(1)) var];  
                    end 
                end 
            end 
        end 
    %------------------------------------------------------------  
    end 
end 
%-------------  work on "eddy" -------------------------------  
% eddy centers can be recorded more than once since the process is  
% repeated for every latitudinal section of v;  
  
% sort eddies for lon and lat  
eddy=sortrows(eddy,2);  
eddy=sortrows(eddy,1);  
% keep only one value per  eddy center  
good=find(diff(eddy(:,1))~=0 | diff(eddy(:,2))~=0);  
eddy=eddy([1;good+1],:);  
% remove zeros row from eddy  
zero_eddy=find(eddy(:,1)==0 & eddy(:,2)==0 & eddy(:,3)==0);  
eddy(zero_eddy,:)=[];  
% Adjust flags for cyclones/anticyclones for the Southe rn Emisphere  
% (it assumes latitudes ranging from -90 to 90)  
eddy(eddy(:,1)<0,3)= -eddy(eddy(:,1)<0,3);  
%----------------------------------------------------------  
  
%-------------  same for "eddy_uv" and "eddy_c" ------------  
% sort eddies for lon and lat  
eddy_uv=sortrows(eddy_uv,2);  
eddy_uv=sortrows(eddy_uv,1);  
% keep only one value per eddy center  
good=find(diff(eddy_uv(:,1))~=0 | diff(eddy_uv(:,2))~=0);  
eddy_uv=eddy_uv([1;good+1],:);  
% remove zeros row from eddy  
zero_eddy_uv=find(eddy_uv(:,1)==0 & eddy_u v(:,2)==0);  
eddy_uv(zero_eddy_uv,:)=[];  
%----------------------------------------------------------  
% sort eddies for lon and lat  
eddy_c=sortrows(eddy_c,2);  
eddy_c=sortrows(eddy_c,1);  
% keep only one value per eddy center  
good=find(diff(eddy_c(:,1))~=0 | diff(eddy_c(:,2))~=0);  
eddy_c=eddy_c([1;good+1],:);  
% remove zeros row from eddy  
zero_eddy_c=find(eddy_c(:,1)==0 & eddy_c(:,2)==0);  
eddy_c(zero_eddy_c,:)=[];  
%----------------------------------------------------------  
end 
 
