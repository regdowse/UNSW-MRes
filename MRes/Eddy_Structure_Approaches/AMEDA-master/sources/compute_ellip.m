function [xbary,ybary,z,a,b,theta,lim]=compute_ellip(xy,grid_ll)
%[xbary,ybary,z,a,b,theta,lim]=compute_ellip(xy {,grid_ll})
%
% Compute the barycentre of a closed polygon 2xN array xy where
%   x or lon is in xy(1,:)
%   y or lat is in xy(2,:)
%   grid_ll =1  (default) if the coordinates are in (longitude,latitude)
%           =0 if coordinates are in 'km'
%
% and fit an ellipse using fitellipse
%
% (xbary,ybary) barycenter from coordinates of the polygon
% (z,a,b,theta) from fitellipse.m where alpha becomes theta here
% (lim) is the vertices number of the polygon output return nan if lim<=5
%
%-------------------------
%  Jan 2016 by B. LE VU
%-------------------------
%
%=========================

% Default grid is (lon,lat)
if nargin==1
    grid_ll = 1;
end

%----------------------------------------
% Size of the polygon
lim = size(xy,2)-1;

%----------------------------------------
% Barycenter computation
xbary = mean(xy(1,1:lim));
ybary = mean(xy(2,1:lim));

%----------------------------------------
% Ellipse fitting
if lim >= 5
    try
        %----------------------------------------
        % Coord = coordinates of the polygon
        if grid_ll

            xs(1) = xbary;
            ys(1) = ybary;

            % initialise
            coord = zeros(2,lim);

            for pt=1:lim

                xs(2) = xy(1,pt);
                ys(2) = xy(2,pt);

                %----------------------------------------
                % distances in km of every point from the barycenter
                coord(1,pt) = sign(diff(xs)) * sw_dist2([ybary ybary],xs);
                coord(2,pt) = sign(diff(ys)) * sw_dist2(ys,[xbary xbary]);
            end

            [~, a, b, theta] = fitellipse(coord,'linear');
            z = [xbary,ybary];
            
        else

            coord = [xy(1,1:lim);xy(2,1:lim)];
            [z, a, b, theta] = fitellipse(coord,'linear');
            
        end

        if a<0 || b<0 || z(1)==0
            z = [NaN NaN];
            a = NaN;
            b = NaN;
            theta = NaN;
        end
        
    catch
        z = [NaN NaN];
        a = NaN;
        b = NaN;
        theta = NaN;
    end
    
else
    
    z = [NaN NaN];
    a = NaN;
    b = NaN;
    theta = NaN;
    
end

