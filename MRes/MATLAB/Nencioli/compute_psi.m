function  psi=compute_psi(u,v,km_dlon,km_dlat)  
% compute_psi.m  (call function ) 
% compute_psi  (u,v,km_dlon,km_dlat) computes the streamfunction (PSI) field  by spatially integrating the 
u- and v -component of velocity within the area around the detected eddy center defined in eddy_dim.m. . 
The underlaying assumption is that in the presence of an eddy the  the velocity field is characterized by 
weak divergence, so  that contours  of PSI are tangential to the velocity vectors.  To reduce the error 
associated with this assumption PSI is computed  integrating u and v along two different paths, and the two 
fields are  then averaged.  
% u and v are NxM matrices of the two com ponent of velocity;  
% km_dlon is Nx1 vector containing the longitudinal spacing in km between grid points. (It can vary with 
latitude)  
% km_dlat is the latitudinal spacing in km between grid points. (It is assumed constant)  
% OUTPUT:  
%-psi is the NxM matri x of the streamfunction field  
% See also Appendix in Nencioli et al. (2009) for further details.  
 
%-------------------------  
%   Ver. 1.2 May.2010  
%   Ver. 1.1 Dec.2009  
%   Authors: Francesco Nencioli, francesco.nencioli@opl.ucsb.edu,  
%            Charles  Dong, cdong@atmos.ucla.edu  
%-------------------------  
  
% Indices for domain size  
lx=size(u,1);  
ly=size(u,2);  
% itegrate first row of v along longitude (first term of eq.A2)  
cx=cumtrapz(v(1,:))*km_dlon(1); % trapezoidal sum  
% integrate first column of u a long latitude (second term of eq.A3)  
cy=cumtrapz(u(:,1))*km_dlat;  
% expand the vectors into matrices to compute PSI  
cx=cx(ones(lx,1),:);  
cy=cy(:,ones(1,ly));  
% create km_lon matrix for spacial integration of v  
km_dlon=repmat(km_dlon,1,size(v,2));  
% compute streamfunction -----------------------  
% PSI from integrating v firts and then u  
psi_xy=( -cx + cumtrapz(u)*km_dlat); %(eq. A2)  
% PSI from integrating u first and then v  
psi_yx=( -cumtrapz(v,2).*km_dlon + cy); %(eq. A3)  
% final PSI as average betw een the two  
psi=(psi_xy+psi_yx)/2;  
%---------------------------------------------  
  
%psi=((cumtrapz(v,2).*km_dlon -cy)-(cumtrapz(u)*km_dlat -cx))/2;  
 
