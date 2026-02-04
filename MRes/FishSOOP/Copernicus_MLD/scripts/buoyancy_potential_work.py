import numpy as np
from scipy import interpolate

"""
Function to calculate the WORK DONE BY THE BUOYANCY FORCE (WB) in
vertically displacing a water parcel from different depths to a reference
depth zref.
"""
def buoyancy_potential_work(rho, z, zref=-10):
	"""
	WB quantifies the water column's vertical homogeneity in terms of the work
	done by the buoyancy force in vertically displacing a water parcel under 
	static instability conditions.

	Inputs:
	* rho: Sigma-0. Potential Density Anomaly profile referred to 0 dbar (kg m^-3)
	* z: Depths in negative values (m)
	* zref: Reference depth (default -10 m). WB(zref)=0 by definition.
	The data must be a column vector and be arranged from the greatest to the shallowest depth.

	Outputs:
	* WB: Vertical profile of work done by buoyancy (J m^-3)
	* z: Depths (m)
	"""
	
	# If there is no data at zref, it is interpolated
	if not np.any(z == zref):
		# Interpolate rho at zref
		f = interpolate.interp1d(z, rho, kind='linear', fill_value='extrapolate')
		rho_zref = f(zref)

		# Construct the vectors z and rho, adding zref and rho_zref
		z_extended = np.insert(z, 0, zref)
		rho_extended = np.insert(rho, 0, rho_zref)

		# Sort vectors z and rho
		order = np.argsort(z_extended)
		z = z_extended[order]
		rho = rho_extended[order]
	
	# Defining variables
	WB = np.zeros_like(z)*np.nan # Work done by buoyancy (J m^-3)
	nz  = rho.shape[0] # Number of data in the vertical
	g = 9.81 # Acceleration of gravity (m s^-2)
	izref = np.nonzero(z==zref) # Index of zref in the vector z
	dz = np.diff(z) # Differences between adjacent elements of z (m), with positive units
	
	# Calculation of the vertical profile of WB
	for it in range(nz): # Loop for every depth
		zref = z[it]
		if np.isnan(zref): # If the depth of analysis is nan, WB is also nan
			WB[it] = np.nan
		else:
			rho_int = rho[it]
			
			# Computation of WB
			S = np.nan*np.zeros_like(z)
			S[it] = 0
			
			for i in range(0,it):
				sr=rho[i:it]+rho[i+1:it+1]
				sdz=dz[i:it]
				S[i] = -0.5*np.nansum(sr*sdz)
			for i in range(it+1,nz):
				sr=rho[it:i]+rho[it+1:i+1]
				sdz=dz[it:i]
				S[i] = 0.5*np.nansum(sr*sdz)

			V = -(g*(z-zref)*rho_int)+(g*S) # Buoyancy potential (J m^-3)
			V = -V # Change the sign of the buoyancy potential
			WB[it] = V[izref] # WB at the depth of analysis (J m^-3)
	
	return WB, z
