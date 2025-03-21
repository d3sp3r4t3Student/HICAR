import xarray as xr
import os
import sys
sys.path.append(os.getcwd())
import HICAR_Domain as hd

###########################################################################
#######################  User options, to be edited  ######################
###########################################################################

#The horizontal resolution of the domain in meters
res = 10

# The target domain, including lat and lon variables named as "lat" and "lon", and
# a DEM labeled as "topo". Optionally, landuse and landmask variables should be specified here.
target_domain_fn = '/mnt/c/Users/sesselma/Code/Code_PhD/HICAR/domains/10m_domains/WGS84_domain_swissaltiregio10m_HEF_wlu_small.nc'
# A domain with extent ~20km beyond the borders of the above target domain. 
# Only lat,lon, and topo are required variables.
large_domain_fn = '/mnt/c/Users/sesselma/Code/Code_PhD/HICAR/domains/10m_domains/WGS84_domain_swissaltiregio10m_HEF_wlu.nc'
# Name of output file
output_domain_fn = '/mnt/c/Users/sesselma/Code/Code_PhD/HICAR/domains/10m_domains/Hintereis_10m_final.nc'

# These are used in the calculation of ridgelines, and can be tuned if the user
# is not satisfied with the deliniation of ridgelines in the output file
# terr_filter = 10
# TPI_thresh = 100

###########################################################################
############################  End user options  ###########################
###########################################################################

dom = xr.open_dataset(target_domain_fn)
dom_rad = xr.open_dataset(large_domain_fn)

dom_out = hd.wholeShebang(dom,dom_rad,res=res)

dom_out.to_netcdf(output_domain_fn)