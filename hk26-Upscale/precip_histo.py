import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.gridspec as gridspec
from matplotlib.pyplot import figure
import xarray as xr
import os
from matplotlib.colors import LogNorm
from matplotlib.colors import ListedColormap
import matplotlib.colors as mcolors
import datetime as dt

# gpm_DS = xr.open_dataset('/gws/nopw/j04/kscale/USERS/emg/data/DYAMOND_Summer/precip/GPM/gpm_DS_15NS_dm.nc')
# precip_15NS=gpm_DS['precipitation_rate']

# G9n1280_DS_precip = xr.open_dataset('/gws/nopw/j04/kscale/USERS/emg/data/DYAMOND_Summer/precip/global_GAL9_n1280/GAL9_n1280_DMn1280GAL9_precip_15NS_dm.nc')
# G9n1280_DS_15NS = G9n1280_DS_precip['precipitation_rate']

# R34p4_DS_precip = xr.open_dataset('/gws/nopw/j04/kscale/USERS/emg/data/DYAMOND_Summer/precip/channel_RAL3_4p4km/RAL3_4p4km_DMn1280GAL9_precip_15NS_dm.nc')
# R34p4_DS_15NS = R34p4_DS_precip['precipitation_rate']

# G9n2560G9DM_DS_precip = xr.open_dataset('/gws/nopw/j04/kscale/USERS/emg/data/DYAMOND_Summer/precip/channel_GAL9_n2560/GAL9_n2560_DMn1280GAL9_precip_15NS_dm.nc')
# G9n2560_DS = G9n2560G9DM_DS_precip['precipitation_rate']

# R3n2560_DS_precip = xr.open_dataset('/gws/nopw/j04/kscale/USERS/emg/data/DYAMOND_Summer/precip/channel_RAL3_n2560/RAL3_n2560_DMn1280GAL9_precip_15NS_dm.nc')
# R3n2560_DS = R3n2560_DS_precip['precipitation_rate']

gpm_DW = xr.open_dataset('/gws/nopw/j04/kscale/USERS/emg/data/DYAMOND_Winter/precip/GPM/gpm_DW_15NS_dm.nc')
precip_gpm_DW=gpm_DW['precipitation_rate']

G9n1280_DW_precip = xr.open_dataset('/gws/nopw/j04/kscale/USERS/emg/data/DYAMOND_Winter/precip/global_GAL9_n1280/GAL9_n1280_DMn1280GAL9_precip_15NS_dm.nc')
G9n1280_DW = G9n1280_DW_precip['precipitation_rate']

G9n2560_DW_precip = xr.open_dataset('/gws/nopw/j04/kscale/USERS/emg/data/DYAMOND_Winter/precip/channel_GAL9_n2560/GAL9_n2560_DMn1280GAL9_precip_15NS_dm.nc')
G9n2560_DW = G9n2560_DW_precip['precipitation_rate']

R3n2560_DW_precip = xr.open_dataset('/gws/nopw/j04/kscale/USERS/emg/data/DYAMOND_Winter/precip/channel_RAL3_n2560/RAL3_n2560_DMn1280GAL9_precip_15NS_dm.nc')
R3n2560_DW = R3n2560_DW_precip['precipitation_rate']

R34p4_DW_precip = xr.open_dataset('/gws/nopw/j04/kscale/USERS/emg/data/DYAMOND_Winter/precip/channel_RAL3_4p4km/RAL3_4p4km_DMn1280GAL9_precip_15NS_dm.nc')
R34p4_DW = R34p4_DW_precip['precipitation_rate']

############ DS PROCESSING #################
#Define logarithmically spaced bins
# min_precip = max(precip_15NS.min(), 0.01)  # Avoid zero in log scale
# max_precip = precip_15NS.max()
# num_bins = 100
# log_bins = np.logspace(np.log10(min_precip), np.log10(max_precip), num_bins)
# # Compute histogram
# bin_counts, bin_edges = np.histogram(precip_15NS, bins=log_bins, density=False)
# bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
# # Compute contribution to mean
# weighted_contributions = bin_centers * bin_counts  # Total precipitation contribution
# total_precipitation = np.sum(weighted_contributions)  # Sum of all precipitation values
# contribution_to_mean = weighted_contributions / total_precipitation  # Normalize
# print('gpm done')

# ############## GLOB10GAL ###################
# min_precipG10 = max(G9n1280_DS_15NS.min(), 0.01)  
# max_precipG10 = G9n1280_DS_15NS.max()
# log_binsG10 = np.logspace(np.log10(min_precipG10), np.log10(max_precipG10), num_bins)

# bin_countsG10, bin_edgesG10 = np.histogram(G9n1280_DS_15NS, bins=log_binsG10, density=False)
# bin_centersG10 = (bin_edgesG10[:-1] + bin_edgesG10[1:]) / 2

# weighted_contributionsG10 = bin_centersG10 * bin_countsG10
# total_precipitationG10 = np.sum(weighted_contributionsG10)
# contribution_to_meanG10 = weighted_contributionsG10 / total_precipitationG10
# print('GLOB10GAL done')

# ################ CTC5GAL #################
# min_precipG5 = max(G9n2560_DS.min(), 0.01) 
# max_precipG5 = G9n2560_DS.max()
# log_binsG5 = np.logspace(np.log10(min_precipG5), np.log10(max_precipG5), num_bins)

# bin_countsG5, bin_edgesG5 = np.histogram(G9n2560_DS, bins=log_binsG5, density=False)
# bin_centersG5 = (bin_edgesG5[:-1] + bin_edgesG5[1:]) / 2

# weighted_contributionsG5 = bin_centersG5 * bin_countsG5  
# total_precipitationG5 = np.sum(weighted_contributionsG5)  
# contribution_to_meanG5 = weighted_contributionsG5 / total_precipitationG5
# print('CTC5GAL done')

# ############### CTC5RAL ##############
# min_precipR5 = max(R3n2560_DS.min(), 0.01) 
# max_precipR5 = R3n2560_DS.max()
# log_binsR5 = np.logspace(np.log10(min_precipR5), np.log10(max_precipR5), num_bins)

# bin_countsR5, bin_edgesR5 = np.histogram(R3n2560_DS, bins=log_binsR5, density=False)
# bin_centersR5 = (bin_edgesR5[:-1] + bin_edgesR5[1:]) / 2

# weighted_contributionsR5 = bin_centersR5 * bin_countsR5  
# total_precipitationR5 = np.sum(weighted_contributionsR5)  
# contribution_to_meanR5 = weighted_contributionsR5 / total_precipitationR5  
# print('CTC5RAL done')

# ############ CTC4RAL ############
# min_precipC4 = max(R34p4_DS_15NS.min(), 0.01) 
# max_precipC4 = R34p4_DS_15NS.max()
# log_binsC4 = np.logspace(np.log10(min_precipC4), np.log10(max_precipC4), num_bins)

# bin_countsC4, bin_edgesC4 = np.histogram(R34p4_DS_15NS, bins=log_binsC4, density=False)
# bin_centersC4 = (bin_edgesC4[:-1] + bin_edgesC4[1:]) / 2

# weighted_contributionsC4 = bin_centersC4 * bin_countsC4  
# total_precipitationC4 = np.sum(weighted_contributionsC4)  
# contribution_to_meanC4 = weighted_contributionsC4 / total_precipitationC4  
# print('CTC4RAL done')

########### DW PROCESSING #################
# Define logarithmically spaced bins
min_precip = max(precip_gpm_DW.min(), 0.01)  # Avoid zero in log scale
max_precip = precip_gpm_DW.max()
num_bins = 100
log_bins = np.logspace(np.log10(min_precip), np.log10(max_precip), num_bins)
# Compute histogram
bin_counts, bin_edges = np.histogram(precip_gpm_DW, bins=log_bins, density=False)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
# Compute contribution to mean
weighted_contributions = bin_centers * bin_counts  # Total precipitation contribution
total_precipitation = np.sum(weighted_contributions)  # Sum of all precipitation values
contribution_to_mean = weighted_contributions / total_precipitation  # Normalize
print('gpm done')

############## GLOB10GAL ###################
min_precipG10W = max(G9n1280_DW.min(), 0.01)  
max_precipG10W = G9n1280_DW.max()
log_binsG10W = np.logspace(np.log10(min_precipG10W), np.log10(max_precipG10W), num_bins)

bin_countsG10W, bin_edgesG10W = np.histogram(G9n1280_DW, bins=log_binsG10W, density=False)
bin_centersG10W = (bin_edgesG10W[:-1] + bin_edgesG10W[1:]) / 2

weighted_contributionsG10W = bin_centersG10W * bin_countsG10W
total_precipitationG10W = np.sum(weighted_contributionsG10W)
contribution_to_meanG10W = weighted_contributionsG10W / total_precipitationG10W
print('GLOB10GAL done')

################ CTC5GAL #################
min_precipG5W = max(G9n2560_DW.min(), 0.01) 
max_precipG5W = G9n2560_DW.max()
log_binsG5W = np.logspace(np.log10(min_precipG5W), np.log10(max_precipG5W), num_bins)

bin_countsG5W, bin_edgesG5W = np.histogram(G9n2560_DW, bins=log_binsG5W, density=False)
bin_centersG5W = (bin_edgesG5W[:-1] + bin_edgesG5W[1:]) / 2

weighted_contributionsG5W = bin_centersG5W * bin_countsG5W  
total_precipitationG5W = np.sum(weighted_contributionsG5W)  
contribution_to_meanG5W = weighted_contributionsG5W / total_precipitationG5W
print('CTC5GAL done')

############### CTC5RAL ##############
min_precipR5W = max(R3n2560_DW.min(), 0.01) 
max_precipR5W = R3n2560_DW.max()
log_binsR5W = np.logspace(np.log10(min_precipR5W), np.log10(max_precipR5W), num_bins)

bin_countsR5W, bin_edgesR5W = np.histogram(R3n2560_DW, bins=log_binsR5W, density=False)
bin_centersR5W = (bin_edgesR5W[:-1] + bin_edgesR5W[1:]) / 2

weighted_contributionsR5W = bin_centersR5W * bin_countsR5W  
total_precipitationR5W = np.sum(weighted_contributionsR5W)  
contribution_to_meanR5W = weighted_contributionsR5W / total_precipitationR5W  
print('CTC5RAL done')

############ CTC4RAL ############
min_precipC4W = max(R34p4_DW.min(), 0.01) 
max_precipC4W = R34p4_DW.max()
log_binsC4W = np.logspace(np.log10(min_precipC4W), np.log10(max_precipC4W), num_bins)

bin_countsC4W, bin_edgesC4W = np.histogram(R34p4_DW, bins=log_binsC4W, density=False)
bin_centersC4W = (bin_edgesC4W[:-1] + bin_edgesC4W[1:]) / 2

weighted_contributionsC4W = bin_centersC4W * bin_countsC4W  
total_precipitationC4W = np.sum(weighted_contributionsC4W)  
contribution_to_meanC4W = weighted_contributionsC4W / total_precipitationC4W
print('CTC4RAL done')


##### PLOTTING #######
print('plotting...')

plt.scatter(bin_centers, contribution_to_mean, color='k', label="GPM", facecolors='none',s=15)
plt.scatter(bin_centersG10W, contribution_to_meanG10W, color='b', label="GLOB10GAL", facecolors='none',s=15)
plt.scatter(bin_centersG5W, contribution_to_meanG5W, color='c', label="CTC5GAL", facecolors='none',s=15)
plt.scatter(bin_centersR5W, contribution_to_meanR5W, color='darkorange', label="CTC5RAL", facecolors='none',s=15)
plt.scatter(bin_centersC4W, contribution_to_meanC4W, color='m', label="CTC4RAL", facecolors='none',s=15)

plt.xscale('log')
# plt.yscale('log')
plt.xlim([1e-2,2e2])
plt.ylim([1e-7,5e-2])
plt.xlabel("Precipitation rate (mm/h)")
plt.ylabel("Contribution to Mean (mm/h)")
plt.title("Precipitation histograms - DYAMOND Winter | diurnal cycle removed")
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
print('figure saved')
plt.savefig('/home/users/emg97/emgPlots/precip_histo_DW_allmodelsGPM_logx_dm.png',dpi=300,bbox_inches='tight')