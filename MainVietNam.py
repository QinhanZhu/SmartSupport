# from climada_petals.hazard.river_flood import RiverFlood
from Functions import intro_haz_riverflood, tracks_ID
from climada.entity import Exposures, LitPop
from climada_petals.hazard import TCSurgeBathtub
import os
import numpy as np
from climada.hazard import Centroids
import time


project_folder = '/Users/qinhan/Documents/IIASA/IIASA-ISF-ETH/GSSP/UltimateOutputs/'  # results and data will be saved in this folder
input_folder = '/Users/qinhan/Documents/IIASA/IIASA-ISF-ETH/GSSP/Inputs/'  # this folder should be the folder containing the input files
# hazard_folder = input_folder + 'pakistan_flood/'
file_identifier = '_v01'
region_id = 704
country = 'VNM'

#%%
class Object(object):
    pass
hist_storms = Object()
hist_storms.lon_min = 102
hist_storms.lat_min = 8.5
hist_storms.lon_max = 110
hist_storms.lat_max = 23.5
hist_storms.start_year = 1983
hist_storms.end_year = 2023
horizon = 2050

EXP_norm_30 = LitPop.from_countries(countries=country, res_arcsec=30, reference_year=2020,fin_mode='norm')

EXP_norm_gdf = EXP_norm_30.gdf

quantile90 = EXP_norm_gdf['value'].quantile(0.90)

value_boolean = EXP_norm_gdf['value'] >= quantile90

# sum(EXP_norm_gdf['value'][value_boolean])

# sum(value_boolean)

lon_value = EXP_norm_gdf.centroid.x[value_boolean].to_numpy()

lat_value = EXP_norm_gdf.centroid.y[value_boolean].to_numpy()

VNM_centr= Centroids.from_lat_lon(EXP_norm_gdf.centroid.y,EXP_norm_gdf.centroid.x)

VNM_centr.set_elevation('/Users/qinhan/Documents/IIASA/IIASA-ISF-ETH/GSSP/Inputs/dem_compress.tif')

boolean_coast = VNM_centr.elevation <= 20

lat_coast = VNM_centr.lat[boolean_coast]

lon_coast = VNM_centr.lon[boolean_coast]

# centr_coast = Centroids.from_lat_lon(lat_coast, lon_coast)

# centr_coast.plot()

# centr_value = Centroids.from_lat_lon(lat_value, lon_value)

# centr_value.plot()


EXP_norm_150 = LitPop.from_countries(countries=country, res_arcsec=150, reference_year=2020,fin_mode='norm')

lat_150 = EXP_norm_150.gdf.centroid.y

lon_150 = EXP_norm_150.gdf.centroid.x

lat_final = np.concatenate((lat_coast,lat_value,lat_150),axis=None)

lon_final = np.concatenate((lon_coast,lon_value,lon_150),axis=None)

final_centr = Centroids.from_lat_lon(lat_final, lon_final)

# final_centr.plot()

ID_tracks_VNM = tracks_ID(hist_storms,'WP')
#%%

from climada.hazard.trop_cyclone import TropCyclone
from climada.hazard.tc_tracks import TCTracks

# fl_HAZ_tc = os.path.join(project_folder, 'HAZ_tc_hist_' + country + file_identifier + '.hdf5')

fl_HAZ_tc = '/Users/qinhan/Documents/IIASA/IIASA-ISF-ETH/CLIMADA_output/HAZ_tc_mix50_VNM_v01.hdf5'

if os.path.isfile(fl_HAZ_tc):

    Haz_tc_hist = TropCyclone.from_hdf5(fl_HAZ_tc)
    
else:
    
    tracks_VNM = TCTracks.from_ibtracs_netcdf(storm_id=ID_tracks_VNM)
    print('tracks selected')
    tracks_VNM.equal_timestep()
    print('equal timestep')
    tracks_VNM.calc_perturbed_trajectories(nb_synth_tracks=24)
    print('probabilistic events generated')
    start = time.time()
    Haz_tc_hist = TropCyclone.from_tracks(tracks_VNM, centroids=final_centr)
    print('HAZ finished')
    end = time.time()
    print(end-start)
    Haz_tc_hist.check()

    Haz_tc_hist.write_hdf5(fl_HAZ_tc)

#%%
Haz_tc_hist.plot_intensity(0)
#%%

fl_HAZ_tc_futr45 = os.path.join(project_folder, 'HAZ_tc_futr45_' + country + file_identifier + '.hdf5')

# fl_HAZ_tc = '/Users/qinhan/Documents/IIASA/IIASA-ISF-ETH/CLIMADA_output/HAZ_tc_mix50_final_MDG_v01.hdf5'

if os.path.isfile(fl_HAZ_tc_futr45):

    Haz_tc_futr45 = TropCyclone.from_hdf5(fl_HAZ_tc_futr45)
    
else:

    Haz_tc_futr45 = Haz_tc_hist.apply_climate_scenario_knu(rcp_scenario=45)

    Haz_tc_futr45.write_hdf5(fl_HAZ_tc_futr45)


#%%

fl_HAZ_tc_futr85 = os.path.join(project_folder, 'HAZ_tc_futr85_' + country + file_identifier + '.hdf5')

# fl_HAZ_tc = '/Users/qinhan/Documents/IIASA/IIASA-ISF-ETH/CLIMADA_output/HAZ_tc_mix50_final_MDG_v01.hdf5'

if os.path.isfile(fl_HAZ_tc_futr85):

    Haz_tc_futr85 = TropCyclone.from_hdf5(fl_HAZ_tc_futr85)
    
else:

    Haz_tc_futr85 = Haz_tc_hist.apply_climate_scenario_knu(rcp_scenario=85)

    Haz_tc_futr85.write_hdf5(fl_HAZ_tc_futr85)

    
#%%

fl_EXP_GA = os.path.join(project_folder, 'EXP_GA_' + country + file_identifier + '.hdf5')

if os.path.isfile(fl_EXP_GA):

    EXP_GA = LitPop.from_hdf5(fl_EXP_GA)

else:
    EXP_GA = LitPop.from_countries(countries=country, res_arcsec=30, reference_year=2020)
    # EXP_GA.gdf['impf_TC'] = np.ones(len(EXP_GA.gdf['value'])).astype(int)
    # EXP_GA.gdf['impf_RF'] = np.ones(len(EXP_GA.gdf['value'])).astype(int)
    # EXP_GA.set_geometry_points()
    EXP_GA.write_hdf5(fl_EXP_GA)

#%%
# sum(EXP_GA.gdf.value)

#%%
# fl_EXP_road = os.path.join(project_folder, 'EXP_road_' + country + file_identifier + '.hdf5')

fl_EXP_edu = os.path.join(project_folder, 'EXP_edu_' + country + file_identifier + '.hdf5')
fl_EXP_health = os.path.join(project_folder, 'EXP_health_' + country + file_identifier + '.hdf5')

# if os.path.isfile(fl_EXP_road):

#     EXP_road = LitPop.from_hdf5(fl_EXP_road)

# else:
#     from climada_petals.entity.exposures.openstreetmap.osm_dataloader import OSMRaw, OSMApiQuery, OSMFileQuery
#     from pathlib import Path
#     OSMRaw().get_data_geofabrik(country, file_format='pbf', save_path=project_folder)
#     MDGFileQuery = OSMFileQuery(Path(project_folder,'madagascar-latest.osm.pbf'))
#     gdf_roads = MDGFileQuery.retrieve_cis('road')
#     gdf_roads.set_crs(epsg=4326)
#     EXP_road = Exposures(gdf_roads)
#     EXP_road.gdf['impf_'] = 1
#     EXP_road.write_hdf5(fl_EXP_road)


if os.path.isfile(fl_EXP_edu):

    EXP_edu = LitPop.from_hdf5(fl_EXP_edu)

else:
    import geopandas as gpd
    import pandas as pd
    
    fl_OSM_school = '/Users/qinhan/Documents/ISF/UI_input/education/edu_OSM_point_VNM.shp'
    fl_GAR_school = '/Users/qinhan/Documents/ISF/UI_input/education/edu_NonOSM_point_VNM.shp'
    
    OSM_school = gpd.read_file(fl_OSM_school)
    
    GAR_school = gpd.read_file(fl_GAR_school)
    
    OSM_school = OSM_school[['geometry','Value']]
    
    GAR_school = GAR_school[['geometry','education']]
    
    OSM_school.rename(columns={'Value': 'value'}, inplace=True)
    
    GAR_school.rename(columns={'education': 'value'}, inplace=True)
    
    school = gpd.GeoDataFrame(pd.concat([OSM_school,GAR_school], ignore_index=True), crs=OSM_school.crs)
    
    EXP_edu = Exposures(school)
    
    print('\n' + '\x1b[1;03;30;30m' + 'EXP_school is now an Exposures:', str(type(EXP_edu)) + '\x1b[0m')
    
    EXP_edu.set_lat_lon()  # set latitude and longitude attributes from geometry
    
    EXP_edu.gdf['impf_'] = np.ones(EXP_edu.gdf.shape[0], int)  # provide impact functions for TC or any other peril
    
    EXP_edu.gdf['value'] = EXP_edu.gdf['value'] * 1000000 * 1.025**5
    
    EXP_edu.value_unit = '$'
    
    EXP_edu.write_hdf5(fl_EXP_edu)

#%%
if os.path.isfile(fl_EXP_health):

    EXP_health = LitPop.from_hdf5(fl_EXP_health)

else:
    import geopandas as gpd
    import pandas as pd
    
    fl_OSM_health = '/Users/qinhan/Documents/ISF/UI_input/health/health_OSM_point_VNM.shp'
    fl_GAR_health = '/Users/qinhan/Documents/ISF/UI_input/health/health_NonOSM_point_VNM.shp'
    
    OSM_health = gpd.read_file(fl_OSM_health)
    
    GAR_health = gpd.read_file(fl_GAR_health)

    OSM_health = OSM_health[['geometry','VALUE']]

    GAR_health = GAR_health[['geometry','health']]

    OSM_health.rename(columns={'VALUE': 'value'}, inplace=True)

    GAR_health.rename(columns={'health': 'value'}, inplace=True)
    
    health = gpd.GeoDataFrame(pd.concat([OSM_health,GAR_health], ignore_index=True), crs=OSM_health.crs)
    
    EXP_health = Exposures(health)
    
    print('\n' + '\x1b[1;03;30;30m' + 'EXP_health is now an Exposures:', str(type(EXP_health)) + '\x1b[0m')
    
    EXP_health.set_lat_lon()  # set latitude and longitude attributes from geometry
    
    EXP_health.gdf['impf_'] = np.ones(EXP_health.gdf.shape[0], int)  # provide impact functions for TC or any other peril
    
    EXP_health.gdf['value'] = EXP_health.gdf['value'] * 1000000 * 1.025**5
    
    EXP_health.value_unit = '$'
    
    EXP_health.write_hdf5(fl_EXP_health)
    
#%%


#%%

# imp_hist_tc_edu.aai_agg*100/sum(EXP_school.gdf.value)

#%%

# EXP_GA.plot_scatter()
from climada.engine import Impact
from climada.entity.impact_funcs import ImpactFuncSet, ImpfTropCyclone
from climada.entity.impact_funcs.trop_cyclone import ImpfSetTropCyclone
from climada.entity import ImpactFuncSet, ImpactFunc
# imp_fun = ImpfTropCyclone.from_emanuel_usa()
# imp_fun.plot();

# imp_fun_set_TC = ImpfSetTropCyclone.from_calibrated_regional_ImpfSet()

# imp_fun = imp_fun_set_TC.get_func(fun_id=6)

# imp_fun[0].id = 1

# imp_fun_set = ImpactFuncSet(imp_fun)

impf_ASIA = ImpactFunc(
    id=1,
    name="Flood Asia JRC Residential noPAA",
    intensity_unit='m/s',
    haz_type='TC',
    intensity=np.array([0, 10, 20, 30, 40, 50, 60, 70, 80]),
    mdd=np.array([0, 0, 0, 0, 0.05, 0.20, 0.45, 0.65, 0.78]),
    paa=np.array([1, 1, 1, 1, 1, 1, 1, 1, 1])
)
# impf_ASIA.plot()
imp_fun_set = ImpactFuncSet([impf_ASIA])
#%%

imp_hist_tc = Impact()
imp_hist_tc.calc(EXP_GA, imp_fun_set, Haz_tc_hist, save_mat=True)

freq_curve_hist_tc = imp_hist_tc.calc_freq_curve() # impact exceedance frequency curve
# freq_curve_hist_tc.plot()

print('Expected average annual historical impact: {:.3e} USD'.format(imp_hist_tc.aai_agg))
#%%

imp_hist_tc_edu = Impact()
imp_hist_tc_edu.calc(EXP_edu, imp_fun_set,Haz_tc_hist,save_mat=True)

freq_curve_hist_tc_edu = imp_hist_tc_edu.calc_freq_curve() # impact exceedance frequency curve
# freq_curve_hist_tc.plot()

print('Expected average annual historical impact on education: {:.3e} USD'.format(imp_hist_tc_edu.aai_agg))
#%%

imp_hist_tc_health = Impact()
imp_hist_tc_health.calc(EXP_health, imp_fun_set,Haz_tc_hist,save_mat=True)

freq_curve_hist_tc_health = imp_hist_tc_health.calc_freq_curve() # impact exceedance frequency curve
# freq_curve_hist_tc.plot()

print('Expected average annual historical impact on health: {:.3e} USD'.format(imp_hist_tc_health.aai_agg))



#%%
# fl_HAZ_ts = os.path.join(project_folder, 'HAZ_ts_hist_' + country + file_identifier + '.hdf5')

fl_HAZ_ts = '/Users/qinhan/Documents/IIASA/IIASA-ISF-ETH/CLIMADA_output/HAZ_ts_mix50_CDEM_VNM_v01.hdf5'

if os.path.isfile(fl_HAZ_ts):

    Haz_ts_hist = TCSurgeBathtub.from_hdf5(fl_HAZ_ts)
    
else:

    topo_path = input_folder + 'dem_compress.tif'
    
    topo_path2 = '/Users/qinhan/Documents/AXA/ECA_Vietnam/Exposure/Coastal_DEM/CDEM.tif'

    Haz_ts_hist = TCSurgeBathtub.from_tc_winds(Haz_tc_hist,topo_path2)

    Haz_ts_hist.haz_type = 'TCSurgeBathtub'

    Haz_ts_hist.write_hdf5(fl_HAZ_ts)


# Haz_ts_hist.plot_intensity(0)

#%%
fl_HAZ_ts45 = os.path.join(project_folder, 'HAZ_ts_futr45_' + country + file_identifier + '.hdf5')

if os.path.isfile(fl_HAZ_ts45):

    Haz_ts_futr45 = TCSurgeBathtub.from_hdf5(fl_HAZ_ts45)
    
else:

    topo_path = input_folder + 'output_SRTM15Plus.tif'
    
    topo_path2 = input_folder + 'elevation.tif'

    Haz_ts_futr45 = TCSurgeBathtub.from_tc_winds(Haz_tc_futr45,topo_path)

    Haz_ts_futr45.write_hdf5(fl_HAZ_ts45)
    

#%%
fl_HAZ_ts85 = os.path.join(project_folder, 'HAZ_ts_futr85_' + country + file_identifier + '.hdf5')

if os.path.isfile(fl_HAZ_ts85):

    Haz_ts_futr85 = TCSurgeBathtub.from_hdf5(fl_HAZ_ts85)
    
else:

    topo_path = input_folder + 'output_SRTM15Plus.tif'
    
    topo_path2 = input_folder + 'elevation.tif'

    Haz_ts_futr85 = TCSurgeBathtub.from_tc_winds(Haz_tc_futr85,topo_path)

    Haz_ts_futr85.write_hdf5(fl_HAZ_ts85)
    

#%%
from climada.engine import Impact
from climada.entity.impact_funcs import ImpactFuncSet, ImpfTropCyclone
from climada.entity.impact_funcs.trop_cyclone import ImpfSetTropCyclone
from climada_petals.entity.impact_funcs.river_flood import flood_imp_func_set
from climada.entity import ImpactFuncSet, ImpactFunc


impf_ASIA = ImpactFunc(
    id=1,
    name="Flood Asia JRC Residential noPAA",
    intensity_unit='m',
    haz_type='TCSurgeBathtub',
    intensity=np.array([0, 0.5, 1., 1.5, 2., 3., 4., 5., 6., 12.]),
    mdd=np.array([0.000, 0.3266, 0.4941, 0.6166, 0.7207, 0.8695,
                         0.9315, 0.9836, 1.0000, 1.0000]),
    paa=np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
)
# impf_ASIA.plot()
impf_set = ImpactFuncSet([impf_ASIA])

imp_hist_ts=Impact()

imp_hist_ts.calc(EXP_GA,impf_set,Haz_ts_hist,save_mat=True)
# imp_hist_ts.plot_scatter_eai_exposure()


freq_curve_hist = imp_hist_ts.calc_freq_curve() # impact exceedance frequency curve
# freq_curve_hist.plot()

print('Expected average annual historical impact: {:.3e} USD'.format(imp_hist_ts.aai_agg))


#%%plotting



#%% Plotting section (run after the above)
import matplotlib.font_manager as font_manager
import matplotlib.colors as colors
import matplotlib.pyplot as plt

upper = plt.cm.Blues(np.arange(256))
lower = np.ones((1,4))
# combine parts of colormap
cmap = np.vstack(( lower, upper ))
# convert to matplotlib colormap
cmap = colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

upper = plt.cm.Reds(np.arange(256))
cmap_red = np.vstack(( lower, upper ))
# convert to matplotlib colormap
cmap_red = colors.ListedColormap(cmap_red, name='myColorMap', N=cmap_red.shape[0])

fontpath = '/System/Library/Fonts/Optima.ttc'

prop = font_manager.FontProperties(fname=fontpath)

plt.rcParams["font.weight"] = "bold"

plt.rcParams["axes.labelweight"] = "bold"

# plt.rcParams.update(
#     {'font.size': 15, 'axes.titlesize': 15, 'legend.fontsize': 15, 'font.family': [prop.get_name(), 'sans-serif']})

plt.rcParams.update(
        {'font.size': 12, 'axes.titlesize': 16, 'legend.fontsize': 16, 'font.family': [prop.get_name(), 'sans-serif']})

#%%

norm = colors.LogNorm(vmin=100, vmax=100000000)

ax = EXP_GA.plot_hexbin(norm=norm, pop_name=False, cmap='RdBu_r', buffer=1)

    # norm = colors.LogNorm(vmin=minmax[0], vmax=minmax[1])

    # ax = EXP.plot_hexbin(pop_name=False, cmap='RdBu_r', buffer=1)

ax.set_title('')

fig = plt.gcf()

fig.set_size_inches(4.5, 8)

plt.subplots_adjust(left=0.15, right=0.85, top=0.95, bottom=0.05)

cbar = fig.axes[1]

cbar.set_ylabel('Asset value ($) per km\u00b2')

#%%

norm = colors.LogNorm(vmin=100, vmax=100000000)

ax = EXP_edu.plot_hexbin(norm=norm, pop_name=False, cmap='RdBu_r', buffer=1)

    # norm = colors.LogNorm(vmin=minmax[0], vmax=minmax[1])

    # ax = EXP.plot_hexbin(pop_name=False, cmap='RdBu_r', buffer=1)

ax.set_title('')

fig = plt.gcf()

fig.set_size_inches(4.5, 8)

plt.subplots_adjust(left=0.15, right=0.85, top=0.95, bottom=0.05)

cbar = fig.axes[1]

cbar.set_ylabel('Asset value ($) per km\u00b2')

#%%

norm = colors.LogNorm(vmin=100, vmax=100000000)

ax = EXP_edu.plot_hexbin(norm=norm, pop_name=False, cmap='RdBu_r', buffer=1)

    # norm = colors.LogNorm(vmin=minmax[0], vmax=minmax[1])

    # ax = EXP.plot_hexbin(pop_name=False, cmap='RdBu_r', buffer=1)

ax.set_title('')

fig = plt.gcf()

fig.set_size_inches(4.5, 8)

plt.subplots_adjust(left=0.15, right=0.85, top=0.95, bottom=0.05)

cbar = fig.axes[1]

cbar.set_ylabel('Asset value ($) per km\u00b2')

#%%

norm = colors.LogNorm(vmin=100, vmax=100000000)

ax = EXP_health.plot_hexbin(norm=norm, pop_name=False, cmap='RdBu_r', buffer=1)

    # norm = colors.LogNorm(vmin=minmax[0], vmax=minmax[1])

    # ax = EXP.plot_hexbin(pop_name=False, cmap='RdBu_r', buffer=1)

ax.set_title('')

fig = plt.gcf()

fig.set_size_inches(4.5, 8)

plt.subplots_adjust(left=0.15, right=0.85, top=0.95, bottom=0.05)

cbar = fig.axes[1]

cbar.set_ylabel('Asset value ($) per km\u00b2')

#%% Hazard tc hist

norm = colors.LogNorm(vmin=0.000000001, vmax=0.0001)

ax = imp_hist_tc.plot_hexbin_eai_exposure(pop_name=False, norm=norm, cmap=cmap, buffer=1, ignore_zero=False)

ax.set_title('')

fig = plt.gcf()

fig.set_size_inches(4.5, 8)

plt.subplots_adjust(left=0.15, right=0.85, top=0.95, bottom=0.05)

cbar = fig.axes[1]

cbar.set_ylabel('Damage ($) per km\u00b2')
#%% Don't run
from matplotlib.colors import SymLogNorm


diff85 = imp_futr_tc85.eai_exp - imp_hist_tc.eai_exp

imp_hist_tc_his.eai_exp = diff85

from  matplotlib.colors import LinearSegmentedColormap
cmap=LinearSegmentedColormap.from_list('rg',["g", "w", "r"], N=256) 


ax = imp_hist_tc_his.plot_hexbin_eai_exposure(pop_name=False, 
                                              norm=SymLogNorm(linthresh=2, linscale=2, vmin=-1000, vmax=10000), 
                                              cmap=cmap, buffer=1, ignore_zero=False)

ax.set_title('')

fig = plt.gcf()

fig.set_size_inches(4.5, 8)

plt.subplots_adjust(left=0.15, right=0.85, top=0.95, bottom=0.05)

cbar = fig.axes[1]

cbar.set_ylabel('Damage ($) per km\u00b2')


#%%
norm = colors.LogNorm(vmin=1, vmax=100000)

ax = imp_hist_tc_edu.plot_hexbin_eai_exposure(pop_name=False, norm=norm, cmap=cmap, buffer=1, ignore_zero=False)

ax.set_title('')

fig = plt.gcf()

fig.set_size_inches(4.5, 8)

plt.subplots_adjust(left=0.15, right=0.85, top=0.95, bottom=0.05)

cbar = fig.axes[1]

cbar.set_ylabel('Damage ($) per km\u00b2')

norm = colors.LogNorm(vmin=1, vmax=100000)

ax = imp_hist_tc_health.plot_hexbin_eai_exposure(pop_name=False, norm=norm, cmap=cmap, buffer=1, ignore_zero=False)

ax.set_title('')

fig = plt.gcf()

fig.set_size_inches(4.5, 8)

plt.subplots_adjust(left=0.15, right=0.85, top=0.95, bottom=0.05)

cbar = fig.axes[1]

cbar.set_ylabel('Damage ($) per km\u00b2')



#%% Hazard flood hist

norm = colors.LogNorm(vmin=1, vmax=100)

ax = imp_hist_ts.plot_hexbin_eai_exposure(pop_name=False, norm=norm, cmap=cmap, buffer=1, ignore_zero=False)

ax.set_title('')

fig = plt.gcf()

fig.set_size_inches(4.5, 8)

plt.subplots_adjust(left=0.15, right=0.85, top=0.95, bottom=0.05)

cbar = fig.axes[1]

cbar.set_ylabel('Damage ($) per km\u00b2')

#%%

