# from climada_petals.hazard.river_flood import RiverFlood
from Functions import intro_haz_riverflood, tracks_ID
from climada.entity import Exposures, LitPop
from climada_petals.hazard import TCSurgeBathtub
import os
import numpy as np

import time


project_folder = '/Users/qinhan/Documents/IIASA/IIASA-ISF-ETH/GSSP/UltimateOutputs/'  # results and data will be saved in this folder
input_folder = '/Users/qinhan/Documents/IIASA/IIASA-ISF-ETH/GSSP/Inputs/'  # this folder should be the folder containing the input files
# hazard_folder = input_folder + 'pakistan_flood/'
file_identifier = '_v01'
country = 'MDG'
#%%
class Object(object):
    pass


hist_storms = Object()
hist_storms.lon_min = 43.26
hist_storms.lat_min = -25.68
hist_storms.lon_max = 50.54
hist_storms.lat_max = -11.86
hist_storms.start_year = 1983
hist_storms.end_year = 2023
horizon = 2050


EXP_norm = LitPop.from_countries(countries=country, res_arcsec=30, reference_year=2020,fin_mode='norm')

EXP_norm_gdf = EXP_norm.gdf

quantile90 = EXP_norm_gdf['value'].quantile(0.99)

value_boolean = EXP_norm_gdf['value'] >= quantile90

lon_value = EXP_norm_gdf.centroid.x[value_boolean].to_numpy()

lat_value = EXP_norm_gdf.centroid.y[value_boolean].to_numpy()

from climada.hazard import Centroids

min_lat, max_lat, min_lon, max_lon = hist_storms.lat_min, hist_storms.lat_max, hist_storms.lon_min, hist_storms.lon_max

as30cent = Centroids.from_pnt_bounds((min_lon, min_lat, max_lon, max_lat), res=0.0083333)

as30cent.set_on_land()

boolean_land = as30cent.on_land

lat_land = as30cent.lat[boolean_land]

lon_land = as30cent.lon[boolean_land]

landas30cent = Centroids.from_lat_lon(lat_land, lon_land)

# landas30cent.set_dist_coast()

# boolean_coast = landas30cent.dist_coast <= 50000

landas30cent.set_elevation('/Users/qinhan/Documents/IIASA/IIASA-ISF-ETH/GSSP/Inputs/output_SRTM15Plus.tif')

boolean_coast = landas30cent.elevation <= 20

lat_coast = landas30cent.lat[boolean_coast]

lon_coast = landas30cent.lon[boolean_coast]

# coast_centr=Centroids.from_lat_lon(lat_coast,lon_coast)

# coast_centr.set_dist_coast()

# coast_centr.plot()

as150cent = Centroids.from_pnt_bounds((min_lon, min_lat, max_lon, max_lat), res=0.0416667)

as150cent.set_on_land()

boolean_land150 = as150cent.on_land

lat_land150 = as150cent.lat[boolean_land150]

lon_land150 = as150cent.lon[boolean_land150]


lat_final = np.concatenate((lat_coast,lat_value,lat_land150),axis=None)

lon_final = np.concatenate((lon_coast,lon_value,lon_land150),axis=None)

final_centr = Centroids.from_lat_lon(lat_final, lon_final)

# final_centr.plot()



ID_tracks_MDG = tracks_ID(hist_storms,'SI')
#%%

from climada.hazard.trop_cyclone import TropCyclone
from climada.hazard.tc_tracks import TCTracks

fl_HAZ_tc = os.path.join(project_folder, 'HAZ_tc_hist_' + country + file_identifier + '.hdf5')

# fl_HAZ_tc = '/Users/qinhan/Documents/IIASA/IIASA-ISF-ETH/CLIMADA_output/HAZ_tc_mix50_final_MDG_v01.hdf5'

if os.path.isfile(fl_HAZ_tc):

    Haz_tc_hist = TropCyclone.from_hdf5(fl_HAZ_tc)
    
else:
    
    tracks_MDG = TCTracks.from_ibtracs_netcdf(storm_id=ID_tracks_MDG)
    print('tracks selected')
    tracks_MDG.equal_timestep()
    print('equal timestep')
    tracks_MDG.calc_perturbed_trajectories(nb_synth_tracks=24)
    print('probabilistic events generated')
    start = time.time()
    Haz_tc_hist = TropCyclone.from_tracks(tracks_MDG, centroids=final_centr)
    print('HAZ finished')
    end = time.time()
    print(end-start)
    Haz_tc_hist.check()

    Haz_tc_hist.write_hdf5(fl_HAZ_tc)
    
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

# sum(EXP_GA.gdf['value'])

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
    from climada_petals.entity.exposures.openstreetmap.osm_dataloader import OSMRaw, OSMApiQuery, OSMFileQuery
    from pathlib import Path
    OSMRaw().get_data_geofabrik(country, file_format='pbf', save_path=project_folder)
    MDGFileQuery = OSMFileQuery(Path(project_folder,'madagascar-latest.osm.pbf'))    
    gdf_education = MDGFileQuery.retrieve_cis('education')
    
    boolean_education_polygon = gdf_education.geometry.type == 'MultiPolygon'
    
    gdf_education = gdf_education[boolean_education_polygon].reset_index()
    
    education_area = gdf_education.set_crs(epsg=4326).to_crs(epsg=3857).area
    
    gdf_education.geometry = gdf_education.centroid
    
    EXP_education = Exposures(gdf_education)
    
    EXP_education.gdf = EXP_education.gdf.to_crs(4326)
    
    EXP_education.gdf['longitude'] = EXP_education.gdf.centroid.x
    
    EXP_education.gdf['latitude'] = EXP_education.gdf.centroid.y
    
    EXP_education.gdf['impf_'] = 1
    
    EXP_education.gdf['value'] = education_area*250
    
    EXP_education.write_hdf5(fl_EXP_edu)
    
if os.path.isfile(fl_EXP_health):
    EXP_health = LitPop.from_hdf5(fl_EXP_health)


else:
    from climada_petals.entity.exposures.openstreetmap.osm_dataloader import OSMRaw, OSMApiQuery, OSMFileQuery
    from pathlib import Path
    OSMRaw().get_data_geofabrik(country, file_format='pbf', save_path=project_folder)
    MDGFileQuery = OSMFileQuery(Path(project_folder,'madagascar-latest.osm.pbf'))    
    gdf_health = MDGFileQuery.retrieve_cis('healthcare')
    
    boolean_health_polygon = gdf_health.geometry.type == 'MultiPolygon'
    
    gdf_health = gdf_health[boolean_health_polygon].reset_index()
    
    health_area = gdf_health.set_crs(epsg=4326).to_crs(epsg=3857).area
    
    gdf_health.geometry = gdf_health.centroid
    
    EXP_health = Exposures(gdf_health)
    
    EXP_health.gdf = EXP_health.gdf.to_crs(4326)
    
    EXP_health.gdf['longitude'] = EXP_health.gdf.centroid.x
    
    EXP_health.gdf['latitude'] = EXP_health.gdf.centroid.y
    
    EXP_health.gdf['impf_'] = 1
    
    EXP_health.gdf['value'] = health_area*250
    
    EXP_health.write_hdf5(fl_EXP_health)

    
#%%

# EXP_GA.plot_scatter()
from climada.engine import Impact
from climada.entity.impact_funcs import ImpactFuncSet, ImpfTropCyclone
from climada.entity.impact_funcs.trop_cyclone import ImpfSetTropCyclone
# imp_fun = ImpfTropCyclone.from_emanuel_usa()
# imp_fun.plot();

imp_fun_set_TC = ImpfSetTropCyclone.from_calibrated_regional_ImpfSet()

imp_fun = imp_fun_set_TC.get_func(fun_id=6)

imp_fun[0].id = 1

imp_fun[0].paa = imp_fun[0].paa

imp_fun_set = ImpactFuncSet(imp_fun)
#%%
imp_hist_tc = Impact()
imp_hist_tc.calc(EXP_GA, imp_fun_set,Haz_tc_hist,save_mat=True)

freq_curve_hist_tc = imp_hist_tc.calc_freq_curve() # impact exceedance frequency curve
# freq_curve_hist_tc.plot()

print('Expected average annual historical impact: {:.3e} USD'.format(imp_hist_tc.aai_agg))
#%%

imp_hist_tc_edu = Impact()
imp_hist_tc_edu.calc(EXP_edu, imp_fun_set,Haz_tc_hist,save_mat=True)

freq_curve_hist_tc_edu = imp_hist_tc_edu.calc_freq_curve() # impact exceedance frequency curve
# freq_curve_hist_tc.plot()

print('Expected average annual historical impact on education: {:.3e} USD'.format(imp_hist_tc_edu.aai_agg))



imp_hist_tc_health = Impact()
imp_hist_tc_health.calc(EXP_health, imp_fun_set,Haz_tc_hist,save_mat=True)

freq_curve_hist_tc_health = imp_hist_tc_health.calc_freq_curve() # impact exceedance frequency curve
# freq_curve_hist_tc.plot()

print('Expected average annual historical impact on health: {:.3e} USD'.format(imp_hist_tc_health.aai_agg))



#%%
fl_HAZ_ts = os.path.join(project_folder, 'HAZ_ts_hist_' + country + file_identifier + '.hdf5')

if os.path.isfile(fl_HAZ_ts):

    Haz_ts_hist = TCSurgeBathtub.from_hdf5(fl_HAZ_ts)
    
else:

    topo_path = input_folder + 'output_SRTM15Plus.tif'
    
    topo_path2 = input_folder + 'elevation.tif'

    Haz_ts_hist = TCSurgeBathtub.from_tc_winds(Haz_tc_hist,topo_path)

    Haz_ts_hist.write_hdf5(fl_HAZ_ts)
    

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

norm = colors.LogNorm(vmin=1, vmax=100000)

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

