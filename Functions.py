import os
from climada_petals.hazard.river_flood import RiverFlood
from climada.hazard.trop_cyclone import TropCyclone
from climada.hazard.tc_tracks import TCTracks
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gamma
from scipy.stats import weibull_min
from sklearn.metrics import r2_score
from scipy.integrate import simps

project_folder = '/Users/qinhan/Documents/IIASA/IIASA-ISF-ETH/GSSP/UltimateOutputs/'  # results and data will be saved in this folder
input_folder = '/Users/qinhan/Documents/IIASA/IIASA-ISF-ETH/GSSP/Inputs/'  # this folder should be the folder containing the input files
hazard_folder = input_folder + 'pakistan_flood/'
country = 'PAK'
basin = 'NI'
hist_years = list(range(1975, 2006))
future_years = list(range(2020, 2051))
number_rdw = 5
horizon = 2050
file_identifier = '_v01'
fl_haz_tc_cc0 = os.path.join(project_folder, 'HAZ_tc_hist_' + country + file_identifier + '.hdf5')

# class Object(object):
#     pass
# hist_storms = Object()
# hist_storms.lon_min = 60.8
# hist_storms.lat_min = 23.6
# hist_storms.lon_max = 80.3
# hist_storms.lat_max = 37.3
# hist_storms.start_year = 1983
# hist_storms.end_year = 2023

def tracks_ID(storm_range, basin):

    lon_min = storm_range.lon_min
    lat_min = storm_range.lat_min
    lon_max = storm_range.lon_max
    lat_max = storm_range.lat_max
    start_year = storm_range.start_year
    end_year = storm_range.end_year

    box = np.array([[lat_min, lat_max], [lon_min, lon_max]])
    #
    year_range = [start_year,end_year]

    # reading ibtracs (csv version)
    ibtracs = pd.read_csv('/Users/qinhan/Documents/IIASA/IIASA-ISF-ETH/GSSP/Inputs/ibtracs.ALL.list.v04r00.csv', skiprows=[1], delimiter=',',
                          usecols=['SID','BASIN', 'NAME', 'ISO_TIME', 'LAT', 'LON', 'USA_WIND'])
    ibtracs[['LAT', 'LON', 'USA_WIND']] = ibtracs[['LAT', 'LON', 'USA_WIND']].apply(pd.to_numeric, errors='coerce')
    ibtracs['ISO_TIME'] = pd.to_datetime(ibtracs['ISO_TIME'])

    # just setting time to the right format to perform comparisons (simpler version is to use int(ibtracs.SID[:4]) which is exactly the year of the event)
    time_bounds = (np.datetime64('{}-01-01T00:00:00.000000000'.format(year_range[0])),
                   np.datetime64('{}-01-01T00:00:00.000000000'.format(year_range[1] + 1)))

    # filters
    basin = ibtracs.BASIN == basin

    region = (ibtracs.LAT >= box[0, 0]) & (ibtracs.LAT <= box[0, 1]) & (ibtracs.LON >= box[1, 0]) & (
                ibtracs.LON <= box[1, 1])

    period = (ibtracs.ISO_TIME >= time_bounds[0]) & (ibtracs.ISO_TIME < time_bounds[1])

    # intensity = (ibtracs.USA_WIND >= min_wind)

    # gathering storm ids
    selection = region & period

    ID = np.unique(ibtracs.SID[selection])

    return list(ID)

def intro_haz_tc(fl_haz_tc, cc):

    haz_tc = TropCyclone()

    if os.path.isfile(fl_haz_tc):

        haz_tc.from_hdf5(fl_haz_tc)

    elif cc == 0:

        sel_ibtracs = TCTracks()

        ID = tracks_ID()

        sel_ibtracs.from_ibtracs_netcdf(storm_id=ID)

        print('tracks selected')

        sel_ibtracs.calc_perturbed_trajectories()

        print('probabilistic events generated')

        sel_ibtracs.equal_timestep()

        print('equal timestep')

        # centr = centr_gen(fl_centroids)

        haz_tc.set_from_tracks(sel_ibtracs)

        print('HAZ finished')

        haz_tc.check()

        haz_tc.write_hdf5(fl_haz_tc)

    elif cc == 45:

        haz_tc_cc0 = intro_haz_tc(fl_haz_tc_cc0, 0)

        haz_tc = haz_tc_cc0.set_climate_scenario_knu(ref_year=horizon,rcp_scenario=45)

        haz_tc.write_hdf5(fl_haz_tc)

    elif cc == 85:

        haz_tc_cc0 = intro_haz_tc(fl_haz_tc_cc0, 0)

        haz_tc = haz_tc_cc0.set_climate_scenario_knu(ref_year=horizon, rcp_scenario=85)

        haz_tc.write_hdf5(fl_haz_tc)

    return haz_tc

def calc_annually_aggregated_impact(imp,total_asset_value):
    
    year_index = []

    for i in imp.event_name:
        
        year_index.append(str(i[0:4])+i[14:])

    # df = pd.DataFrame({'year':year_index,'impact':imp_hist_tc.at_event,'impact_public':imp_hist_tc_health.at_event+imp_hist_tc_edu.at_event})

    df = pd.DataFrame({'year':year_index,'impact':imp.at_event})

    # df = pd.DataFrame({'year':year_index,'impact':imp_hist_tc_health.at_event+imp_hist_tc_edu.at_event})

    year_no_duplicate = df.year.drop_duplicates().reset_index(drop=True)


    imp_year = []

    for i in range(len(year_no_duplicate)):
        year_i = year_no_duplicate[i]
        boolean_year = df.year == year_i
        impact_year_i = sum(df.impact[boolean_year])
        imp_year.append(impact_year_i)


    df_final = pd.DataFrame({'year':year_no_duplicate,'impact_year':imp_year})

    boolean_zero = df_final['impact_year']>0

    df_final_non_zero = df_final[boolean_zero]


    df_final_non_zero = df_final_non_zero.sort_values(by=['impact_year'])

    df_final_non_zero['impact_year'] = df_final_non_zero['impact_year']/total_asset_value

    df_final_non_zero['frequency'] = [1/len(df_final_non_zero)]*len(df_final_non_zero)

    cumsum = df_final_non_zero['frequency'].cumsum()

    df_final_non_zero['cumsum'] = cumsum

    return df_final_non_zero

def gamma_CDF(damage_rate, return_period):

# Fit gamma cumulative distribution
    params = gamma.fit(damage_rate, floc=0)  # Assuming your X coordinates are positive
    shape, loc, scale = params
    
    # Generate x values for plotting
    x_values = np.linspace(0, max(damage_rate), 1000)
    
    # Calculate y values using gamma cumulative distribution function
    y_values = gamma.cdf(x_values, shape, loc=loc, scale=scale)
    
    y_predicted = gamma.cdf(damage_rate, shape, loc=loc, scale=scale)
    r_squared = r2_score(return_period, y_predicted)   
    # Plot data and fitted curve
    plt.plot(damage_rate, return_period, 'o', label='Data')
    plt.plot(x_values, y_values, label='Fitted Gamma CDF')
    
    plt.xlabel('X')
    plt.ylabel('Cumulative Probability (Y)')
    plt.title(f'Fitting Gamma Cumulative Distribution to Data\nR-squared: {r_squared:.2f}')
    plt.legend()
    plt.grid(True)
    plt.show()
    
    return params
    
def weibull_CDF(damage_rate, return_period):
# Fit weibull cumulative distribution
    params = weibull_min.fit(damage_rate, floc=0)  # Assuming your X coordinates are positive
    c, loc, scale = params

    # Generate x values for plotting
    x_values = np.linspace(0, max(damage_rate), 1000)

    # Calculate y values using Weibull cumulative distribution function
    y_values = weibull_min.cdf(x_values, c, loc=loc, scale=scale)

    y_predicted = weibull_min.cdf(damage_rate, c, loc=loc, scale=scale)
    r_squared = r2_score(return_period, y_predicted)

    # Plot data and fitted curve
    plt.plot(damage_rate, return_period, 'o', label='Data')
    plt.plot(x_values, y_values, label='Fitted Weibull CDF')
    
    plt.xlabel('Percentage of total asset damages (%)')
    plt.ylabel('Non-exceedence frequency')
    plt.title(f'Fitting Weibull Cumulative Distribution to Data\nR-squared: {r_squared:.2f}')
    plt.legend()
    plt.grid(True)
    plt.show()
    
    return params

def exponential_CDF(damage_rate, return_period):
    # Define the model function
    def model_function(x, b):
        return 1 - np.exp(-b * x)
    
    
    integral = simps(1-return_period, damage_rate)
    

    # Generate x values for plotting
    x_values = np.linspace(0, max(damage_rate), 1000)
    
    # Calculate y values using the model function
    y_values = model_function(x_values, 1/integral)
    
    y_predicted = model_function(damage_rate, 1/integral)
    r_squared = r2_score(return_period, y_predicted)
    
    # Plot data and fitted curve
    plt.plot(damage_rate, return_period, 'o', label='Data')
    plt.plot(x_values, y_values, label='Fitted Exponential CDF')
    
    plt.xlabel('X')
    plt.ylabel('Cumulative Probability (Y)')
    plt.title(f'Fitting 1 - exp(-bx) Cumulative Distribution to Data\nR-squared: {r_squared:.2f}')
    plt.legend()
    plt.grid(True)
    plt.show()
    return 1/integral



def intro_haz_riverflood(fl_Haz_riverflood,time):

    if os.path.isfile(fl_Haz_riverflood):
        HAZ_ens = RiverFlood.from_hdf5(fl_Haz_riverflood)
    else:

        fl_final_rcp6085_haz = os.path.join(input_folder, 'selectedmodels')

        with open(fl_final_rcp6085_haz) as f:
            final_rcp6085_haz = f.readlines()

        final_rcp6085_haz = [x.strip() for x in final_rcp6085_haz]

        if time == 'historical':
            years = hist_years
        else:
            years = future_years

        i = 0

        for fl_HAZ_rf in final_rcp6085_haz:
            if time == 'historical':
                fl_HAZ_rf = fl_HAZ_rf.replace('rcp60_2005soc', 'historical_histsoc')
            else:
                fl_HAZ_rf = fl_HAZ_rf.replace('rcp60_2005soc', 'rcp85_2005soc')

            fl_HAZ_rf_dph = os.path.join(hazard_folder + 'flddph_' + fl_HAZ_rf)

            fl_HAZ_rf_frc = os.path.join(hazard_folder + 'fldfrc_' + fl_HAZ_rf)

            if os.path.isfile(fl_HAZ_rf_dph):

                i += 1

                HAZ_rf = RiverFlood.from_nc(dph_path=fl_HAZ_rf_dph, frc_path=fl_HAZ_rf_frc, years=years)

                if i == 1:
                    temp_frc = HAZ_rf.fraction

                    temp_int = HAZ_rf.intensity

                else:
                    temp_frc += HAZ_rf.fraction

                    temp_int += HAZ_rf.intensity

        import copy

        HAZ_ens = copy.deepcopy(HAZ_rf)

        frac = temp_frc / i
        inte = temp_int / i

        HAZ_ens.fraction = frac
        HAZ_ens.intensity = inte

        # HAZ_ens_hist.plot_intensity(event=0, smooth = False)

        HAZ_ens.write_hdf5(fl_Haz_riverflood)

        # HAZ_ens.set_flooded_area(save_centr=True)
        #
        # HAZ_ens.set_flood_volume(save_centr=True)

    # HAZ_ens.write_hdf5(fl_HAZ_city)

    return HAZ_ens
