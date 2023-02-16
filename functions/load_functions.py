import sys
sys.path.insert(0, '/glade/u/home/yuwandi/functions/')
from functions_statistics import mvr_filter, corr_nan, linear_regress, detrend

from functions_io import read_nc_data_universal, save_to_2d, save_to_3d, save_to_4d, add_dimension_to_filename, read_mls_h2o, read_mls_t, read_obs_index, read_cmip6_index

from functions_plot import add_order,print_lat, plot_p, get_color, add_z_axis, mon_arr

from functions_models import mvr_trend_model, pmc_0d_pressure_model

from functions_basics import geo_avg, cal_anomaly, select_month, normalize, find_index,moving_average,convert_cftime_to_int,get_lat_lim

from functions_physics import SWV,cold_point_tropopause, number_density, get_cooling_power, power_spectrum,get_dtdy_term

from functions_ssp import ssp_read_data,ssp_mvr_trend_model, ssp_run_mean

from functions_grid import GRID_3d, GRID_2d, GRID_3d_sameloc