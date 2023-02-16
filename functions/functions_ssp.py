def ssp_read_data(filename_suffix,varlist,input_dimension, output_dimension,**kw):
    '''
    read data under different ssp scenarios
    in input: everything the same with read_nc_data_universal, except input filename_suffix
        instead of filename. For example, 'QRS_TOT_QRL_TOT_coarse_4d'
    output: 
        a dictionary of data under different scenarios and different runs
    '''
    from functions_io import read_nc_data_universal
    scenarios = {'SSP1-2.6':1,
            'SSP2-4.5':5,
            'SSP3-7.0':1,
            'SSP5-8.5':5}
    result_scenarios = {}
    for scenario in scenarios:
        result_runs = {}
        for i_run in range(scenarios[scenario]):
            run_num = i_run+1
            filename = '/glade/work/yuwandi/cmip6-waccm/{}-WACCM.00{}_{}.nc'\
                         .format(str(scenario).zfill(2),run_num,filename_suffix)
            result_runs['00{}'.format(run_num)] = read_nc_data_universal(filename,
                                        varlist,input_dimension, output_dimension,**kw)
        result_scenarios[scenario] = result_runs
    return result_scenarios


def ssp_mvr_trend_model(filename_suffix,varlist,input_dimension, output_dimension,**kw):
    '''
    loop through different runs of differnt ssps, call read_nc_data_universal,read_index
            and mvr_trend_model functions
    input: 
        filename_suffix:the suffix of a filename
            For example, 'QRS_TOT_QRL_TOT_coarse_4d'
        varlist: a list of variable names 
        input_dimension: input dimension, 1,2,3, or 4
        output_dimension: output dimension, 1,2,3, or 4
        **kw:
            lat_range_n: range of averaging latitude
            lat_range_s: range of averaging latitude,lat_range_s<lat_range_n
            lev_high: higher bound of pressure
            lev_low: lower bound of pressure 
                    if output_dimension == 1, use lev_high == lev_low to specify a level for 1-d output
    '''
    from functions_io import read_nc_data_universal
    from functions_models import mvr_trend_model
    from functions_statistics import read_index
    scenarios = {'SSP1-2.6':1,
            'SSP2-4.5':5,
            'SSP3-7.0':1,
            'SSP5-8.5':5}
    result_scenarios = {}
    for scenario in scenarios:
        result_runs = {}
        for i_run in range(scenarios[scenario]):
            run_num = i_run+1
            filename = '/glade/work/yuwandi/cmip6-waccm/{}-WACCM.00{}_{}.nc'\
                         .format(str(scenario).zfill(2),run_num,filename_suffix)
            var = read_nc_data_universal(filename,
                    varlist,input_dimension, output_dimension,**kw)
            index = read_index(scenario,run_num)
            result_vars = {}
            for var_name in varlist:
                result_vars[var_name] = mvr_trend_model(var[var_name],index,**kw)
            result_runs['00{}'.format(run_num)] = result_vars
        result_scenarios[scenario] = result_runs
    return result_scenarios

def ssp_run_mean(scenario_in):
    '''
    calculate the mean value of each scenario 
    input: 
        scenario_in: a dictionary of scenarios, containing a dictionary of runs 
    e.g.: 
    scenario_output = ssp_read_data(filename_suffix,varlist,input_dimension, output_dimension)
        var_out = ssp_run_mean(scenario_output)
    output:
        a dict of scenario mean 
    '''
    import numpy as np
    scenarios = {'SSP1-2.6':1,
            'SSP2-4.5':5,
            'SSP3-7.0':1,
            'SSP5-8.5':5}
    scenario_out = {}
    for scenario in scenarios:
        if type(scenario_in[scenario]['001']) == np.ndarray:
            var_out = np.zeros([scenarios[scenario],*(scenario_in[scenario]['001'].shape)])
            for ii in range(scenarios[scenario]):
                var_out[ii] = scenario_in[scenario]['00{}'.format(1+ii)]
            scenario_out[scenario] = np.nanmean(var_out,axis=0)
        elif type(scenario_in[scenario]['00{}'.format(1)]) == dict:
            run_dict = {}
            for var in scenario_in[scenario]['001']:
                if type(scenario_in[scenario]['001'][var]) == np.ndarray:
                    var_out = np.zeros([scenarios[scenario],*(scenario_in[scenario]['001'][var].shape)])
                    for ii in range(scenarios[scenario]):
                        var_out[ii] = scenario_in[scenario]['00{}'.format(1+ii)][var]
                    run_dict[var] = np.nanmean(var_out,axis=0)
                elif type(scenario_in[scenario]['001'][var]) == dict:
                    run_subdict = {}
                    for var_sub in scenario_in[scenario]['001'][var]:
                        var_out = np.zeros([scenarios[scenario],*(scenario_in[scenario]['001'][var][var_sub].shape)])
                        for ii in range(scenarios[scenario]):
                            var_out[ii] = scenario_in[scenario]['00{}'.format(1+ii)][var][var_sub]
                        run_subdict[var_sub] = np.nanmean(var_out,axis=0)
                    run_dict[var] = run_subdict
        scenario_out[scenario] = run_dict
    return scenario_out
                