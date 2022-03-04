"""
This script contains multiple methods aiming to parse SCONE .par files, to modify
them and apply them on a controller file.
"""
import pandas as pd
import scone_recursive_parser as srp


def apply_par_file_to_series_with_new_range(par_file, series, min_max_ratio):
    '''Extracts parameter values from a par file and applies them on
    a controller PANDAS Series, with new min, max (but same std).
    Parameters
    ----------
    par_file: (string) path to par file.
    series: (PANDAS Series) recursive series containing controller data.
    min_max_ratio: (float) ratio such that new max/min = par_value*(1 +/- min_max_ratio*sign(par_value)).

    Returns
    -------
    series_out: (PANDAS Series) controller series with new values extracted from par files.
    '''
    par_df = import_par_df(par_file)
    par_df = create_new_parameter_distribution(par_df, min_max_ratio)

    key_df = generate_key_df(par_df)
    threshold_df = generate_threshold_df(par_df)
    series_out = srp.series_deep_copy(series)

    key_df.apply(lambda x: apply_par_key_to_series(series_out, x['key_name'],
                                                   x['new_range_value'], x['states'],
                                                   x['target'], x['source']), axis=1)

    threshold_df.apply(lambda x: srp.set_unique_key(series_out, x['parameter'],
                                                    x['new_range_value']), axis=1)

    return series_out


def apply_par_file_to_series_with_fixed_values(par_file, series):
    '''Extracts parameter values from a par file and applies them on
    a controller PANDAS Series, without allowing parameter variation on the 
    new controller file.
    Parameters
    ----------
    par_file: (string) path to par file.
    series: (PANDAS Series) recursive series containing controller data.

    Returns
    -------
    series_out: (PANDAS Series) controller series with new values extracted from par files.
    '''
    par_df = import_par_df(par_file)
    
    key_df = generate_key_df(par_df)
    threshold_df = generate_threshold_df(par_df)
    series_out = srp.series_deep_copy(series)
    
    key_df.apply(lambda x: apply_par_key_to_series(series_out, x['key_name'],
                                                x['value'], x['states'], 
                                                x['target'], x['source']), axis=1)    

    threshold_df.apply(lambda x: srp.set_unique_key(series_out, x['parameter'], 
                                                    x['value']), axis =1)
    
    return series_out


def import_par_df(par_file):
    '''Parses a SCONE parameter files to a PANDAS DataFrame
    Parameters
    ----------
    par_file: (string) path to par file.

    Returns
    -------
    par_df: (PANDAS DataFrame) DataFrame containing parsed par file.
    '''
    par_df = pd.read_csv(par_file, delim_whitespace = True, names = ['parameter','first_value', 'second_value', 'std'])
    par_df.drop(par_df[par_df.parameter.str.contains('.offset')].index, inplace=True)
    par_df['value'] = par_df['second_value']
    par_df.drop(columns=['first_value','second_value'],inplace=True)
    
    return par_df


def create_new_parameter_distribution(par_df, min_max_range):
    '''Takes a par DataFrame and creates extra columns with new min/max ranges (but same std).
    Parameters
    ----------
    par_df: (PANDAS DataFrame) DataFrame containing parsed par file.
    min_max_ratio: (float) ratio such that new max/min = par_value*(1 +/- min_max_ratio*sign(par_value)).

    Returns
    -------
    new_par_df: (PANDAS DataFrame) DataFrame containing extra columns with new std and min/max ranges.
    '''
    new_par_df = par_df.copy()
    # Create new std, min and max
    # new max/min = value*(1 +/- min_max_range)
    new_par_df['new_min'] = new_par_df['value'].apply(lambda x: x - abs(x) * min_max_range)
    new_par_df['new_max'] = new_par_df['value'].apply(lambda x: x + abs(x) * min_max_range)
    new_par_df['new_std'] = new_par_df['std'].apply(lambda x: x)
    new_par_df['new_range_value'] = new_par_df.apply(lambda x: '{0:.8f}~{1:.8f}<{2:.8f},{3:.8f}>'.format(x.value,x.new_std,x.new_min,x.new_max),axis = 1)
    return new_par_df


def states_code_to_string(state_code):
    '''Takes a string with state code and returns the corresponding states string.
    Parameters
    ----------
    state_code: (string) code for gait state where each element is in ['1'..'5'].

    Returns
    -------
    TYPE: string containing the names of the states based off state_code.
    '''
    state_dict = {1 : 'Landing',
                 2 : 'Swing',
                 3 : 'Liftoff',
                 4 : 'LateStance',
                 5 : 'EarlyStance'}
    state_list = []
    for code_place in range(1,len(state_code)):
        if state_code[code_place] == '1':
            state_list.append(state_dict[code_place])
    return " ".join(state_list)


def apply_par_key_to_series(series, key_name, key_value, states, target, source = None):
    '''Sets a new key value on controller PANDAS Series, based on different filters.
    Parameters
    ----------
    series: (PANDAS Series) recursive series containing controller data.
    key_name: (string) name of the key which value is desired.
    key_value: (string or PANDAS Series) new value of the key from the controller block verifying desired conditions.
    states: (string) string containing desired controller states.
    target: (string) desired target muscle.
    source: (string) optional, desired source muscle. The default is None.

    Returns
    -------
    None.
    '''
    if source != None and srp.get_key(series, key_name, states, target, source) is None and srp.get_key(series, 'dof', states, target, None) != None:
        srp.set_key(series, key_name, key_value, states, target, None)
    else:
        srp.set_key(series, key_name, key_value, states, target, source)


def generate_key_df(par_df):
    '''Takes as input a raw DataFrame from a par file and pre-processes
    its column values and only extracts rows that are not threshold-related.
    Parameters
    ----------
    par_df: (PANDAS DataFrame) raw DataFrame from a par file.

    Returns
    -------
    key_df: (PANDAS DataFrame) pre-processed, filtered DataFrame.
    '''
    # extract key_df
    key_df = par_df[par_df.parameter.str.startswith('S')].copy()
    # extract states_code, muscles, key_name columns from parameter column.
    key_df['states_code'],key_df['muscles'], key_df['key_name'] = key_df['parameter'].str.split('.').str
    # split muscles column into target & source muscles.
    key_df['target'],key_df['source'] = key_df['muscles'].str.split('-', n = 2).str
    key_df = key_df.where(pd.notnull(key_df), None)
    # extract states name from states code.
    key_df['states'] = key_df['states_code'].apply(states_code_to_string)
    # drop unused columns
    key_df.drop(columns=['muscles','states_code'],inplace=True)
    
    return key_df


def generate_threshold_df(par_df):
    '''Extracts from a raw par file DataFrame the rows related to thresholds.
    Parameters
    ----------
    par_df: (PANDAS DataFrame) raw DataFrame from a par file.

    Returns
    -------
    threshold_df: (PANDAS DataFrame) DataFrame, containing only threshold-related rows.
    '''
    threshold_df = par_df[~par_df.parameter.str.startswith('S')].copy()

    return threshold_df
