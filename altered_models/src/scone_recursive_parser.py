"""
This file contains multiple functions to handle SCONE controller files, namely
parsing SCONE files to PANDAS Series, filtering, getter and setter methods,
exporting to SCONE files, etc.
"""
import os, re, shutil, pprint
import pandas as pd
import numpy as np
from collections import OrderedDict
from copy import deepcopy


def series_deep_copy(series):
    '''Performs deep copying of an input PANDAS Series.
    Parameters
    ----------
    series: (PANDAS Series) input series.

    Returns
    -------
    series_out: (PANDAS Series) deepcopy of series.
    '''
    series_out = pd.Series(dtype=object)
    for idx in series.index:
        if isinstance(series[idx], pd.Series):
            series_out[idx] = series_deep_copy(series[idx])
        else:
            series_out[idx] = deepcopy(series[idx])
    return series_out


def bracketed_split(string, delimiter, strip_brackets=False):
    '''Splits a string by the delimiter unless it is inside brackets.
    e.g.
        list(bracketed_split('abc,(def,ghi),jkl', delimiter=',')) == ['abc', '(def,ghi)', 'jkl']
    Parameters
    ----------
    string: (str) input string.
    delimiter: (character or string) delimiter to split by.
    strip_brackets: (bool) optional, whether to leave brackets after splitting or not. The default is False.

    Yields
    ------
    current_string: generator object such that list(current_string) is a list containing split input string.
    '''
    string = string.strip()
    openers = '[{('
    closers = ']})'
    opener_to_closer = dict(zip(openers, closers))
    opening_bracket = dict()
    current_string = ''
    depth = 0
    for c in string:
        if c in openers:
            depth += 1
            opening_bracket[depth] = c
            if strip_brackets and depth == 1:
                continue
        elif c in closers:
            assert depth > 0, f"You exited more brackets that we have entered in string {string}"
            assert c == opener_to_closer[opening_bracket[
                depth]], f"Closing bracket {c} did not match opening bracket {opening_bracket[depth]} in string {string}"
            depth -= 1
            if strip_brackets and depth == 0:
                continue
        if depth == 0 and c == delimiter:
            yield current_string
            current_string = ''
        else:
            current_string += c
    assert depth == 0, f'You did not close all brackets in string {string}'
    yield current_string


def dict_recursive_split(string):
    '''Performs a recursive split of a SCONE controller file contained in a string,
    where all lines are separated with '\n'. Returns a recursive dictionnary,
    containing the fields of the controller file.
    Parameters
    ----------
    string: (string) SCONE Controller file in a string.

    Returns
    -------
    dict_out: (dict) recursive dictionnary containing the fields of the controller file.
    '''
    dict_out = OrderedDict()
    repeated_container_index = 0
    repeated_comment_index = 0
    # split controller string into lines leaving blocks between brackets whole.
    first_level_list = list(bracketed_split(string, delimiter='\n'))
    # Remove trailing white characters from each string.
    first_level_list = list(map(str.strip, first_level_list))
    # Loop over each element of the split string.
    for first_level_element in first_level_list:
        # '{' indicates a container.
        if '{' in first_level_element:
            # split controller block lines again, leaving control blocks whole
            second_level_list = list(bracketed_split(first_level_element, delimiter=' ', strip_brackets=True))
            # remove trailing white space.
            second_level_list = list(map(str.strip, second_level_list))
            # assert controller block integrity.
            assert len(second_level_list) == 2
            # call the same method recursively.
            recursive_dict = dict_recursive_split(second_level_list[1])
            # if at the bottom of rescurisve depth, i.e. no further elements are detected, append the original list.
            if recursive_dict == {}:
                dict_out[second_level_list[0] + '#' + str(repeated_container_index)] = second_level_list[1]
            # otherwise, append the obtained dict.
            else:
                dict_out[second_level_list[0] + '#' + str(repeated_container_index)] = recursive_dict
            # increment index for multiple containers of the same name.
            repeated_container_index += 1

        # Append definitions of parameters to the dict. Parameter values can also be between brackets.
        elif '=' in first_level_element:
            second_level_list = list(bracketed_split(first_level_element, delimiter='='))
            second_level_list = list(map(str.strip, second_level_list))
            assert len(second_level_list) == 2
            dict_out[second_level_list[0]] = second_level_list[1]
        # Append commented lines from the controller file to the dict.
        elif '#' in first_level_element:
            dict_out['Comment#' + str(repeated_comment_index)] = first_level_element
            repeated_comment_index += 1
            # recursively return output dict.
    return dict_out


def series_recursive_split(string):
    '''Performs a recursive split of a SCONE controller file contained in a string,
    where all lines are separated with '\n'. Returns a recursive PANDAS Series,
    containing the fields of the controller file.
    Parameters
    ----------
    string: (string) SCONE Controller file in a string.

    Returns
    -------
    series_out: (PANDAS Series) recursive series containing the fields of the controller file.
    '''
    series_out = pd.Series(dtype=object)
    repeated_container_index = 0
    repeated_comment_index = 0
    # split controller string into lines leaving blocks between brackets whole.
    first_level_list = list(bracketed_split(string, delimiter='\n'))
    # Remove trailing white characters from each string.
    first_level_list = list(map(str.strip, first_level_list))
    # Loop over each element of the split string.
    for first_level_element in first_level_list:
        # '{' indicates a container.
        if '{' in first_level_element:
            # split controller block lines again, leaving control blocks whole
            second_level_list = list(bracketed_split(first_level_element, delimiter=' ', strip_brackets=True))
            # remove trailing white space.
            second_level_list = list(map(str.strip, second_level_list))
            # assert controller block integrity.
            assert len(second_level_list) == 2
            # call the same method recursively.
            recursive_series = series_recursive_split(second_level_list[1])
            # if at the bottom of rescurisve depth, i.e. no further elements are detected, append the original list.
            if recursive_series.empty:
                series_out[second_level_list[0] + '#' + str(repeated_container_index)] = second_level_list[1]
            # otherwise, append the obtained Series.
            else:
                series_out[second_level_list[0] + '#' + str(repeated_container_index)] = recursive_series
            # increment index for multiple containers of the same name.
            repeated_container_index += 1

        # Append definitions of parameters to the dict. Parameter values can also be between brackets.
        elif '=' in first_level_element:
            second_level_list = list(bracketed_split(first_level_element, delimiter='='))
            second_level_list = list(map(str.strip, second_level_list))
            assert len(second_level_list) == 2
            series_out[second_level_list[0]] = second_level_list[1]
        # Append commented lines from the controller file to the dict.
        elif '#' in first_level_element:
            series_out['Comment#' + str(repeated_comment_index)] = first_level_element
            repeated_comment_index += 1
    # recursively return output Series.
    return series_out


def auto_indent_output_str(output_str):
    '''Takes as input a non indented output string and indents it automatically
    so that it is better readeable when exported into a file.
    Example:
        output_str = 'abc {\n
                      def\n
                      }'
        indented_str= 'abc {\n
                       \tdef\n
                       }'
    Parameters
    ----------
    output_str: (string) string containing multiple lines and brackets.

    Returns
    -------
    indented_str: (string) string where lines are properly indented according to bracket levels.
    '''
    indented_str = ""
    tabdepth = 0
    # split string by lines and loop through them.
    for line in output_str.split("\n"):
        # remove trailing whitespaces or tabulations before auto-indent.
        line = line.lstrip()
        depth_changed = False
        indentedline = ''.join(['\t' * tabdepth, line])
        # adapt line depth according to brackets
        for c in line:
            if c == '{':
                tabdepth += 1

            elif c == '}':
                tabdepth -= 1
                depth_changed = True
        # indent line in case of depth change.
        if depth_changed:
            indentedline = ''.join(['\t' * tabdepth, line])
        # append indented line to output string.
        indented_str += indentedline + '\n'
    return indented_str


def export_series_as_file_str(series):
    '''Exports a recursive PANDAS Series representing a SCONE controller
    file to a string.
    Parameters
    ----------
    series: (PANDAS Series) series containing information about controller file.

    Returns
    -------
    str_out: (string) string equivalent to series controller file, where lines are separated by '\n'.
    '''
    str_out = ''
    series.dropna(inplace=True)
    # loop through elements of series.
    for element_idx in series.index:
        # element is also a PANDAS Series --> Call this method recursively on this element.
        if isinstance(series[element_idx], pd.Series):
            str_out += element_idx.split('#', 1)[0] + ' {\n' + export_series_as_file_str(series[element_idx]) + '}\n'
        # Comment --> line starts with '#'
        elif 'Comment' in element_idx:
            str_out += series[element_idx] + '\n'
        # element is string --> parameter definition.
        elif isinstance(series[element_idx], str):
            str_out += element_idx + ' = ' + series[element_idx] + '\n'
    # return recursively-appended string.
    return str_out


def series_to_scone_file(series, file_path):
    '''Exports a controller Series to a SCONE file.
    Parameters
    ----------
    series: (PANDAS Series) series containing information about controller file.
    file_path: (string) export path of SCONE file.

    Returns
    -------
    None.
    '''
    # auto indent series.
    output_str = auto_indent_output_str(export_series_as_file_str(series))
    # save to file.
    with open(file_path, "w") as text_file:
        text_file.write("{0}".format(output_str))


def controller_file_to_string(controller_file_path):
    '''Extracts controller file and returns a string containing the
    content of that file.
    Parameters
    ----------
    controller_file_path: (string) path to controller file.

    Returns
    -------
    string: string containing controller file content, with lines separated by '\n'.
    '''
    with open(controller_file_path, "r") as controller_file:
        controller_data = controller_file.readlines()
    return ''.join(controller_data)


def are_states_equal(state_str_1, state_str_2):
    '''Returns whether two state strings from a controller file are equivalent.
    Example: are_states_equal('LateStance LiftOff', 'LiftOff LateStance') == True
    Parameters
    ----------
    state_str_1: (string) string containing controller states.
    state_str_2: (string) string containing controller states.

    Returns
    -------
    bool: True if two states are equivalent.
    '''
    states_list_1 = list(filter(None, re.split('\W+', state_str_1)))
    states_list_1.sort()
    states_list_2 = list(filter(None, re.split('\W+', state_str_2)))
    states_list_2.sort()
    return states_list_1 == states_list_2


def find_by_states(series, desired_states):
    '''Returns a Series containing all the controller blocks acting during
    the desired states.
    Parameters
    ----------
    series: (PANDAS Series) recursive series containing controller data.
    desired_states: (string) string containing desired controller states.

    Returns
    -------
    found_series: (PANDAS Series) series containing the applicable controller blocks.
    '''
    found_series = pd.Series(dtype=object)
    if 'states' in series.index:
        if are_states_equal(desired_states, series['states']):
            found_series = found_series.append(pd.Series([series]), ignore_index=True)
    else:
        for k in series:
            if isinstance(k, pd.Series):
                found_series = found_series.append(find_by_states(k, desired_states), ignore_index=True)
    return found_series


def find_by_target_source(series, target, source=None):
    '''Filters input series and returns a Series containing all the
    controller blocks based on desired target and source muscle.
    If source is None, source is considered the same as target muscle.
    Parameters
    ----------
    series: (PANDAS Series) recursive series containing controller data.
    target: (string) desired target muscle.
    source: (string) optional, desired source muscle. The default is None.

    Returns
    -------
    found_series: (PANDAS Series) series containing the applicable controller blocks.
    '''
    found_series = pd.Series(dtype=object)
    if 'target' in series.index and target in series['target']:
        if source is None:
            if 'source' not in series.index:
                found_series = found_series.append(pd.Series([series]), ignore_index=True)
        elif 'source' in series.index and source in series['source']:
            found_series = found_series.append(pd.Series([series]), ignore_index=True)
    else:
        for k in series:
            if isinstance(k, pd.Series):
                found_series = found_series.append(find_by_target_source(k, target, source), ignore_index=True)
    return found_series


def find_by_unique_key(series, key_name):
    '''Recursively filters input series and returns a Series containing
    controller block which have key_name as part of them. It is assumend that
    the key will be present at only one occurence. Otherwise, the method includes
    the first occurence
    Example: find_by_unique_key({'a' : 0, 'b' : {'a' : 1, 'c' : 9}},'a')
                == { {'a' : 0, 'b' : {'a' : 1, 'c' : 9}} }
    Parameters
    ----------
    series: (PANDAS Series) recursive series containing controller data.
    key_name: (string) desired key to filter with.

    Returns
    -------
    found_series: (PANDAS Series) series containing the applicable controller blocks.
    '''
    found_series = pd.Series(dtype=object)
    if key_name in series.index:
        found_series = found_series.append(pd.Series([series]), ignore_index=True)
    else:
        for k in series:
            if isinstance(k, pd.Series):
                found_series = found_series.append(find_by_unique_key(k, key_name), ignore_index=True)
    return found_series


def find_by_states_target_source(series, desired_states, target, source=None):
    '''Filters input series and returns a Series containing all the
    controller blocks based on desired gait states, desired target and source muscle.
    If source is None, source is considered the same as target muscle.
    Parameters
    ----------
    series: (PANDAS Series) recursive series containing controller data.
    desired_states: (string) string containing desired controller states.
    target: (string) desired target muscle.
    source: (string) optional, desired source muscle. The default is None.

    Returns
    -------
    found_target_source: (PANDAS Series) series containing the applicable controller blocks.
    '''
    found_states = find_by_states(series, desired_states)
    found_target_source = find_by_target_source(found_states, target, source)
    return found_target_source


def find_by_key(series, key_name):
    '''Filters input series and returns a Series containing controller block
    which have key_name as part of them. Multiple occurences of the same key can be handled
    Example: find_by_unique_key({'a' : 0, 'b' : -1,'a' : 1,'a')
                == { {'a' : 0, 'a' : 1} }
    Parameters
    ----------
    series: (PANDAS Series) recursive series containing controller data.
    key_name: (string) desired key to filter with.

    Returns
    -------
    found_series: (PANDAS Series) series containing the applicable controller blocks.
    '''
    found_series = pd.Series(dtype=object)
    for element in series:
        if key_name in element.index:
            found_series = found_series.append(pd.Series([element]), ignore_index=True)
    return found_series


def get_key(series, key_name, states, target, source=None):
    '''Filters input series based on desired gait states, desired target
    and source muscle. If source is None, source is considered the same as target
    muscle. The method returns the value of the desired key.
    Parameters
    ----------
    series: (PANDAS Series) recursive series containing controller data.
    key_name: (string) name of the key which value is desired.
    states: (string) string containing desired controller states.
    target: (string) desired target muscle.
    source: (string) optional, desired source muscle. The default is None.

    Returns
    -------
    string or PANDAS Series: value of the key from the controller block verifying desired conditions.
    '''
    found_by_states_target_source = find_by_states_target_source(series, states, target, source)
    found_by_key = find_by_key(found_by_states_target_source, key_name)
    assert found_by_key.size <= 1, "Error: multiple elements verifying input conditions.\n{0}".format(
        str(found_by_key.to_json(indent=4)))

    if found_by_key.size == 0:
        return None
    else:
        return found_by_key[0][key_name]


def set_key(series, key_name, key_value, states, target, source=None):
    '''Filters input series based on desired gait states, desired target
    and source muscle. If source is None, source is considered the same as target
    muscle. The method then sets the value of the desired key.
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
    found_by_states_target_source = find_by_states_target_source(series, states, target, source)
    found_by_key = find_by_key(found_by_states_target_source, key_name)

    assert found_by_key.size >= 1, "Error: multiple elements verifying input conditions.\n{0}".format(
        str(found_by_key.to_json(indent=4)))
    assert found_by_key.size != 0, "Error: Key not found for input conditions.\n{0}".format(
        ' '.join([key_name, states, target, str(source)]))
    found_by_key[0][key_name] = str(key_value)


def set_key_percent(series, key_name, key_value, states, target, source=None):
    '''Filters input series based on desired gait states, desired target
    and source muscle. If source is None, source is considered the same as target
    muscle. The method then sets the percent value of the desired key.
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
    found_by_states_target_source = find_by_states_target_source(series, states, target, source)
    found_by_key = find_by_key(found_by_states_target_source, key_name)

    assert found_by_key.size >= 1, "Error: multiple elements verifying input conditions.\n{0}".format(
        str(found_by_key.to_json(indent=4)))
    assert found_by_key.size != 0, "Error: Key not found for input conditions.\n{0}".format(
        ' '.join([key_name, states, target, str(source)]))
    found_by_key[0][key_name] = str(key_value/100 * float(found_by_key[0][key_name].split('~')[0]))


def get_unique_key(series, key_name):
    '''Recursively filters input series and returns the value of the
    filtered key. It is assumend that  the key will be present at only one
    occurence. Otherwise, the method fails.
    Parameters
    ----------
    series: (PANDAS Series) recursive series containing controller data.
    key_name: (string) desired key to filter with.

    Returns
    -------
    string or PANDAS Series: value of the key from the controller block verifying desired conditions.
    '''
    found_by_key = find_by_unique_key(series, key_name)

    assert found_by_key.size <= 1, "Error: multiple elements verifying input conditions.\n{0}".format(
        str(found_by_key.to_json(indent=4)))

    if found_by_key.size == 0:
        return None
    else:
        return found_by_key[0][key_name]


def set_unique_key(series, key_name, key_value):
    '''Recursively filters input series and sets the new value of the
    filtered key. It is assumend that  the key will be present at only one
    occurence. Otherwise, the method fails.
    Parameters
    ----------
    series: (PANDAS Series) recursive series containing controller data.
    key_name: (string) desired key to filter with.
    key_value: (string or PANDAS Series) new value of the key from the controller block verifying desired conditions.

    Returns
    -------
    None.
    '''
    found_by_key = find_by_unique_key(series, key_name)

    assert found_by_key.size <= 1, "Error: multiple elements verifying input conditions.\n{0}".format(
        str(found_by_key.to_json(indent=4)))
    assert found_by_key.size != 0, "Error: Key not found for input conditions.\n"

    found_by_key[0][key_name] = str(key_value)


def generate_controllers_with_key_combinations(series, key_dict_list, output_folder_path, scone_base_files_path,
                                                 main_file_name):
    '''Takes as input a list of dictionnaries, each one defining multiple altered values of a reflex parameter
    and combines them together to generate the folder structure used to run SCONE simulations for each
    combination of values.
    key_dict_list is a list of dictionaries, each one defining an altered parameter for two targets
    with pourcent key_values ex:
    key_dict_0 = {'key_name': 'KV', 'states': 'LateStance EarlyStance Liftoff', 'target': ['soleus', 'gastroc'],
             'key_values': range(0,100,10)}
    key_dict_1 = {'key_name': 'KV', 'states': 'Swing Landing', 'target': ['soleus', 'gastroc'],
             'key_values': range(0,100,10)}
    key_dict_list = [key_dict_0, key_dict_1]
    Parameters
    ----------
    series: (PANDAS Series) recursive series containing controller data.
    key_dict_list: (list) list of dictionaries, each one defining multiple altered values of a reflex parameter.
    output_folder_path: (string) export folder path.
    scone_base_files_path: (string) path to SCONE base files, necessary to run SCONE simulations.
    main_file_name: (string) name of the main SCONE file.

    Returns
    -------
    None.
    '''
    prefix = output_folder_path.split('/')[-1]
    value_ranges = [key_dict['key_values'] for key_dict in key_dict_list]
    n_dict = len(key_dict_list)
    mesh = np.array(np.meshgrid(*value_ranges))
    value_combinations = mesh.T.reshape(-1, n_dict)

    if not os.path.exists(output_folder_path):
        os.mkdir(output_folder_path)

    for value_combination in value_combinations:
        combination_series = series_deep_copy(series)
        combination_dict_list = key_dict_list.copy()
        combination_name = prefix
        for dict_id in range(n_dict):
            key_dict = combination_dict_list[dict_id]
            key_value = value_combination[dict_id]
            key_dict['key_values'] = key_value
            for target in key_dict['target']:
                set_key_percent(combination_series, key_dict['key_name'], key_value, key_dict['states'], target,
                        key_dict['source'])
            combination_name = '_'.join([combination_name, str(key_value)])

        combination_folder_path = output_folder_path + '/' + combination_name
        if not os.path.exists(combination_folder_path):
            os.mkdir(combination_folder_path)

        # export new scone controller file
        series_to_scone_file(combination_series, '{}/gait_controller.scone'.format(combination_folder_path))
        # copy scone base files
        shutil.copytree(scone_base_files_path, combination_folder_path + '/base_files')
        # generate txt file for changed parameters

        with open('{}/changed_parameters.txt'.format(combination_folder_path), 'w') as fp:
            fp.write(pprint.pformat(combination_dict_list))

        combination_main_file_path = combination_folder_path + '/base_files/' + main_file_name
        set_unique_key_on_file(combination_main_file_path, 'signature_prefix', '"{}"'.format(combination_name))


def generate_scone_main_files0(series, key_dict, output_folder_path, scone_base_files_path,
                                               main_file_name):
    '''Takes as input a dictionnary defining multiple altered values of a biomechanical parameter.
    and generate the folder structure used to run SCONE simulations for each combination of values.
    key_dict is a dictionary defining an altered parameter for two targets with pourcent key_values ex:
    key_dict = {'key_name': 'max_isometric_force', 'target': ['soleus', 'gastroc'], 'key_values': range(0,100,10)}
    Parameters
    ----------
    series: (PANDAS Series) recursive series containing controller data.
    key_dict: (dict) dictionary defining multiple altered values of a biomechanical parameter.
    output_folder_path: (string) export folder path.
    scone_base_files_path: (string) path to SCONE base files, necessary to run SCONE simulations.
    main_file_name: (string) name of the main SCONE file.

    Returns
    -------
    None.
    '''
    prefix = output_folder_path.split('/')[-1]
    value_combinations = key_dict['key_values']

    if not os.path.exists(output_folder_path):
        os.mkdir(output_folder_path)

    for value_combination in value_combinations:
        combination_series = series_deep_copy(series)
        combination_dict = key_dict.copy()
        combination_name = prefix
        key_dict = combination_dict
        key_value = value_combination
        key_dict['key_values'] = key_value
        combination_name = '_'.join([combination_name, str(key_value)])

        combination_folder_path = output_folder_path + '/' + combination_name
        if not os.path.exists(combination_folder_path):
            os.mkdir(combination_folder_path)

        # export new scone controller file
        series_to_scone_file(combination_series, '{}/gait_controller.scone'.format(combination_folder_path))
        # copy scone base files
        shutil.copytree(scone_base_files_path, combination_folder_path + '/base_files')

        combination_main_file_path = combination_folder_path + '/base_files/' + main_file_name
        set_unique_key_on_file(combination_main_file_path, 'signature_prefix', '"{}"'.format(combination_name))

        scone_main_file_path = combination_folder_path + '/base_files/' + main_file_name
        set_scone_main_file_gain(scone_main_file_path, prefix, key_dict['target'], key_value)


def generate_scone_main_files(series, key_dict_list, output_folder_path, scone_base_files_path,
                                               main_file_name):
    '''Takes as input a dictionnary defining multiple altered values of a biomechanical parameter.
    and generate the folder structure used to run SCONE simulations for each combination of values.
    key_dict is a dictionary defining an altered parameter for two targets with pourcent key_values ex:
    key_dict = {'key_name': 'max_isometric_force', 'target': ['soleus', 'gastroc'], 'key_values': range(0,100,10)}
    Parameters
    ----------
    series: (PANDAS Series) recursive series containing controller data.
    key_dict: (dict) dictionary defining multiple altered values of a biomechanical parameter.
    output_folder_path: (string) export folder path.
    scone_base_files_path: (string) path to SCONE base files, necessary to run SCONE simulations.
    main_file_name: (string) name of the main SCONE file.

    Returns
    -------
    None.
    '''
    prefix = output_folder_path.split('/')[-1]

    value_ranges = [key_dict['key_values'] for key_dict in key_dict_list]
    n_dict = len(key_dict_list)
    mesh = np.array(np.meshgrid(*value_ranges))
    value_combinations = mesh.T.reshape(-1, n_dict)

    if not os.path.exists(output_folder_path):
        os.mkdir(output_folder_path)

    for value_combination in value_combinations:
        combination_series = series_deep_copy(series)
        combination_dict_list = key_dict_list.copy()
        combination_name = prefix
        param_names = []
        targets = []
        key_values = []
        for dict_id in range(n_dict):
            key_dict = combination_dict_list[dict_id]
            key_value = value_combination[dict_id]
            key_dict['key_values'] = key_value
            combination_name = '_'.join([combination_name, str(key_value)])
            param_names.append(key_dict['key_name'])
            targets.append(key_dict['target'])
            key_values.append(key_value)

        combination_folder_path = output_folder_path + '/' + combination_name
        if not os.path.exists(combination_folder_path):
            os.mkdir(combination_folder_path)

        # export new scone controller file
        series_to_scone_file(combination_series, '{}/gait_controller.scone'.format(combination_folder_path))
        # copy scone base files
        shutil.copytree(scone_base_files_path, combination_folder_path + '/base_files')

        combination_main_file_path = combination_folder_path + '/base_files/' + main_file_name
        set_unique_key_on_file(combination_main_file_path, 'signature_prefix', '"{}"'.format(combination_name))

        scone_main_file_path = combination_folder_path + '/base_files/' + main_file_name
        set_scone_main_file_gain(scone_main_file_path, param_names, targets, key_values)


def set_scone_main_file_gain(scone_main_file, param_names, targets, key_values):
    """Sets new key value of a biomechanical parameter in scone main file.
    Parameters
    --------
    """
    file = open(scone_main_file, 'r')
    lines = file.readlines()
    new_lines = []
    for l in range(len(lines)):
        new_lines.append(lines[l])
        if len(lines[l].split()) > 1 and lines[l].split()[-1] == '"*_tx;*_ty;*_u"':
            new_lines.append('\n')
            new_lines.append('\t\t\tProperties {\n')
            for p in range(len(param_names)):
                for i in range(len(targets[p])):
                    new_lines.append('\t\t\t\t'+targets[p][i]+'_l {' + param_names[p] + '.factor = ' + str(key_values[p]/100) + '}\n')
                    new_lines.append('\t\t\t\t'+targets[p][i]+'_r {' + param_names[p] + '.factor = ' + str(key_values[p]/100) + '}\n')
            new_lines.append('\t\t\t\t}\n')
    file.close()

    # write new file
    file = open(scone_main_file, 'w')
    for l in new_lines:
        file.write(l)
    file.close()

def set_scone_main_file_gain0(scone_main_file, param_name, targets, key_value):
    """Sets new key value of a biomechanical parameter in scone main file.
    Parameters
    --------
    """
    file = open(scone_main_file, 'r')
    lines = file.readlines()
    new_lines = []
    for l in range(len(lines)):
        new_lines.append(lines[l])
        if len(lines[l].split()) > 1 and lines[l].split()[-1] == '"*_tx;*_ty;*_u"':
            new_lines.append('\n')
            new_lines.append('\t\t\tProperties {\n')
            new_lines.append('\t\t\t\t'+targets[0]+'_l {' + param_name + '.factor = ' + str(key_value/100) + '}\n')
            new_lines.append('\t\t\t\t'+targets[0]+'_r {' + param_name + '.factor = ' + str(key_value/100) + '}\n')
            new_lines.append('\t\t\t\t'+targets[1]+'_l {' + param_name + '.factor = ' + str(key_value/100) + '}\n')
            new_lines.append('\t\t\t\t'+targets[1]+'_r {' + param_name + '.factor = ' + str(key_value/100) + '}\n')
            new_lines.append('\t\t\t\t}\n')
    file.close()

    # write new file
    file = open(scone_main_file, 'w')
    for l in new_lines:
        file.write(l)
    file.close()


def scone_file_to_list(scone_file_path):
    '''Reads a SCONE file and returns its lines as a list.
    Parameters
    ----------
    scone_file_path: (string) path to SCONE file.

    Returns
    -------
    scone_data: (list) list containing SCONE file lines.
    '''
    with open(scone_file_path, "r") as scone_file:
        scone_data = scone_file.readlines()
    return scone_data


def scone_list_to_file(scone_list, scone_file_path):
    '''Exports a list of strings to a file.
    Parameters
    ----------
    scone_list: (list) list containing SCONE file lines.
    scone_file_path: (string) path of output SCONE file.

    Returns
    -------
    None.
    '''
    with open(scone_file_path, "w") as scone_file:
        scone_file.writelines(scone_list)


def change_key_in_scone_list(scone_list, key_name, key_value):
    '''Sets a new value of key from a file list and returns the modified list.
    In case of multiple occurences, only the first one will be changed.
    Parameters
    ----------
    scone_list: (list) list containing SCONE file lines.
    key_name: (string) name of the desired key to change.
    key_value: (string) new value to be set.

    Returns
    -------
    new_scone_list: (list) modified scone_list.
    '''
    new_scone_list = scone_list.copy()
    for line_id in range(len(new_scone_list)):
        key_found = re.search('^\s+{}'.format(key_name), new_scone_list[line_id])
        if key_found != None:
            new_scone_list[line_id] = key_name + ' = ' + key_value
            break
    return new_scone_list


def set_unique_key_on_file(file_path, key_name, key_value):
    '''Modifies a key value directly on a scone file. In case of
    multiple occurences, only the first one will be changed.
    Parameters
    ----------
    file_path: (string) path of the file to be modified.
    key_name: (string) name of the desired key to change.
    key_value: (string) new value to be set.

    Returns
    -------
    None.
    '''
    scone_list = scone_file_to_list(file_path)
    new_scone_list = change_key_in_scone_list(scone_list, key_name, key_value)
    scone_list_to_file(new_scone_list, file_path)
