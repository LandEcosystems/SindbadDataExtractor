#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 18:57:14 2021

@author: @dr-ko
"""
import json
import sys
import os
import copy
import utils.shared_utils as shut

sys.path.append(os.getcwd())


def add_unc_qc_vars(var, unc_qc_var_dic_in, _unc_qc, _unc_qc_var):
    """
    Generate the variable and information for QC and RANDUNC
     - assumes FLUXNET is following $VAR_QC or $VAR_RANDUNC as associated variables
     - assumes bounds of QC to be [0,1]
     - assumes bounds and units of RANDUNC to be the same as the parent $VAR
    """
    unc_qc_var_dic = unc_qc_var_dic_in.copy()

    if _unc_qc == _unc_qc_var:
        unc_qc_var_dic["sourceVariableName"] = f'{var}_{_unc_qc}'
    else:
        unc_qc_var_dic["sourceVariableName"] = _unc_qc_var

    unc_qc_var_dic["nameShort"] = f'{var}_{_unc_qc}'
    if len(unc_qc_var_dic_in["weightVar"]) > 0:
        unc_qc_var_dic["weightVar"] = unc_qc_var_dic_in["weightVar"]

    # unc_qc_var_dic["sourceVariableName"] = f'{var}_{_unc_qc}'

    if _unc_qc == 'QC':
        unc_qc_var_dic["bounds"] = [0, 1]
        unc_qc_var_dic["source2sindbadUnit"] = 1
        unc_qc_var_dic["variableUnit"] = "adimensional"
        unc_qc_var_dic["sourceVariableUnit"] = "adimensional"
        unc_qc_var_dic["isEnergy"] = False
        unc_qc_var_dic["isCarbon"] = False
        unc_qc_var_dic["isWater"] = False
        unc_qc_var_dic[
            "nameLong"] = "quality flag for " + unc_qc_var_dic_in["nameLong"]
    else:
        unc_qc_var_dic[
            "nameLong"] = "random uncertainty for " + unc_qc_var_dic_in[
                "nameLong"]

    # remove the unnecessary fields in the variable info for QC and RANDUNC
    # for unwant in ['QC', 'RANDUNC']:
    #     unc_qc_var_dic.pop(unwant)
    return unc_qc_var_dic


def get_cliff_dirname(fnversion, temp_reso):
    """
    Get the folder structure of cliff input and output based on a simple convention
    """
    if fnversion == 'FLUXNET2015':
        opath = 'FLUXNET2015'
    elif fnversion == 'LaThuile':
        opath = 'LaThuile'
    else:
        opath = fnversion
    opath = fnversion

    if temp_reso == 'daily':
        otime = 'DD'
    else:
        otime = 'HR'
    return f'{opath}.{otime}'


def setup_out_dir(exp_config):
    """
    setup and create the output fields and directories
    """
    out_main_dir = exp_config['output_dir_path']
    run_name = exp_config['run_name']
    out_sub_path = f"{run_name}.{get_cliff_dirname(exp_config['FLUXNET_version'], exp_config['temporal_resolution'])}"

    exp_config['output_dir_path'] = {}
    exp_config['output_dir_path']['main'] = os.path.join(out_main_dir, out_sub_path)
    opath = out_sub_path.split('.')[-2]
    otime = out_sub_path.split('.')[-1]
    exp_config['output_dir_path']['nc_file'] = os.path.join(
        exp_config['output_dir_path']['main'], 'data')
    os.makedirs(exp_config['output_dir_path']['nc_file'], exist_ok=True)
    # figures
    if exp_config['diagnostic_plots']:
        exp_config['output_dir_path']['figs'] = os.path.join(
            exp_config['output_dir_path']['main'], 'figs_diagno')
        os.makedirs(exp_config['output_dir_path']['figs'], exist_ok=True)
    # configuration output
    exp_config['output_dir_path']['expName'] = f'{run_name}.{opath}.{otime}'
    exp_config['output_dir_path']['info_config'] = os.path.join(
        exp_config['output_dir_path']['main'], 'config_info')
    exp_config['output_dir_path']['site_info'] = os.path.join(
        exp_config['output_dir_path']['main'], 'config_info','site_info')
    exp_config['output_dir_path']['log'] = os.path.join(exp_config['output_dir_path']['main'],
                                                'logs')
    os.makedirs(exp_config['output_dir_path']['log'], exist_ok=True)
    os.makedirs(exp_config['output_dir_path']['info_config'], exist_ok=True)
    os.makedirs(exp_config['output_dir_path']['site_info'], exist_ok=True)
    return exp_config


def merge_var_info(base_info, sel_var_info):
    """
    merge the base_info in each extractor with the corresponding changes from each variable within the extractor's json. 
    - The values in the base_info are default, and these are changed if they exist for each variable. 
    - This avoids repeating the blocks of same information for each variable in the json
    """
    out_info = {}
    for b_field in list(base_info.keys()):
        out_info[b_field] = base_info[b_field]
        if b_field in sel_var_info:
            out_info[b_field] = sel_var_info[b_field]
    return out_info


def get_selected_list(sel, full):
    """
    Get the list of selected variables or datasets
    - assumes fieldnames of the dictionary as the potential
    - returns a subset if the sel_datasets, or sel_variables in exp[json] provides a list
    - if a string is provided:
        - returns full list if it finds all
        - remove the variables if except is provided eg. 
        - e.g., 'all except var1 var2' returns the full list minus var1, var2
    """
    olist = full
    if isinstance(sel, list):
        for _sel in sel:
            if _sel.strip() == '':
                sel.remove(_sel)
        if "all" in sel:
            return olist
        else:
            return sel
    else:
        if sel == "all":
            return olist
        else:
            sel_list = [x.strip() for x in sel.strip().split(" ")]
            if "except" in sel_list:
                remove_list = [
                    x.strip()
                    for x in sel.split("except")[-1].strip().split(" ")
                ]
                return (list(set(olist) - set(remove_list)))
            elif "all" in sel_list:
                return olist
            else:
                return sel_list


def add_suffix_to_sel_variables(var_list, suffix):
    "adds the suffix to the variables that is provided as var_suffix in the json for each extractor"
    if len(suffix) > 0:
        return [_l + suffix for _l in var_list]
    else:
        return var_list


def get_variable_info(_config):
    _config_out = _config.copy()
    variables_all = {}
    for data_key in _config["sel_datasets"]:
        data_info = _config['dataset'][data_key]
        variables = {}
        var_sfx = data_info["var_suffix"].strip()
        if len(var_sfx) > 0:
            var_sfx = "_" + var_sfx
        var_info = data_info["variables"]
        sel_vars = get_selected_list(data_info["sel_variables"],
                                     list(data_info["variables"].keys()))
        data_info['sel_variables'] = sel_vars
        for var_base in data_info["sel_variables"]:
            var = var_base + var_sfx
            merged_var_info = merge_var_info(data_info["base_info"],
                                             var_info[var_base])
            variables[var] = merged_var_info.copy()

            variables[var]["nameLong"] = var_info[var_base]["nameLong"]
            # if _config["temporal_resolution"] == "daily":
            #     variables[var]["nameLong"] = var_info[var_base]["nameLong"]
            # else:
            #     variables[var]["nameLong"] = var_info[var_base]["nameLong"]

            unc_qc = ['RANDUNC', 'QC']
            for _unc_qc in unc_qc:
                if len(variables[var][_unc_qc].strip()) > 0:
                    unc_qc_var = variables[var][_unc_qc]
                    var_field_name = var_base + "_" + _unc_qc + var_sfx  #TA_QC
                    variables[var_field_name] = add_unc_qc_vars(
                        var_base, merged_var_info.copy(), _unc_qc, unc_qc_var)
                    variables[var][_unc_qc] = variables[var_field_name][
                        "sourceVariableName"]
            if _config["temporal_resolution"] == "daily":
                field_daystats = {
                    "DayMin": "daily minimum",
                    "DayMax": "daily maximum",
                    "DayTime": "daytime mean",
                    "DaySum": "daily sum",
                    "wDayMean": "weighted daily mean"
                }
                for fi_da, fi_prx in field_daystats.items():
                    if fi_da in var_info[var_base] and var_info[var_base][
                            fi_da]:
                        var_field_name = var_base + "_" + fi_da + var_sfx
                        variables[var_field_name] = merged_var_info.copy()
                        if fi_da == "wDayMean":
                            if len(variables[var_field_name]
                                   ["weightVar"]) == 0:
                                variables[var_field_name][
                                    "weightVar"] = var_base + "_QC"
                            variables[var_field_name][
                                "nameLong"] = f'{fi_prx} {var_info[var_base]["nameLong"]} [using {variables[var_field_name]["weightVar"]}]'
                        else:
                            variables[var_field_name][
                                "nameLong"] = f'{fi_prx} {var_info[var_base]["nameLong"]}'

                        for _unc_qc in unc_qc:
                            if len(variables[var][_unc_qc].strip()) > 0:
                                unc_qc_var = variables[var][_unc_qc]
                                var_field_name_qc = var_base + "_" + _unc_qc + "_" + fi_da + var_sfx
                                variables[var_field_name_qc] = add_unc_qc_vars(
                                    var_base, variables[var_field_name],
                                    _unc_qc, unc_qc_var)
                                variables[var_field_name][_unc_qc] = variables[
                                    var_field_name_qc]["sourceVariableName"]

            variables_all.update(variables)
        data_info["sel_variables"] = add_suffix_to_sel_variables(
            sel_vars, var_sfx)

        for var, var_det in variables.items():
            for unwant in "wDayMean DayMin DayMax DaySum DayTime".split():
                if unwant in var_det:
                    variables[var].pop(unwant)
        _config_out['dataset'][data_key].pop('base_info')
        _config_out["dataset"][data_key]['variables'] = variables
    return variables_all, _config_out


def get_extractor_configuration(conf_file):
    with open(conf_file, 'r') as config_file:
        extractor_rf = config_file.read()
    extractor_conf = json.loads(extractor_rf)
    return extractor_conf


def check_units_consistency(full_info):
    with open(
            os.path.join(full_info['output_dir_path']['info_config'],
                         'units_check.log'), 'w') as unf:
        for name in full_info['sel_datasets']:
            dataset = full_info['dataset'][name]

            unf.write(
                f'dataset: {name} | variables: {dataset["sel_variables"]}\n')

            print(f'dataset: {name} | variables: {dataset["sel_variables"]}')
            for var in list(dataset['variables'].keys()):
                # for var in dataset['sel_variables']:
                # print(var, dataset['variables'])
                info = dataset['variables'][var]
                src_unit = info["sourceVariableUnit"]
                tar_unit = info["variableUnit"]
                unit_scalar = info["source2sindbadUnit"]
                skip_check_unit = info["skipUnitcheck"]
                if skip_check_unit == False:
                    unf.write(
                        f"[UNIT-CHECK]{name}::{var}::src_unit: {src_unit}, tar_unit: {tar_unit}, unit_scalar: {unit_scalar}\n"
                    )
                    print(
                        f"[UNIT-CHECK]{name}::{var}::src_unit: {src_unit}, tar_unit: {tar_unit}, unit_scalar: {unit_scalar}"
                    )
                    if src_unit != tar_unit:
                        if float(unit_scalar) == 1:
                            sys.exit(
                                f'dataset: {name}: {var} unit is {src_unit} in source [sourceVariableUnit] and {tar_unit} in target [variableUnit], but the unit scalar [source2sindbadUnit] is {unit_scalar}. Fix inconsistency by setting the units and coversion correctly or using skipUnitcheck: true'
                            )
                else:
                    unf.write(
                        f"[UNIT-SKIP]:{name}::{var}::src_unit: {src_unit}, tar_unit: {tar_unit}, unit_scalar: {unit_scalar}\n"
                    )
                    print(
                        f"[UNIT-SKIP]:{name}::{var}::src_unit: {src_unit}, tar_unit: {tar_unit}, unit_scalar: {unit_scalar}"
                    )
            print('-' * 150)
            unf.write('-' * 150 + '\n')
    return


def check_duplicates_in_list(listOfElems):
    ''' Check if given list contains any duplicates '''
    if len(listOfElems) == len(set(listOfElems)):
        return False
    else:
        return True


def check_duplicate_variables(full_info):
    for name in full_info['sel_datasets']:
        dataset = full_info['dataset'][name]
        if check_duplicates_in_list(list(dataset['variables'].keys())):
            sys.exit(
                f"dataset: {name}: {dataset} has duplicates in the list of target variables (keys). Cannot continue: \n {sorted(list(dataset['variables'].keys()))}"
            )
    return


def check_duplicate_datasets(sel_dataset, check_info="datasets"):
    if check_duplicates_in_list(sel_dataset):
        sys.exit(
            f"Selected {check_info} in exp_config are duplicates. '{check_info}' list and keys of {check_info} need to be unique. Cannot continue: \n {sorted(sel_dataset)}"
        )
    else:
        print(f"No duplicates found in selected {check_info}: \n {sel_dataset}")
    return

def get_gapfill_settings(_config):
    '''
    check the gap filling settings
    '''
    _config['sel_gapfills'] = get_selected_list(
        _config['sel_gapfills'], list(_config['gap_fill'].keys()))
    if check_duplicates_in_list(_config['sel_gapfills']):
        sys.exit(
            f"Selected gap_fills in exp_config are duplicates. 'sel_gapfills' list and keys of gap_fill need to be unique. Cannot continue: \n {sorted(_config['sel_gapfills'])}"
        )

    sgf=_config['sel_gapfills']

    if len(sgf) == 0:
        _config['do_gap_fill'] = False
    if len(sgf) == 1 and (sgf[0].lower().strip() in ['none', 'no', 'nothing', '']):
        _config['do_gap_fill'] = False
    else:
        _config['do_gap_fill'] = True
    
    if _config['do_gap_fill']:
        for gfik in sgf:
            gfiv = _config['gap_fill'][gfik]
            check_fields = ['source', 'target']
            for _cf in check_fields:
                if gfiv[_cf] not in _config['sel_datasets']:
                    sys.exit(
                    f"Cannot do gap filling for {gfik}: {_cf} is set as {gfiv[_cf]}. But, {gfiv[_cf]} is not in selected datasets. \n {_config['sel_datasets']} \n Change sel_datasets or gap filling {_cf} or set sel_gapfills to any of [{' | '.join(['none', 'no', 'nothing'])}].")

    return _config

def get_inp_config(exp_config_path):
    with open(exp_config_path, 'r') as config_file:
        exp_config_rf = config_file.read()
    exp_config = json.loads(exp_config_rf)
    return exp_config


def get_exp_configuration(exp_config_path,
                          fn_version=None,
                          temporal_resolution=None,
                          output_dir_path=None,
                          run_name="fluxnetBGI"):
    # get input configuration
    print(fn_version, temporal_resolution, output_dir_path)
    exp_config_r = get_inp_config(exp_config_path)

    # replace version and time scale configs if passed
    if fn_version is not None:
        exp_config_r["FLUXNET_version"] = fn_version
    if temporal_resolution is not None:
        exp_config_r["temporal_resolution"] = temporal_resolution
    if output_dir_path is not None:
        exp_config_r["output_dir_path"] = output_dir_path
    exp_config_r["run_name"] = run_name
    print(exp_config_r["FLUXNET_version"], exp_config_r["temporal_resolution"], exp_config_r["output_dir_path"])
    exp_config = copy.deepcopy(exp_config_r)

    # setup directory
    exp_config = setup_out_dir(exp_config)
    exp_config['sel_datasets'] = get_selected_list(
        exp_config['sel_datasets'], list(exp_config['dataset'].keys()))

    ## check if the selected datasets have duplicates
    check_duplicate_datasets(exp_config['sel_datasets'])

    exp_config['resample_output'] = get_selected_list(
        exp_config['resample_output'], ['daily', 'monthly', 'annual'])

    ## check if the selected datasets have duplicates
    check_duplicate_datasets(exp_config['resample_output'], check_info="resample_output")

    # get extractor information for each of the selected datasets
    for extractor in exp_config['sel_datasets']:
        extractor_dict = get_extractor_configuration(os.path.join(exp_config['config_dir_path'], exp_config['dataset'][extractor]))
        exp_config['dataset'][extractor] = extractor_dict

    # get variable information for each of the selected extractor
    variables, exp_config = get_variable_info(exp_config)

    ## check if the selected and generated variables are duplicates
    check_duplicate_variables(exp_config)

    # check if the units are consistent with scalars
    check_units_consistency(exp_config)

    # get gap fill settings
    exp_config = get_gapfill_settings(exp_config)
    # save configuration
    shut.save_json(
        os.path.join(exp_config['output_dir_path']['info_config'],
                     "variable_info.json"), variables)
    shut.save_json(
        os.path.join(exp_config['output_dir_path']['info_config'], "exp_info.json"),
        exp_config)
    return exp_config


if __name__ == '__main__':
    import inspect, os
    print('---------------------------------------------------')
    print(
        'Utilities: ',
        os.path.dirname(
            os.path.abspath(inspect.getfile(inspect.currentframe()))))
    print('Config generator: ', inspect.getfile(inspect.currentframe()))
    print('---------------------------------------------------')
    print(__doc__)
