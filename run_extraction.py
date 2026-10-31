#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 18:57:14 2021

@author: sujan
"""
from multiprocessing import set_start_method, Pool
# set_start_method("spawn")
from multiprocessing import get_context


import sys
import os
import utils.shared_utils as shut
import utils.setup_config as set_conf
# At module level, so that pool workers started by spawn also have it.
from utils.workflow_extraction import compile_site_data
sys.path.append(os.getcwd())

def get_site_list(cubepath, version, use_site_list_from_file=False):
    import numpy as np
    site_list_path = os.path.join(exp_config['output_dir_path']['info_config'], f"site_list_{version}.csv")
    print(site_list_path)
    if os.path.exists(site_list_path) and use_site_list_from_file:
        site_list = [site.strip() for site in open(site_list_path).readlines()]
    elif exp_config.get("registry"):
        import utils.sites as sites
        reg = sites.load_registry(exp_config["registry"]["path"])
        site_list = sites.sites_with(exp_config["registry"].get("flag", "all"), reg)
        np.savetxt(site_list_path, site_list, delimiter=",", fmt="%s")
    else:
        import fluxcom.providers as prov
        ep = prov.EddyProvider(version=version, cubepath=cubepath)
        site_list = sorted(ep.sites)
        np.savetxt(site_list_path, site_list, delimiter=",", fmt="%s")
    return site_list

class DataCompiler:    
    def __init__(self, exp_settings=None):
        self.exp_settings = exp_settings
    def process(self, site):
        site_info = compile_site_data(site, self.exp_settings)
        settings = self.exp_settings
        if isinstance(site_info, dict):
            shut.save_json(os.path.join(settings['output_dir_path']['site_info'], f'{site}_info_{settings["FLUXNET_version"]}.json'), site_info)
        return 'done'


if __name__ == '__main__':
    #%% Retrieve parsed arguments
    if len(sys.argv) > 1:
        exp_config_path = sys.argv[1]
        if not os.path.isabs(exp_config_path):
            sys.exit(f'The provided path for exp config: {exp_config_path} is not absolute. Provide full absolute path to the configuration file.')
        if len(sys.argv) > 2:
            fn_versions =[sys.argv[2]]
        else:
            fn_versions=[]
        if len(sys.argv) > 3:
            tm_scales =[sys.argv[3]]
        else:
            tm_scales=[]
        if len(sys.argv) > 4:
            output_dir_path =sys.argv[4]
        else:
            output_dir_path=""
        if len(sys.argv) > 5:
            fluxcom_cube_path =sys.argv[5]
        else:
            fluxcom_cube_path="/Net/Groups/BGI/work_4/scratch/fluxcom/sitecube_proc/v0.07/site_cube"
        if len(sys.argv) > 6:
            run_name =sys.argv[6]
        else:
            run_name = "fluxnetBGI"
    else:
        sys.exit(f'run_extraction cannot extract data without the path of the main extraction config as the only argument. Provide absolute path to the configuration file.')
    # get basic config to check the sources and loop
    exp_config_r = set_conf.get_inp_config(exp_config_path)

    if len(fn_versions) == 0:
        fn_versions = exp_config_r["FLUXNET_version"]
        if isinstance(fn_versions, list):
            fn_versions=[_fn.strip() for _fn in fn_versions if len(_fn) != 0]
        else:
            fn_tmp = fn_versions.split(" ")
            fn_versions=[_fn.strip() for _fn in fn_tmp if len(_fn) != 0]

    if len(tm_scales) == 0:
        tm_scales = exp_config_r["temporal_resolution"]
        if isinstance(tm_scales, list):
            tm_scales=[_tm.strip() for _tm in tm_scales if len(_tm) != 0]
        else:
            tm_tmp = tm_scales.split(" ")
            tm_scales=[_tm.strip() for _tm in tm_tmp if len(_tm) != 0]
    if len(output_dir_path) == 0:
        output_dir_path = exp_config_r["output_dir_path"]

    # loop through sites and get full configurations
    for tm_scale in tm_scales:
        for fn_version in fn_versions:
            exp_config = set_conf.get_exp_configuration(exp_config_path, fn_version=fn_version, temporal_resolution=tm_scale, output_dir_path=output_dir_path, run_name=run_name)
            exp_config["fluxcom_cube_path"] = fluxcom_cube_path

            # #%% Check sites to process
            sites_in_version = get_site_list(exp_config["fluxcom_cube_path"], exp_config["FLUXNET_version"])
            if exp_config["sites_toRun"] != ['all']:
                list_sites= exp_config["sites_toRun"]
                sites = [site for site in sites_in_version if site in list_sites]
            else:
                sites = sites_in_version
                list_sites = sites_in_version
            # if f'sites_{fn_version}' not in locals():
            #     sites_in_version = get_site_list(exp_config["fluxcom_cube_path"], exp_config["FLUXNET_version"])
            #     if exp_config["sites_toRun"] != ['all']:
            #         list_sites= exp_config["sites_toRun"]
            #         sites = [site for site in sites_in_version if site in list_sites]
            #     else:
            #         sites = sites_in_version
            #         list_sites = sites_in_version
            # else:
            #     sites = eval(f'sites_{fn_version}')
            
            if len(sites) == 0:
                print(f'No sites out of {list_sites} to process for {fn_version} and temporal resolution {tm_scale} which has {sites_in_version}. Skipping to the next configuration.')
                continue
            if len(sites) < len(list_sites):
                print(f'Only {len(sites)} out of {len(list_sites)} sites to process for {fn_version} and temporal resolution {tm_scale} which has {sites_in_version}. Skipping the sites not in the list: {[site for site in list_sites if site not in sites]}.')

            print(f'{exp_config["FLUXNET_version"]}: {sites}')

            exp_settings = {'exp_settings' : exp_config}
            data_process_run = DataCompiler(**exp_settings)
            if  exp_config["njobs_parallel"] > 1:
                p = Pool(exp_settings['exp_settings']['njobs_parallel'], maxtasksperchild=1)
                site_info_dict = p.map(data_process_run.process, sites)
                p.close()
                p.join()
            else:
                site_info_dict ={}
                for site in sites:
                    print (f'Processing site: {site} for version: {fn_version} and temporal resolution: {tm_scale}')
                    site_info_dict[site] = compile_site_data(site, **exp_settings)
            print('The extraction is complete. Back to run_extract.py.')

            print('----------------------------\n'
                'End of model data processing\n'
                '----------------------------')
