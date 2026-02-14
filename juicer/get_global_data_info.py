import os
import xarray as xr
import json
import numpy as np
fnsites= np.sort(['AU-Tum', 'CA-TP1','CN-Du1' ,'ES-LMa', 'PT-Mi1', 'DE-Hai','IT-SRo', "ZA-Kru", "BR-Sa1", "FI-Hyy", 'RU-Zot'])
fnvars = ['SW_IN','TA','P','VPD','NETRAD','LW_IN','TA_DayTime','VPD_DayTime','TA_DayMax','TA_DayMin','WS','RH','PA']

fnsites= np.sort(['AU-Tum', "BR-Sa1","BR-Sa3", "US-Wkg", 'CA-TP1','CN-Cha', 'DE-Hai' ,'ES-LJu','IT-MBo', "IT-BCi", "ZA-Kru", "FI-Hyy", 'RU-Fyo', 'US-Ton', 'US-Ha1'])
fnvars = ['SW_IN','TA','P','VPD','NETRAD','LW_IN','WS','RH','PA']
run_name = "fluxnetBGI2021"
out_dir_name = f"{run_name}_siteCube-v0.03_iter-1"
data_sets={
    "CRUNCEP": {
        "path":"/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d50_daily/CRUNCEP/v8/Data/variable",
        "variables":["tair","qair","rain","swdown","lwdown","tmax","tmin","wind","press"],
        "filepatt": "variable.CRUNCEP.v8.720.360.year.nc",
        "syear": "1901",
        "eyear": "2020",
        "t_res": "daily"
    },
    "CRUJRA": {
        "path":"/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d50_daily/CRUJRA/v2_2/Data/variable",
        "variables":"lwdown press qair rain swdown tair tair_day tmax tmin vpd_day wind".split(),
        "filepatt": "variable.CRUJRA_v2_2.720.360.year.nc",
        "t_res": "daily"
    },
    "WFDEI": {
        "path":"/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d50_daily/WFDEI/v1/Data/variable",
        "variables":"LWdown  Precip_CRU  Precip_GPCC  PSurf  Qair  Rainf_CRU  Rainf_GPCC  Rn  Snowf_CRU  Snowf_GPCC  SWdown  Tair  Wind".split(),
        "filepatt": "variable.year.nc",
        "t_res": "daily"
    },
    "ERAinterim": {
        "path":"/Net/Groups/BGI/data/DataStructureMDI/DATA/grid/Global/0d50_daily/ERAinterim/v2/Data/variable",
        "variables":["Tair","VPDday_1MJ","Precip","SWdown","SWnet","LWdown","Tmax","Tmin","Wind","RH","PSurf"],
        "filepatt": "variable.720.360.year.nc",
        "t_res": "daily"
        },
    "ERA5": {
        "path":"/Net/Groups/data_BGC/era5/e1/0d25_hourly/variable/year/",
        "variables":["t2m","q","tp","ssrd","ssr","strd","ws10","rH","sp"],
        "filepatt": "variable.*.01.year.nc",
        "t_res": "hourly"
        },
    "fluxnetBGI2021.LaThuile.DD":{
        "path": f"/Net/Groups/BGI/work_3/cliff/input/{out_dir_name}/{run_name}.LaThuile.DD/data/daily",
        "sites": fnsites, 
        "variables": fnvars,
        "filepatt": "*.nc",
        "t_res": "daily"
    },
    "fluxnetBGI2021.LaThuile.HR":{
        "path": f"/Net/Groups/BGI/work_3/cliff/input/{out_dir_name}/{run_name}.LaThuile.HR/data/hourly",
        "sites": fnsites, 
        "variables": fnvars,
        "filepatt": "*.nc",
        "t_res": "hourly"
    },
    "fluxnetBGI2021.FLUXNET2015.DD":{
        "path": f"/Net/Groups/BGI/work_3/cliff/input/{out_dir_name}/{run_name}.FLUXNET2015.DD/data/daily",
        "sites": fnsites, 
        "variables": fnvars,
        "filepatt": "*.nc",
        "t_res": "daily"
    },
    "fluxnetBGI2021.FLUXNET2015.HR":{
        "path": f"/Net/Groups/BGI/work_3/cliff/input/{out_dir_name}/{run_name}.FLUXNET2015.HR/data/hourly",
        "sites": fnsites, 
        "variables": fnvars,
        "filepatt": "*.nc",
        "t_res": "hourly"
    }

    }

out_dict = {}
sel_data = ["ERAinterim", "CRUNCEP"]
sel_data = [f"{run_name}.{data}" for data in ["FLUXNET2015.DD", "FLUXNET2015.HR", "LaThuile.DD", "LaThuile.HR"]]
# sel_data = list(data_sets.keys())
# sel_data = ["ERA5"]
year_ref = '2000'
# out_dict['reference_year_global_data'] = year_ref
# for data, dset in data_sets.items():
out_dict['datasets']={}
do_all_fn = False
for data in sel_data:
    out_dict['datasets'][data]={}
    dset = data_sets[data]
    out_dict['datasets'][data]['path']=os.path.abspath(os.path.join(dset['path'],'../../'))
    out_dict['datasets'][data]['variables'] = {}
    # print(data,dset)
    if "fluxnetBGI2021" not in data:
        for var_name in dset['variables']:
            filename = dset['filepatt'].replace('variable', var_name).replace('year',year_ref)
            datapath = dset['path'].replace('variable', var_name).replace('year',year_ref)
            datafile = os.path.join(datapath, filename)
            out_dict['datasets'][data]['variables'][var_name] = {}
            if data == 'ERA5':
                d_check = os.path.join(datapath, '../')
            else:
                d_check = datapath
            yrs = []
            for ff in os.listdir(d_check):
                print(d_check, ff)
                if ff.endswith('.nc'):
                    yradd = int(ff.split('.')[-2])
                else:
                    yradd = int(ff)
                
                if yradd > 1850:
                    yrs=np.append(yrs, yradd)
            
            out_dict['datasets'][data]['variables'][var_name]["year_start"] =  str(int(yrs.min()))
            out_dict['datasets'][data]['variables'][var_name]["year_stop"] = str(int(yrs.max()))
            if data == 'ERA5':
                dat=xr.open_mfdataset(datafile, decode_times=False, concat_dim=['time'])
            else:
                dat=xr.open_dataset(datafile, decode_times=False)
            out_dict['datasets'][data]['variables'][var_name]["ref_year_data"] = year_ref
            out_dict['datasets'][data]['variables'][var_name]["file_path"] = datafile
            out_dict['datasets'][data]['variables'][var_name]["units_in_nc_file"] = dat[var_name].attrs["units"]
            out_dict['datasets'][data]['variables'][var_name]["long_name"] = dat[var_name].attrs["long_name"]
            out_dict['datasets'][data]['variables'][var_name]["temporal_resolution"] = dset["t_res"]
            out_dict['datasets'][data]['variables'][var_name]["nan_perc_5"] = str(np.round(dat[var_name].reduce(np.nanpercentile, q=5).values,8))
            out_dict['datasets'][data]['variables'][var_name]["nan_perc_95"] = str(np.round(dat[var_name].reduce(np.nanpercentile, q=95).values,8))
            out_dict['datasets'][data]['variables'][var_name]["nan_max"] = str(np.round(dat[var_name].max(skipna=True).values,8))
            out_dict['datasets'][data]['variables'][var_name]["nan_mean"] = str(np.round(dat[var_name].mean(skipna=True).values,8))
            out_dict['datasets'][data]['variables'][var_name]["nan_min"] = str(np.round(dat[var_name].min(skipna=True).values,8))            # datval = dat[var_name].values
            # out_dict['datasets'][data]['variables'][var_name]["nan_perc_5"] = str(np.nanpercentile(datval, 5))
            # out_dict['datasets'][data]['variables'][var_name]["nan_perc_95"] = str(np.nanpercentile(datval, 95))
            # out_dict['datasets'][data]['variables'][var_name]["nan_max"] = str(np.nanmax(datval))
            # out_dict['datasets'][data]['variables'][var_name]["nan_min"] = str(np.nanmin(datval))
            # print(dat[var_name].attrs["units"], dat[var_name].attrs["long_name"], out_dict['datasets'][data]['variables'][var_name]["nan_max"], out_dict['datasets'][data]['variables'][var_name]["nan_min"])
            dat.close()
            with open('global_data_info.json', 'w') as fp:
                json.dump(out_dict, fp, indent=2, sort_keys=True)
    else:
        if 'DD' in data:
            timesc='daily'
        else:
            timesc='hourly'
        filenames = [f'{site}.1989.2020.{timesc}.nc' for site in dset['sites']]
        filepaths = [os.path.join(dset['path'], f'{site}.1989.2020.{timesc}.nc') for site in dset['sites']]
        if do_all_fn:
            filepaths = []
            sel_sites = []
            for ff in os.listdir(dset['path']):
                print(dset['path'], ff)
                if ff.endswith('.nc'):
                    fileadd = ff
                    filepaths=np.append(filepaths, os.path.join(dset['path'],fileadd))
                    sel_sites = np.append(sel_sites, ff.split('.')[0])
            dset['sites'] = ["all"]
            # dset['sites'] = sel_sites

        print(dset['sites'], filepaths)        
        # kera

        dat = xr.open_mfdataset(filepaths, concat_dim=['latitude', 'longitude'], data_vars=dset['variables'])
        print (dat, filepaths)
        for var_name in dset['variables']:
            out_dict['datasets'][data]['variables'][var_name] = {}
            out_dict['datasets'][data]['variables'][var_name]["year_start"] = '1989'
            out_dict['datasets'][data]['variables'][var_name]["year_stop"] = '2020'
            out_dict['datasets'][data]['variables'][var_name]["file_path"] = ','.join(dset['sites'])
            out_dict['datasets'][data]['variables'][var_name]["temporal_resolution"] = timesc
            out_dict['datasets'][data]['variables'][var_name]["units_in_nc_file"] = dat[var_name].attrs["units"]
            out_dict['datasets'][data]['variables'][var_name]["long_name"] = dat[var_name].attrs["long_name"]
            # dat_5 = dat[var_name].reduce(np.nanpercentile, q=5).values
            # print(dat_5)
            out_dict['datasets'][data]['variables'][var_name]["nan_perc_5"] = str(np.round(dat[var_name].reduce(np.nanpercentile, q=5).values,8))
            out_dict['datasets'][data]['variables'][var_name]["nan_perc_95"] = str(np.round(dat[var_name].reduce(np.nanpercentile, q=95).values,8))
            out_dict['datasets'][data]['variables'][var_name]["nan_max"] = str(np.round(dat[var_name].max(skipna=True).values,8))
            out_dict['datasets'][data]['variables'][var_name]["nan_mean"] = str(np.round(dat[var_name].mean(skipna=True).values,8))
            out_dict['datasets'][data]['variables'][var_name]["nan_min"] = str(np.round(dat[var_name].min(skipna=True).values,8))
            # datval = dat[var_name].values
            # out_dict['datasets'][data]['variables'][var_name]["nan_perc_5"] = str(np.nanpercentile(datval, 5))
            # out_dict['datasets'][data]['variables'][var_name]["nan_perc_95"] = str(np.nanpercentile(datval, 95))
            # out_dict['datasets'][data]['variables'][var_name]["nan_max"] = str(np.nanmax(datval))
            # out_dict['datasets'][data]['variables'][var_name]["nan_min"] = str(np.nanmin(datval))
            out_dict['datasets'][data]['variables'][var_name]["ref_year_data"] = '1989-2020'
            with open('global_data_info.json', 'w') as fp:
                json.dump(out_dict, fp, indent=2, sort_keys=True)
        dat.close()


    with open('global_data_info.json', 'w') as fp:
        json.dump(out_dict, fp, indent=2, sort_keys=True)





