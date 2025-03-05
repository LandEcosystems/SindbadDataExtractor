import os
import csv
def get_cliff_dirname(fnversion, temp_reso):
    """
    Get the folder structure of cliff input and output based on a simple convention
    """
    if fnversion == 'FLUXNET2015':
        opath = 'BRK15'
    elif fnversion == 'LaThuile':
        opath = 'LTL07'
    else:
        opath = fnversion

    if temp_reso == 'daily':
        otime = 'DD'
    else:
        otime = 'HR'
    return f'fluxnetBGI2021.{opath}.{otime}'


output_path='/Net/Groups/BGI/scratch/skoirala/fluxNet_0.06/'
## declare an array variable
fn_versions=("WarmWinter2020", "AmerifluxFLUXNET", "ICOSFLUXNET", "ICOSDrought2018", "FLUXNET2015", "CCBY", "CCBY_v0.05", "CCBY_v0.04")
site_list=[]
for fn_v in fn_versions:
    dname = get_cliff_dirname(fn_v, 'hourly') 
    fname = f'site_list_{fn_v}.csv'
    site_file = os.path.join(output_path, dname, 'config_info', fname)
    print(site_file)
    if os.path.exists(site_file):
        with open(site_file, newline='') as csvfile:
            spamreader = csv.reader(csvfile, delimiter=' ')
            for row in spamreader:
                site = row[0]
                if site not in site_list:
                    site_list.append(site)

site_list.sort()
with open('all_site_list.txt', 'w') as f:
    f.write(f'{len(site_list)} sites \n')
    f.write(', '.join(site_list))

