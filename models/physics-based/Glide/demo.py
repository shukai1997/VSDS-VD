import os
from multiprocessing import Pool

#Please adjust your file path

def get_pdbid(path):
    for file in os.listdir(path):
        if 'optimal.pdb' in file:
            return file

def submit_task(uniprot):
    path = '/home/shukai/expanding_boundary/data/targets_inactives'
    path1 = os.path.join(path , uniprot)
    pdb_file = get_pdbid(path1)[:4]

    if os.path.exists(os.path.join(path1, '%s_optimal.zip'%pdb_file)):

        cmdline = 'cd %s &&'%path1
        cmdline += 'bash vsw_1.sh'
        os.system(cmdline)

pool = Pool(48)
pool.map(submit_task, os.listdir('/home/shukai/expanding_boundary/data/targets_inactives'))
pool.close()
pool.join()