import os
import pandas as pd


def get_pdb_id(path):
    for file in os.listdir(path):
        if 'optimal.pdb' in file:
            return file[:4]

path = '/home/shukai/expanding_boundary/models/EquiScore-main/shukai_data/glide'
for file in os.listdir(path):
    path1 = os.path.join(path , file)
    pdb_file = os.path.join(path1 , '%s_optimal.pdb'%get_pdb_id(path1))
    path_active = os.path.join(path1 , 'actives_protein_noligand')
    inactives_sdf = '/home/shukai/expanding_boundary/data/rtmscore/targets_actives/%s/docking_molecules_10poses.sdf'%file 
    cmdline = 'python ./get_pocket/get_pocket.py --docking_result %s \
                --recptor_pdb %s \
                --single_sdf_save_path %s \
                --pocket_save_dir %s'%(inactives_sdf , pdb_file , os.path.join(path_active , 'tmp_sdfs') , os.path.join(path_active , 'tmp_pockets'))
    os.system(cmdline)



path = '/home/shukai/expanding_boundary/models/EquiScore-main/shukai_data/glide'
for file in os.listdir(path):
    path1 = os.path.join(path , file)
    path_active = os.path.join(path ,file , 'actives_protein_noligand')
    cmdline = 'python Screening.py --ngpu 1 --test --test_path %s \
                --test_name tmp_pockets \
                --MASTER_PORT 29506 \
            --pred_save_path  %s'%(path_active  , os.path.join(path_active , 'result.csv'))
    os.system(cmdline)
