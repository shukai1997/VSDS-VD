import pandas as pd
import os
import subprocess
import glob
from multiprocessing import Pool
from wrapt_timeout_decorator import timeout
os.environ["TA_LICENSE"] = "/home/shukai/expanding_boundary/models/pyhsical-based-docking-tools/Surflex/tools/sybylx2.1.1/AdminTools11.6"

#Please adjust your file path and sorfware path

def pdb_to_mol2(pdb_file:str, output_file:str):
    # if os.path.exists(output_file):
    #     print(f'{output_file}: Exist!')
    #     return output_file
    obabel = 'LD_PRELOAD=/opt/openbabel/3.1.1/lib64/libcoordgen.so.2 /opt/openbabel/3.1.1/bin/obabel'
    cmd = f'{obabel} {pdb_file} -O {output_file}'
    process = subprocess.Popen(
        cmd,
        shell=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    p = process.communicate()
    if os.path.exists(output_file):
        print(f'{output_file}: Have created!')
    else:
        print(f'{output_file}: not created!')
        print(p)

def get_pocket_tripossybylx(protein_prepared_file:str,ligand_for_center_file:str, output_file:str, verbose:bool=True):
    # if os.path.exists(output_file):
    #     print(f'{output_file}: Exists!')
    #     return output_file
    command =f'/home/shukai/expanding_boundary/models/pyhsical-based-docking-tools/Surflex/tools/sybylx2.1.1/sybylx2.1.1/bin/linuxem64/surflex-dock.exe proto {ligand_for_center_file} {protein_prepared_file} {output_file}'
    process = subprocess.Popen(
        command,
        shell=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    p = process.communicate()
    if verbose:
        if os.path.exists(output_file):
            print('successfully prepare the pocket')
        else:
            print('failed to prepare the pocket')
            print(command)
            print(p)
    return output_file

# protein_file = '/home/shukai/expanding_boundary/models/pyhsical-based-docking-tools/Surflex/test_surflex/6PYR_optimal.pdb'
# protein_prepared_file = '/home/shukai/expanding_boundary/models/pyhsical-based-docking-tools/Surflex/test_surflex/6PYR_optimal.mol2'
# pdb_to_mol2(protein_file, protein_prepared_file)
# ligand_for_center_file = '/home/shukai/expanding_boundary/models/pyhsical-based-docking-tools/Surflex/test_surflex/crystal_ligand.mol2'
# output_file = '/home/shukai/expanding_boundary/models/pyhsical-based-docking-tools/Surflex/test_surflex/6PYR_prepare.mol2'
# get_pocket_tripossybylx(protein_prepared_file,ligand_for_center_file, output_file, verbose=True)

def write_file(output_file, outline):
    buffer = open(output_file, 'w')
    buffer.write(outline)
    buffer.close()
    if os.path.exists(output_file):
        print(f'{output_file}  Have been written!')
    else:
        print(f'{output_file}: Failed, not created!')

@timeout(240)
def sybylx_dock(ligand_for_docking:str, protomol_file:str, protein_prepared_file:str,output_file:str, n:int=100,verbose:bool=True):
    # if os.path.exists(output_file):
    #     print(f'{output_file}: Exists!')
    #     return output_file
    list_file = os.path.join(os.path.dirname(output_file),'list')
    write_file(outline=f'{ligand_for_docking}', output_file=list_file)
    command =f'/home/shukai/expanding_boundary/models/pyhsical-based-docking-tools/Surflex/tools/sybylx2.1.1/sybylx2.1.1/bin/linuxem64/surflex-dock.exe -ndock_final {n} -multistart 10 -maxrot 100 -protodthresh 15.0 +fastsearch -maxconfs 20 dock_list {list_file} {protomol_file} {protein_prepared_file} {output_file}'
    process = subprocess.Popen(
        command,
        shell=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    p = process.communicate()
    if verbose:
        if os.path.exists(output_file):
            print('successfully Surflex docking')
        else:
            print('failed to Surflex docking')
            print(command)
            print(p)
    
    return output_file

# ligand_for_docking = '/home/shukai/expanding_boundary/models/pyhsical-based-docking-tools/Surflex/test_surflex/crystal_ligand.mol2'
# protomol_file = '/home/shukai/expanding_boundary/models/pyhsical-based-docking-tools/Surflex/test_surflex/6PYR_prepare.mol2-marked-protein.pdb'
# protein_prepared_file = '/home/shukai/expanding_boundary/models/pyhsical-based-docking-tools/Surflex/test_surflex/6PYR_optimal.mol2'
# output_file = '/home/shukai/expanding_boundary/models/pyhsical-based-docking-tools/Surflex/test_surflex/result.log'
# sybylx_dock(ligand_for_docking, protomol_file, protein_prepared_file, output_file, 100, True)

def split_mol2(multi_mol2 , output_file):
    obabel = 'LD_PRELOAD=/opt/openbabel/3.1.1/lib64/libcoordgen.so.2 /opt/openbabel/3.1.1/bin/obabel'
    cmdline = '%s %s -O %s -m'%(obabel , multi_mol2 , output_file)
    os.system(cmdline)

# multi_mol2 = '/home/shukai/expanding_boundary/models/pyhsical-based-docking-tools/Surflex/test_surflex/result.log-results.mol2'
# output_file = '/home/shukai/expanding_boundary/models/pyhsical-based-docking-tools/Surflex/test_surflex/docking_results.mol2'
# split_mol2(multi_mol2 , output_file)

def log2csv(log_file , csv_file):
    f = open(log_file,'r')
    cons = f.readlines()
    num = len(cons)-1
    data = pd.DataFrame()
    data['score'] = [float(cons[1].split(' ')[5])]
    data.to_csv(csv_file , index=False)
    return num

# log_file = '/home/shukai/expanding_boundary/models/pyhsical-based-docking-tools/Surflex/test_surflex/result.log'
# csv_file = '/home/shukai/expanding_boundary/models/pyhsical-based-docking-tools/Surflex/test_surflex/result.csv'
# log2csv(log_file , csv_file)

def get_pdb_file(path):
    for file in os.listdir(path):
        if len(file.split('.')) == 2 and file.split('.')[1] == 'pdb':
            return file[:4]


def one_dragon_service(uniprot):

    for sdf_name in os.listdir('/home/shukai/expanding_boundary/models/pyhsical-based-docking-tools/Surflex/data/%s/inactives'%uniprot):

        path_output = '/home/shukai/expanding_boundary/models/pyhsical-based-docking-tools/Surflex/data/%s/inactives/%s'%(uniprot, sdf_name)
        if not os.path.exists(os.path.join(path_output , 'docking_result1.mol2')):
            print(path_output)
            ligand_for_docking_sdf = os.path.join(path_output , '%s.sdf'%sdf_name)
            ligand_for_docking_mol2 = os.path.join(path_output , '%s.mol2'%sdf_name)
            sdf_to_mol2(ligand_for_docking_sdf , ligand_for_docking_mol2)

            path_root = '/home/shukai/expanding_boundary/models/pyhsical-based-docking-tools/Surflex/data/%s'%uniprot
            pdb_id = get_pdb_file(path_root)
            protomol_file = os.path.join(path_root , '%s_prepare.mol2-protomol.mol2'%pdb_id)
            protein_prepared_file = os.path.join(path_root , '%s_optimal.mol2'%pdb_id)

            
            output_file = os.path.join(path_output , 'result.log')
            try:
                sybylx_dock(ligand_for_docking_mol2, protomol_file, protein_prepared_file, output_file, 100, True)
           
                multi_mol2 = os.path.join(path_output , 'result.log-results.mol2')
                output_file = os.path.join(path_output , 'docking_result.mol2')
                split_mol2(multi_mol2 , output_file)
                log_file = os.path.join(path_output , 'result.log')
                csv_file = os.path.join(path_output , 'result.csv')
                log2csv(log_file , csv_file)    
                
                for i in range(11,101):
                    os.remove(os.path.join(path_output , 'docking_result%s.mol2'%i))
                os.remove(os.path.join(path_output , 'result.log-results.mol2'))
            
            except:
                pass


def sdf_to_mol2(sdf_file:str, mol2_file:str):

    obabel = 'LD_PRELOAD=/opt/openbabel/3.1.1/lib64/libcoordgen.so.2 /opt/openbabel/3.1.1/bin/obabel'
    cmd = f'{obabel} {sdf_file} -O {mol2_file}'
    process = subprocess.Popen(
        cmd,
        shell=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    p = process.communicate()

    if os.path.exists(mol2_file):
        print(f'{mol2_file}: Have created!')
    else:
        print(f'{mol2_file}: not created!')
        print(p)

import random

def shuffle_list(input_list):
    shuffled_list = input_list[:]
    random.shuffle(shuffled_list)
    return shuffled_list

if __name__ == '__main__':

    #受体准备

    #uniprot = 'O00329'
    pool = Pool(48)
    pool.map(one_dragon_service, shuffle_list(os.listdir('/home/shukai/expanding_boundary/models/pyhsical-based-docking-tools/Surflex/data')))
    pool.close()
    pool.join()




