import pandas as pd
import os
import subprocess
import glob
import csv
from wrapt_timeout_decorator import timeout

def prepare_receptor_local(protein_file:str, output_file:str,verbose:bool=True):
    output_dir = os.path.dirname(output_file)
    command = f'cd {output_dir} && /home/ruolan/shukai/ledock/tools/LeDock/lepro {protein_file}'
    if os.path.exists(output_file):
        print(f'{output_file}: Exist!')
        return output_file
    proc = subprocess.Popen(
        command,
        shell=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    p = proc.communicate()
    if verbose:
        if os.path.exists(output_file):
            print('successfully prepare the protein')
        else:
            print('failed to prepare the protein')
            print(command)
            print(p)
    return output_file

def pdb_to_mol2(pdb_file:str, output_file:str):

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

def sdf_to_mol2(pdb_file:str, output_file:str):
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

def get_mol_center_mol2(ligand_mol2):
    x = os.popen(
        "cat %s |sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p |awk '{{print $3}}' | awk '{{x+=$1}} END {{print x/(NR-2)}}'" % ligand_mol2).read()
    y = os.popen(
        "cat %s |sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p |awk '{{print $4}}' | awk '{{y+=$1}} END {{print y/(NR-2)}}'" % ligand_mol2).read()
    z = os.popen(
        "cat %s |sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p |awk '{{print $5}}' | awk '{{z+=$1}} END {{print z/(NR-2)}}'" % ligand_mol2).read()
    print(f'Center of {os.path.basename(ligand_mol2)} is: {x},{y},{z}')
    return float(x.strip()), float(y.strip()), float(z.strip())

def get_mol_size_mol2(ligand_mol2,grid_plus:float = 15,grid_multi:float = 1):
    x_max = os.popen(
        "cat %s |sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p |awk '{{print $3}}' |awk 'BEGIN {max = 0} {if ($1 > max) max = $1} END {{print max}}'" % ligand_mol2).read()
    x_min = os.popen(
        "cat %s |sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p |awk '{{print $3}}' |awk '{if ($1 != \"\") {print $1}}' |awk 'BEGIN {min = 9999999} {if ($1+0 < min+0) min = $1 fi} END {{print min}}'" % ligand_mol2).read()
    x_size = float(x_max.strip()) - float(x_min.strip())*grid_multi + grid_plus
    y_max = os.popen(
        "cat %s |sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p |awk '{{print $4}}' |awk 'BEGIN {max = 0} {if ($1 > max) max = $1} END {{print max}}'" % ligand_mol2).read()
    y_min = os.popen(
        "cat %s |sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p |awk '{{print $4}}' |awk '{if ($1 != \"\") {print $1}}' |awk 'BEGIN {min = 9999999} {if ($1+0 < min+0) min = $1 fi} END {{print min}}'" % ligand_mol2).read()
    y_size = float(y_max.strip()) - float(y_min.strip())*grid_multi + grid_plus
    z_max = os.popen(
        "cat %s |sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p |awk '{{print $5}}' |awk 'BEGIN {max = 0} {if ($1 > max) max = $1} END {{print max}}'" % ligand_mol2).read()
    z_min = os.popen(
        "cat %s |sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p |awk '{{print $5}}' |awk '{if ($1 != \"\") {print $1}}' |awk 'BEGIN {min = 9999999} {if ($1+0 < min+0) min = $1 fi} END {{print min}}'" % ligand_mol2).read()
    z_size = float(z_max.strip()) - float(z_min.strip())*grid_multi + grid_plus
    print(f'Size of {os.path.basename(ligand_mol2)} is: {x_size},{y_size},{z_size}')
    return x_size, y_size, z_size

def dir_check(dir:str, create_dir:bool = True):
    if os.path.exists(dir):
        print(f'{dir}   Exists')
    else:
        if create_dir:
            os.system(f'mkdir -p {dir}')
            print(f'{dir}  Create Already!')
        else: 
            print(f'{dir}  Not Exists Yet')     
def cp(file_path:str, target_path:str):
    os.system(f'cp -r {file_path} {target_path}')

def write_file(output_file, outline):
    buffer = open(output_file, 'w')
    buffer.write(outline)
    buffer.close()
    if os.path.exists(output_file):
        print(f'{output_file}  Have been written!')
    else:
        print(f'{output_file}: Failed, not created!')

@timeout(200)
def ledock_dock_local(protein_prepared_file:str, ligand_for_center:str, ligand_for_docking:str, output_file:str, n:int,  verbose:bool=True):

    x,y,z = get_mol_center_mol2(ligand_for_center)
    a,b,c = get_mol_size_mol2(ligand_for_center, grid_plus=10)
    o = round(x - a/2, 3) 
    p = round(x + a/2, 3)
    q = round(y - b/2, 3)
    r = round(y + b/2, 3)
    s = round(z - c/2 ,3)
    t = round(z + c/2 ,3)
    outline1 = f'''{os.path.basename(ligand_for_docking)}'''
    write_file(os.path.join(os.path.dirname(output_file),'ligands'),outline1)
    outline2 = f'''Receptor
{protein_prepared_file}

RMSD
1.0

Binding pocket
{o} {p}
{q} {r}
{s} {t}

Number of binding poses
{n}

Ligands list
ligands

END


'''
    write_file(os.path.join(os.path.dirname(output_file),'dock.in'),outline2)
    command = f'cd {os.path.dirname(output_file)} && /home/ruolan/shukai/ledock/tools/LeDock/ledock dock.in '
    proc = subprocess.Popen(
        command,
        shell=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    p = proc.communicate()
    if verbose:
        if os.path.exists(output_file):
            print('successfully ledock')
        else:
            print('failed ledock')
            print(command)
            print(p)
    return output_file


import re

def extract_float(input_string):
    result = re.search(r'[-+]?\d*\.\d+', input_string)
    if result:
        return float(result.group())
    else:
        return None
def get_affinity(file):
    path_dir = os.path.dirname(file)
    f = open(file)
    cons = f.readlines()
    affinity = extract_float(cons[1])
    data = pd.DataFrame()
    data['affinity'] = [affinity]
    data.to_csv(os.path.join(path_dir , 'affinity.csv') ,index=False)


def split_molecules(input_file):
    # 读取输入文件
    with open(input_file, 'r') as file:
        content = file.readlines()
    
    file_directory = os.path.dirname(input_file)
    molecule_count = 0
    molecule = []

    for line in content:
        if line.strip() != 'END':
            molecule.append(line)
        else:
            molecule_file = os.path.join(file_directory, f'molecule_{molecule_count}.pdb')
            with open(molecule_file, 'w') as output_file:
                output_file.write(''.join(molecule))
            molecule_file_mol2 = os.path.join(file_directory, f'molecule_{molecule_count}.mol2')
            pdb_to_mol2(molecule_file,molecule_file_mol2)
            os.remove(molecule_file)
            molecule_count += 1
            molecule = []
    return molecule_count

def get_pdb_file(path):
    for file in os.listdir(path):
        if len(file.split('.')) == 2 and file.split('.')[1] == 'pdb':
            return file[:4]


def one_dragon_service(file):
# uniprot = 'O00329'
    #pdb_id = get_pdb_file('/home/ruolan/shukai/ledock/data/%s'%uniprot)
    # protein_file = '/home/ruolan/shukai/ledock/data/%s/%s_optimal.pdb'%(uniprot , pdb_id)
    uniprot = 'Q86TI2'
    protein_prepared_file = '/home/ruolan/shukai/ledock/data/%s/pro.pdb'%uniprot
    # prepare_receptor_local(protein_file , protein_prepared_file)

    
    path_ligand_root = '/home/ruolan/shukai/ledock/data/%s/inactives'%uniprot
    # for file in os.listdir(path_ligand_root):

    ligand_for_center = '/home/ruolan/shukai/ledock/data/%s/crystal_ligand.mol2'%uniprot
    ligand_for_docking_sdf = '/home/ruolan/shukai/ledock/data/%s/inactives/%s/%s.sdf'%(uniprot , file , file)
    ligand_for_docking_mol2 = '/home/ruolan/shukai/ledock/data/%s/inactives/%s/%s.mol2'%(uniprot , file , file)
    sdf_to_mol2(ligand_for_docking_sdf , ligand_for_docking_mol2)
    output_file = '/home/ruolan/shukai/ledock/data/%s/inactives/%s/%s.dok'%(uniprot , file , file)
    if not os.path.exists('/home/ruolan/shukai/ledock/data/%s/inactives/%s/molecule_0.mol2'%(uniprot , file)):
        #try:
            #ledock_dock_local(protein_prepared_file , ligand_for_center,ligand_for_docking_mol2, output_file,n=10)

        get_affinity(output_file)
        split_molecules(output_file)
        #except:
         #   pass

import random

def shuffle_list(input_list):
    shuffled_list = input_list[:]
    random.shuffle(shuffled_list)
    return shuffled_list

if __name__ == '__main__':

    #受体准备

    from multiprocessing import Pool
    pool = Pool()
    pool.map(one_dragon_service, os.listdir('/home/ruolan/shukai/ledock/data/Q86TI2/inactives'))
    pool.close()
    pool.join()
    # one_dragon_service('O00329')


