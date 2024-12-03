import subprocess
import os
import glob

#Please adjust your file path and sorfware path


def mol2_to_sdf(mol2_file:str, output_file:str):
    if os.path.exists(output_file):
        print(f'{output_file}: Exist!')
        return output_file
    obabel = 'LD_PRELOAD=/opt/openbabel/3.1.1/lib64/libcoordgen.so.2 /opt/openbabel/3.1.1/bin/obabel'
    cmd = f'{obabel} {mol2_file} -O {output_file}'
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

def pdb_to_mol2(pdb_file:str, output_file:str):
    if os.path.exists(output_file):
        print(f'{output_file}: Exist!')
        return output_file
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

def write_file(output_file, outline):
    buffer = open(output_file, 'w')
    buffer.write(outline)
    buffer.close()
    if os.path.exists(output_file):
        print(f'{output_file}  Have been written!')
    else:
        print(f'{output_file}: Failed, not created!')
def cavity_prepare_xiasha(protein_mol2_file:str, ref_ligand_sdf_file:str, output_file:str, verbose:bool=True):
    prm_file = os.path.join(os.path.dirname(output_file),"cavity.prm")
    if os.path.exists(output_file):
        print(f'{output_file}  Exists!')
        return output_file,prm_file
    outfile = f'''RBT_PARAMETER_FILE_V1.00
TITLE {os.path.basename(protein_mol2_file).split('-')[0]}_dock

RECEPTOR_FILE {protein_mol2_file}
### RECEPTOR_FLEX 3.0

##################################################################
### CAVITY DEFINITION: REFERENCE LIGAND METHOD
##################################################################
SECTION MAPPER
    SITE_MAPPER RbtLigandSiteMapper
    REF_MOL {ref_ligand_sdf_file}
        RADIUS 6.0
        SMALL_SPHERE 1.5
##      LARGE_SPHERE 4.0
        MAX_CAVITIES 1
        MIN_VOLUME 100
        VOL_INCR 0.0
        GRIDSTEP 0.5
END_SECTION

#################################
#CAVITY RESTRAINT PENALTY
#################################
SECTION CAVITY
    SCORING_FUNCTION RbtCavityGridSF
    WEIGHT 1.0
END_SECTION
'''
    write_file(prm_file, outfile)
    command = f'cd {os.path.dirname(output_file)} && /home/shukai/expanding_boundary/models/pyhsical-based-docking-tools/rDock/rDock_2022_src/bin/rbcavity -was -d -r {prm_file}'
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
          print(f'{output_file} create Successfully!')
      else:
          print(f'{output_file} failed')
          print(command)
          print(p)
    return output_file,prm_file

def rdock_dock_xiasha(prm_file:str, ligand_for_docking_sdf_file:str, output_file:str, n:int, verbose:bool=False):
    # if os.path.exists(output_file):
    #     print(f'{output_file}  Exists!')
    #     return output_file
    command = f'cd {os.path.dirname(output_file)} && /home/shukai/expanding_boundary/models/pyhsical-based-docking-tools/rDock/rDock_2022_src/bin/rbdock -i \
                {ligand_for_docking_sdf_file} -o {output_file[:-4]} -r {prm_file} -p \
                 /home/shukai/expanding_boundary/models/pyhsical-based-docking-tools/rDock/rDock_2022_src/data/scripts/dock.prm -n {n}'
    proc = subprocess.Popen(
        command,
        shell=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    proc.communicate()
    cmd = f'/home/shukai/expanding_boundary/models/pyhsical-based-docking-tools/rDock/rDock_2022_src/bin/sdreport -l {output_file[:-3]}sd > {output_file}'
    proc = subprocess.Popen(
        cmd,
        shell=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    p = proc.communicate()
    if verbose:
      if os.path.exists(output_file):
          print(f'rdock Successfully!')
      else:
          print(f'rdock failed')
          print(command)
          print(p)
    #multi_sd_to_mol2(f'{output_file[:-3]}sd',f'{output_file[:-4]}_rdock_.mol2')
    return output_file





def one_dragon_service(uniprot):
    pdbid_dir = '/home/shukai/expanding_boundary/models/pyhsical-based-docking-tools/rDock/data/%s'%uniprot
    prm_file = os.path.join(pdbid_dir , 'cavity.prm')
    for sdf_name in os.listdir(os.path.join(pdbid_dir , 'inactives')):
        ligand_docking_prepared_file = os.path.join(pdbid_dir , 'inactives' , sdf_name , '%s.sdf'%sdf_name )
        score_file = os.path.join(pdbid_dir , 'inactives' , sdf_name , 'score.txt' )
        score_sd = os.path.join(pdbid_dir , 'inactives' , sdf_name , 'score.sd' )
        if not os.path.exists(score_sd):
            print(score_sd)
            # try:
            rdock_dock_xiasha(prm_file, ligand_docking_prepared_file,score_file, 10, verbose=True)
           # except:
                #pass

import random

def shuffle_list(input_list):
    shuffled_list = input_list[:]
    random.shuffle(shuffled_list)
    return shuffled_list
if __name__ == '__main__':


    from multiprocessing import Pool
    # pool = Pool(48)
    # pool.map(one_dragon_service, os.listdir('/home/shukai/expanding_boundary/models/pyhsical-based-docking-tools/rDock/data'))
    # pool.close()
    # pool.join()
    path = '/home/shukai/expanding_boundary/models/pyhsical-based-docking-tools/rDock/data'
    for file in shuffle_list(os.listdir(path)):
        one_dragon_service(file)
