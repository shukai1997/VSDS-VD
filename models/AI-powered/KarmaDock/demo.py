import os

def get_pdb_file(path):
    for file in os.listdir(path):
        if len(file.split('.')) == 2 and file.split('.')[1] == 'pdb' and 'optimal' in file:
            return file[:4]
        
import random

def shuffle_list(input_list):
    shuffled_list = input_list[:]
    random.shuffle(shuffled_list)
    return shuffled_list


path = '/home/shukai/expanding_boundary/test/KarmaDock-main/chemdiv_vsDataset'
path_output = '/home/shukai/expanding_boundary/test/KarmaDock-main/chemdiv_vsDataset_result'
for uniprot in shuffle_list(os.listdir(path)):
    path1 =os.path.join(path , uniprot)
    pdb_id = get_pdb_file(path1)
    ligand_smi = os.path.join(path1 , 'chemdiv_decoys.smi')
    protein_file = os.path.join(path1 , '%s_optimal.pdb'%pdb_id)
    crystal_ligand_file = os.path.join(path1 , 'crystal_ligand.mol2')
    out_dir = os.path.join(path_output , uniprot)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    cmd = 'python -u virtual_screening.py --ligand_smi %s \
        --protein_file %s \
        --crystal_ligand_file %s \
        --out_dir %s \
        --score_threshold 0 --batch_size 64 --random_seed 2023 --out_uncoorected  --out_corrected'%(ligand_smi , protein_file , crystal_ligand_file , out_dir)
    os.system(cmd)