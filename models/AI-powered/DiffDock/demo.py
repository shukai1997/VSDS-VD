import os
import random

def shuffle_list(input_list):
    shuffled_list = input_list.copy() 
    random.shuffle(shuffled_list)  
    return shuffled_list

path_result = '/home/shukai/expanding_boundary/models/DiffDock-main/results_chemdivVsDataset'

path = '/home/shukai/expanding_boundary/models/DiffDock-main/data/targets_chemdiv_vsDataset'
for file in shuffle_list(os.listdir(path)):
    uniprot = file.split('.')[0]
    path1 = os.path.join(path_result, uniprot)
    if not os.path.exists(path1):
        cmdline = 'cd /home/shukai/expanding_boundary/models/DiffDock-main &&'
        cmdline += 'python -m inference --protein_ligand_csv data/targets_chemdiv_vsDataset/%s.csv --out_dir results_chemdivVsDataset/%s --inference_steps 20 --samples_per_complex 10 \
        --batch_size 10 --actual_steps 18 --no_final_step_noise'%(uniprot, uniprot)
        os.system(cmdline)
