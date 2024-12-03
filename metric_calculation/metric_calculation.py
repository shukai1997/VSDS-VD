import os
import pandas as pd
import numpy as np
import argparse
os.environ["BABEL_LIBDIR"] = "/home/shukai/miniconda3/envs/pbcnet/lib/openbabel/3.1.1"
from oddt import metrics

#Please pip install Posebuster before use the code

def docking_result_dealt(ligand_true, ligand_predict , protein_file):
    cmd = 'bust %s -l %s -p %s --outfmt csv --output docking_result.csv'%(ligand_predict, ligand_true , protein_file)
    os.system(cmd)
    
def vs_result_dealt(vs_result_csv_path):
    data = pd.read_csv(vs_result_csv_path)
    score = list(data['score'])
    x = np.array(score)
    y = data['label']
    y = np.array(y)
    ef_0_5 = metrics.enrichment_factor(y, x, percentage=0.5, pos_label=1, kind='fold')
    ef_1 = metrics.enrichment_factor(y, x, percentage=1, pos_label=1, kind='fold')
    ef_5 = metrics.enrichment_factor(y, x, percentage=5, pos_label=1, kind='fold')
    return  ef_0_5, ef_1, ef_5 

if __name__ == "__main__":

    argparser = argparse.ArgumentParser()
    
    argparser.add_argument('--task_type', 
                        default='docking',
                        help='the task type')    
    
    argparser.add_argument('--ligand_true', 
                        default='./ligand_true.sdf',
                        help='the crystal ligand')
    argparser.add_argument('--ligand_predict', 
                        default='./ligand_predict.sdf',
                        help='the predicted ligand conformation')
    argparser.add_argument('--protein_file', 
                            default='./protein.pdb',
                            help='the protein file')
    argparser.add_argument('--vs_result_csv_path', 
                            default='./vs_result.csv',
                            help='the vs result')

    args = argparser.parse_args()


    if args.task_type == 'docking':
        docking_result_dealt(args.ligand_true, args.ligand_predict , args.protein_file)
    else:
        vs_result_dealt(args.vs_result_csv_path)