import os

def get_pdb_id(path):
    for file in os.listdir(path):
        if 'optimal.pdb' in file:
            return file[:4]

path = '/mnt/cfs/users/luohao/shukai/CarsiDock/chemdiv_vsDataset'
for uniprot in os.listdir(path):
    pdb_id = get_pdb_id(os.path.join(path , uniprot))
    protein_path = '/mnt/cfs/users/luohao/shukai/CarsiDock/chemdiv_vsDataset/%s/%s_optimal.pdb'%(uniprot , pdb_id)
    ligand_path = '/mnt/cfs/users/luohao/shukai/CarsiDock/chemdiv_vsDataset/%s/crystal_ligand.sdf'%uniprot
    actives = '/mnt/cfs/users/luohao/shukai/CarsiDock/chemdiv_vsDataset/%s/chemdiv_decoys_round2.sdf'%uniprot
    outdir = '/mnt/cfs/users/luohao/shukai/CarsiDock/chemdiv_vsDataset/%s/chemdiv_inactives'%uniprot
    cmdline = 'python run_screening.py --cuda_convert --pdb_file %s --reflig %s --ligands %s --output_dir %s'%(protein_path , ligand_path , actives,outdir)
    os.system(cmdline)
    