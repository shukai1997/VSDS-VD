import os

protein_path ='/home/shukai/expanding_boundary/result/posebuster/poses/O00329/6PYR_optimal.pdb'
ligand_path = '/home/shukai/expanding_boundary/models/CarsiDock-main/shukai_data_vsResult_dealt/O00329/actives/actives.sdf'
refer_ligand_path = '/home/shukai/expanding_boundary/result/posebuster/poses/O00329/crystal_ligand.sdf'
output_file = '/home/shukai/expanding_boundary/models/CarsiDock-main/shukai_data_vsResult_dealt/O00329/actives/actives.csv'
command = f'python3 rtmscore.py -p {protein_path} \
        -l {ligand_path} \
        -rl {refer_ligand_path} \
        -o  {output_file}\
        -gen_pocket -c 10.0 -m ../trained_models/rtmscore_model1.pth'
os.system(command)