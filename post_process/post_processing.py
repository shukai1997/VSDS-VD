import copy
import numpy as np
import rmsd
from rdkit import Chem, RDLogger
from rdkit import Geometry
from rdkit.Chem import AllChem, rdMolTransforms
from scipy.optimize import differential_evolution
import argparse
import os



def set_rdkit_mol_position(rdkit_mol, position):
    for j in range(position.shape[0]):
        rdkit_mol.GetConformer().SetAtomPosition(j,
                                            Geometry.Point3D(*position[j]))
    return rdkit_mol

def get_torsions(mol_list):
    # print('USING GEOMOL GET TORSIONS FUNCTION')
    atom_counter = 0
    torsionList = []
    for m in mol_list:
        torsionSmarts = '[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]'
        torsionQuery = Chem.MolFromSmarts(torsionSmarts)
        matches = m.GetSubstructMatches(torsionQuery)
        for match in matches:
            idx2 = match[0]
            idx3 = match[1]
            bond = m.GetBondBetweenAtoms(idx2, idx3)
            jAtom = m.GetAtomWithIdx(idx2)
            kAtom = m.GetAtomWithIdx(idx3)
            for b1 in jAtom.GetBonds():
                if (b1.GetIdx() == bond.GetIdx()):
                    continue
                idx1 = b1.GetOtherAtomIdx(idx2)
                for b2 in kAtom.GetBonds():
                    if ((b2.GetIdx() == bond.GetIdx())
                            or (b2.GetIdx() == b1.GetIdx())):
                        continue
                    idx4 = b2.GetOtherAtomIdx(idx3)
                    # skip 3-membered rings
                    if (idx4 == idx1):
                        continue
                    if m.GetAtomWithIdx(idx4).IsInRing():
                        torsionList.append(
                            (idx4 + atom_counter, idx3 + atom_counter, idx2 + atom_counter, idx1 + atom_counter))
                        break
                    else:
                        torsionList.append(
                            (idx1 + atom_counter, idx2 + atom_counter, idx3 + atom_counter, idx4 + atom_counter))
                        break
                break

        atom_counter += m.GetNumAtoms()
    return torsionList

def torsional_align(rdkit_mol, pred_conf, rotable_bonds):
    for rotable_bond in rotable_bonds:
        diheral_angle = GetDihedral(pred_conf, rotable_bond)
        SetDihedral(rdkit_mol.GetConformer(0), rotable_bond, diheral_angle)
    return rdkit_mol

def position_align_np(rdkit_mol, refer_mol, algo='kabsch'):
    A = rdkit_mol.GetConformer().GetPositions()
    B = refer_mol.GetConformer().GetPositions()
    B_center = rmsd.centroid(B)
    A -= rmsd.centroid(A)
    B -= B_center
    rmsd.quaternion_rotate
    if algo == 'kabsch':
        U = rmsd.kabsch(A, B)
    else: # quaternion
        U = rmsd.quaternion_rotate(A, B)
    A = np.dot(A, U)
    A += B_center
    set_rdkit_mol_position(rdkit_mol=rdkit_mol, position=A)

def GetDihedral(conf, atom_idx):
    return rdMolTransforms.GetDihedralRad(conf, atom_idx[0], atom_idx[1], atom_idx[2], atom_idx[3])
def SetDihedral(conf, atom_idx, new_vale):
    rdMolTransforms.SetDihedralRad(conf, atom_idx[0], atom_idx[1], atom_idx[2], atom_idx[3], new_vale)


def add_suffix_to_sdf_path(sdf_path , method):
   
    directory, filename = os.path.split(sdf_path)
    name, extension = os.path.splitext(filename)
    
  
    new_filename = name + '_%s'%method + extension
    
    new_path = os.path.join(directory, new_filename)
    
    return new_path

def correct_one(mol_true_path, mol_predict_path, method):
    # set pos
    mol = Chem.MolFromMol2File(mol_true_path)
    pos_pred = Chem.MolFromMolFile(mol_predict_path).GetConformer().GetPositions()
    raw_mol = copy.deepcopy(mol)
    pred_mol = copy.deepcopy(mol)
    pred_mol = set_rdkit_mol_position(pred_mol, pos_pred)
    # FF
    if method == 'ff':
        raw_mol = set_rdkit_mol_position(raw_mol, pos_pred)
        try:
            AllChem.MMFFOptimizeMolecule(raw_mol, maxIters=10)
        except:
            print('FF optimization failed')
    else:
        # get torsion_bonds
        rotable_bonds = get_torsions([pred_mol])
        # torsional align
        raw_mol = torsional_align(rdkit_mol=raw_mol, pred_conf=pred_mol.GetConformer(), rotable_bonds=rotable_bonds)
        # postion align
        # position_align_mol(rdkit_mol=mol, refer_mol=pred_mol)
        position_align_np(rdkit_mol=raw_mol, refer_mol=pred_mol)
        
    corrected_file = add_suffix_to_sdf_path(mol_predict_path , method)
    
    Chem.MolToMolFile(pred_mol, corrected_file)
    return raw_mol, pred_mol


if __name__ == "__main__":

    argparser = argparse.ArgumentParser()
    argparser.add_argument('--true_ligand', 
                        default='./true_ligand.mol2',
                        help='the crystal ligand')
    argparser.add_argument('--predicted_ligand', 
                        default='./predicted_ligand.sdf',
                        help='the predicted ligand conformation')
    argparser.add_argument('--method', 
                            default='align',
                            help='the optimal method')

    args = argparser.parse_args()


    correct_one(args.true_ligand, args.predicted_ligand, args.method)
   