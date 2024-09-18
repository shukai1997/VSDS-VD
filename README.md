# VSDS-VD: benchmarking AI-powered docking methods from the perspective of virtual screening

![](https://github.com/shukai1997/VSDS-VD/blob/main/Figure1.jpg)

## Contents

- [Overiew](#overview)
- [Obtain Datasets](#obtain-datasets)
- [Docking and Scoring](#docking-and-scoring)
- [Post-processing](#post-processing)
- [Metrics Calculation](#demo--reproduction-ligand-docking-on-pdbbind-core-set)

## Overview

In this work, we conducted a comprehensive benchmark analysis of four AI-powered, four physics-based docking tools, and two AI-powered re-scoring methods. We initially constructed DTEBV-D, on which the re-docking experiments reveal that KarmaDock and CarsiDock surpassed all physics-based tools on docking accuracy while all physics-based tools significantly outperformed AI-based methods on structural rationality. The VS results on DTEBV-D highlights the effectiveness of RTMScore as a re-score function and Glide-based methods achieved the highest enrichment factors (EFs) among all docking tools. We additionally constructed DRSM-D that more closely resembles real VS scenarios,where the employed AI-based tools obviously outperformed Glide. Finally, we proposed a hierarchical VS strategy that could efficiently and accurately enrich active molecules in real large-scale VS projects.

## Obtain Datasets

The prepared benchmarking set is available at https://zenodo.org/records/13684010.

## Docking and Scoring

Four AI-powered docking tools (CarsiDock, KarmaDock, DiffDock, and FlexPose), four physics-based tools (Glide, Surflex, rDock, and LeDock) and two AI-powered re-scoring methods (RTMScore and EquiScore) were employed on three benchmarks. *The demo scripts are provided.*

#### AI-powered docking tools

The AI-powered docking  and rescoring protocols can be used from their specific Github repository.

```
CarsiDock : https://github.com/carbonsilicon-ai/CarsiDock
KarmaDock : https://github.com/schrojunzhang/KarmaDock 
FlexPose : https://github.com/tiejundong/FlexPose
DiffDock : https://github.com/gcorso/DiffDock
RTMScore : https://github.com/sc8668/RTMScore
EquiScore : https://github.com/CAODH/EquiScore
```

#### physics-based tools

The use of physics-based tools requires the installation of the corresponding version of the software, and *the demo scripts are provided.*

```
Glide : schrodinger v2022-1
LeDock : LeDock v1.0
Surflex : surflex v2.706
rDock : rDock r2021.1
```

## Post-processing

Two post-processing methods, i.e., Force Field optimization (FF) and RDKit-conformation alignment (Align) are employed in this work to explore their role in re-docking experiments and VS tasks. FF optimization can be easily implemented by MMFFOptimizeMolecule, an Application Programming Interface (API) of RDKit, which uses the MMFF94 Force Field to optimize the predicted ligand conformations with a maximum of ten steps. Align optimization involves aligning the chemically plausible RDKit conformations, which lack accurate orientation and torsional angles, with the predicted binding poses and replacing the initially predicted pose with the appropriately aligned RDKit conformations. Specifically, the ligands' torsional angles were first defined and then assigned to the RDKit conformations. Finally, the modified RDKit conformations were aligned with the predicted binding poses using the Kabsch algorithm.

e.g.

```
cd /root/VSDS-VD/post_process 
python -u post_processing.py --true_ligand ./true_ligand.mol2 --predicted_ligand ./predicted_ligand.sdf --method align
```

## Metrics Calculation

In re-docking experiments, we analyzed the docking accuracy and physical plausibility of the top-ranked conformation generated by different docking methods. The docking accuracy was measured by the success rate of generating binding conformations with an RMSD of less than 2 Å from the ground-truth conformations. The physical plausibility of predicted conformation was analyzed by PoseBusters, a Python package that utilizes the well-established cheminformatics toolkit RDKit to conduct a series of standard quality assessments on docked conformations. In VS tasks, the enriching power of different protocols is mainly reflected by EFs at different thresholds (0.5%, 1%, and 5%).

e.g. 
(*remember to install posebuster and oddt package*)

Docking

```
cd /root/VSDS-VD/metric_calculation 
python -u metric_calculation.py --task_type docking --ligand_true ./ligand_true.sdf --ligand_predict ./ligand_predict.sdf --protein_file ./protein.pdb
```

Virtual Screening

```
cd /root/VSDS-VD/metric_calculation 
python -u metric_calculation.py --task_type vs --vs_result_csv_path ./vs_result.csv

