# Sr50-AvrSr50

## 1. Initial structural hypothesis generation

### Molecular docking simulations
The predicted structures of Sr50 and AvrSr50 can be found in [files](https://github.com/s-kyungyong/Sr50-AvrSr50/tree/main/files). These structures were submitted as a receptor and a ligand, respectively, to [ZDOCK](https://zdock.umassmed.edu/), [HDOCK](http://hdock.phys.hust.edu.cn/) and [ClusPro](https://cluspro.bu.edu/login.php). The default parameters were used, and no restraints were given. 

### Identification of the candidate poses
All docking models can be downloaded from [Zenodo](). There are two more input files needed to run this step: Tomborski_2022_Sr50_group_entropy.txt and Sr50.LRR.targets.txt. These files contain entropy calculation for each residue, which comes from [this publication](https://apsjournals.apsnet.org/doi/full/10.1094/MPMI-07-22-0154-R), and the LRR residues that compise inner concave of Sr50's LRR domain based on the predicted Sr50 structure. These files are available in [files](https://github.com/s-kyungyong/Sr50-AvrSr50/tree/main/files). Then, run the following script in [scripts](https://github.com/s-kyungyong/Sr50-AvrSr50/tree/main/scripts).
```
python find_candidate.py
```
The paths to the input files are currently hard-coded and need to be modified accordingly. This script will print out three clusters as below. 
```
3: Cluspro-43.pdb, Hdock-6.pdb, Zdock-47.pdb, Zdock-5.pdb
4: Cluspro-90.pdb, Zdock-13.pdb
7: Hdock-26.pdb, Zdock-1.pdb
```
We chose the models with the lowest numbers (e.g., Zdock-**1**.pdb) and relaxed the docking models with [amber](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/beta/relax_amber.ipynb). These models are available in Zenodo.  
