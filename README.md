# Sr50-AvrSr50

## 1. Initial structural hypothesis generation

### Molecular docking simulations
The predicted structures of Sr50 and AvrSr50 can be found in [files](https://github.com/s-kyungyong/Sr50-AvrSr50/tree/main/files). These structures were submitted as a receptor and a ligand, respectively, to [ZDOCK](https://zdock.umassmed.edu/), [HDOCK](http://hdock.phys.hust.edu.cn/) and [ClusPro](https://cluspro.bu.edu/login.php). The default parameters were used, and no restraints were set. 

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
We chose the models with the lowest numbers (e.g., Zdock-**1**.pdb from Cluster 7). To remove any abnormal features, such as steric clashes, the representative structures were remodeled with [ColabFold](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb). The following parameters were used. 

```
num_relax: 0
template_mode: custom
model_type: alphafold2_ptm
```

For 'query_sequence', concatnate Sr50 and mature AvrSr50 sequences with 200 X's between them.
```
MNIVTGAMGSLIPKLGELLMDEYKLHKRIKKDVEFLKKELESMHAALIKVGEVPRDQLDRQVKLWADEVRELSYNMEDVVDKFLVRVDGDGIQQPHDNSGRFKELKNKMIGLFKKGRNHHRIADAIKEIKEQLQEVAARRDRNKVAVPNPMEPITIDPCLRALYAEATELVGIYGKRDEELMRLLSMEGDDASNKRLKKVSIVGFGGLGKTTLARAVYDKIKGDFDCRAFVPVGQNPDMKKVLRDILIDLGNPHSDLAILDDKQLVKKLHDFLENKRYLVIIDDIWDEMLWEGINFAFSNRNNLGSRLITTTRNFDVSKSCCLSADDSIYKMKPLSTDDSRRLFHKRIFPDAGGCPSEFQQVSEDILKKCGGVPLAIITIASALASGQHVKPKHEWDILLQSLGSGVTKDNSLVEMRRILSFSYYNLPSHLKTCLLYLCIYPEDSTIGRDRLIWKWVAEGFVHHGDQGTSLFLVGLNYFNQLINRSMIQPIYDELGQVHACRVHDMVLDLICNFSHEAKFVNVLDGTGNSISSQSNVRRLSLQNKMEDHQAKPLTNIMSMSRVRSITIFPPAVSIMPSLSMFEVLRVLDLSNCDLGKSSSLQLNLKGVGHLIHLRYLDLQGTQISELPTEIGNLQFLEVLDLDNNYELDELPSTLFKLRRLIYLNVMLYKVVPTPGVLQNMTSIEVLRGVLVSLNIIAQELGNLTRLRELKICFKDGNLDSYKLFVKSLGNLHHIESLSISYNSKETSFELMDLLGERWVPPVHLREFVSWMPSQLSALRGWIKRDPSHLSNLSELILWPVKEVQQEDVEIIGGLLSLRRLWIKSTHQTQRLLVIRADGFRCMMDFELNCGSAAQIMFEPGALPRAEVLVFSLGVRVAQEDGNCGFDLGLQGNLLSLRHDVFVRIYCGGARVGEAKEAEAAVRHALEAHPNHPPIDIEMTPYIAEGARDDDLCEENXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXARSLVKIDWSGSEYTILGANHYEEPNTGAAAQFPGTMTVDDGRSPYIVRKLRNSSGKRFYVFTGHPQQPIVWNPHEEIEIQFNRKFLIAVLTEFEADSQVFNHFARRQHR
```

For the custom template, upload each docking model. We first made the following changes in PyMOL to use it as a template. 
```
alter (chain B), ID = ID + 15305
alter (chain B), resi = int(resi) + 1156
alter chain B, chain="A"
save zk1t.pdb
```
For the repesentative docking models, Zdock-1.pdb, Zdock-5.pdb and Zdock-13, the modified templates, zk1t.pdb, zk5t.pdb and zk15.pdb, are available in [files](https://github.com/s-kyungyong/Sr50-AvrSr50/tree/main/files). 

After the prediction, the best model is then re-indexed in PyMOL as below.
```
select effector, (resi 1156-1266)
alter (effector), chain="B"
alter (chain B), resi=int(resi)-1156
save Zdock1_c8745_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.reindexed.pdb
```
Once the model is modified, it can be relaxed with [amber](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/beta/relax_amber.ipynb). Acess these prediction files through Zenodo().  

## 2. Refinement of structural hypotheses

We used [AlphaLink2](https://github.com/Rappsilber-Laboratory/AlphaLink2) and [ColabDock](https://github.com/JeffSHF/ColabDock) to refine our structural hypotheses. 

For AlphaLink2, we used [our customized Colab notebook](), which has a slight modification from [the original one](https://colab.research.google.com/github/Rappsilber-Laboratory/AlphaLink2/blob/main/notebooks/alphalink2.ipynb#scrollTo=z4QQBdHdv4yK). Here, all that need to be done is to define restraints. 

For ColabDock, we also have [a cuscomized Colab notebook](), which comes from this [this notebook](https://colab.research.google.com/github/JeffSHF/ColabDock/blob/dev/ColabDock.ipynb). As AlphaLink2, the restraints can be manually defined as variables. This notebook will ask for an input structure. ColabDock.input.pdb in [files](https://github.com/s-kyungyong/Sr50-AvrSr50/tree/main/files) can be uploaded. This PDB file contains the LRR domain of Sr50 (521-956) and the mature protein of AvrSr50 with its chain ID assigned to B. We failed to run the notebook with the free T4 GPUs, and high memory GPUs will be needed for this job. 


### Model II

Model II were generated with the following restraints for Sr50 and AvrSr50. The position for AvrSr50 is based on the mature protein.
```
R688  Q99
K711  Q99
K824  E95
R904  E95
```


### Model III
```
R904  Q68

R688  Q99
K711  Q99
K824  E95
R904  E95
```


```

Model II 0.372 
