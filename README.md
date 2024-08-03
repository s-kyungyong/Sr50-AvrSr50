# Engineering the plant intracellular immune receptor Sr50 to restore recognition of the AvrSr50 escape mutant

Here, the descriptions detail the method we used to derive our models used as our structural hypotheses for [our recent preprint](). You can also access other files through [Zenodo]().

# 1. Initial structural hypothesis generation (Model I)

### Molecular docking simulations
The predicted structures of Sr50 and AvrSr50 can be found in [files](https://github.com/s-kyungyong/Sr50-AvrSr50/tree/main/files). These structures were submitted as a receptor and a ligand, respectively, to [ZDOCK](https://zdock.umassmed.edu/), [HDOCK](http://hdock.phys.hust.edu.cn/) and [ClusPro](https://cluspro.bu.edu/login.php) online servers. The default parameters were used, and no restraints were set. 

### Identification of the candidate poses
All docking models can be downloaded from [Zenodo](): 'All_docking_models.zip'. There are two more input files needed to run this step: Tomborski_2022_Sr50_group_entropy.txt and Sr50.LRR.targets.txt. These files contain entropy calculation for each residue, which comes from [this publication](https://apsjournals.apsnet.org/doi/full/10.1094/MPMI-07-22-0154-R), and the LRR residues that compise inner concave of Sr50's LRR domain based on the predicted Sr50 structure. These files are available in [files](https://github.com/s-kyungyong/Sr50-AvrSr50/tree/main/files). Then, run the following script in [scripts](https://github.com/s-kyungyong/Sr50-AvrSr50/tree/main/scripts).
```
python find_candidate.py
```
The paths to the input files are currently hard-coded and need to be modified accordingly. This script will print out three clusters as below. 
```
3: Cluspro-43.pdb, Hdock-6.pdb, Zdock-47.pdb, Zdock-5.pdb (Alternative model II)
4: Cluspro-90.pdb, Zdock-13.pdb (Alternative model I)
7: Hdock-26.pdb, Zdock-1.pdb (Model I)
```
We chose the models with the lowest numbers (e.g., Zdock-**1**.pdb from Cluster 7). To remove any abnormal features, such as steric clashes, the representative structures were remodeled with [ColabFold](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb). The following parameters were used. For this step, the multimer model (v3) also works ok. Sometimes they behave differetnly. You can play with the parameter. Make sure that num_relaxation is set to 0. Amber relaxation cannot happen when there are unknown residues (X).

```
num_relax: 0
template_mode: custom
model_type: alphafold2_ptm
```

For 'query_sequence', concatnate Sr50 and mature AvrSr50 sequences with 200 X's between them.
```
MNIVTGAMGSLIPKLGELLMDEYKLHKRIKKDVEFLKKELESMHAALIKVGEVPRDQLDRQVKLWADEVRELSYNMEDVVDKFLVRVDGDGIQQPHDNSGRFKELKNKMIGLFKKGRNHHRIADAIKEIKEQLQEVAARRDRNKVAVPNPMEPITIDPCLRALYAEATELVGIYGKRDEELMRLLSMEGDDASNKRLKKVSIVGFGGLGKTTLARAVYDKIKGDFDCRAFVPVGQNPDMKKVLRDILIDLGNPHSDLAILDDKQLVKKLHDFLENKRYLVIIDDIWDEMLWEGINFAFSNRNNLGSRLITTTRNFDVSKSCCLSADDSIYKMKPLSTDDSRRLFHKRIFPDAGGCPSEFQQVSEDILKKCGGVPLAIITIASALASGQHVKPKHEWDILLQSLGSGVTKDNSLVEMRRILSFSYYNLPSHLKTCLLYLCIYPEDSTIGRDRLIWKWVAEGFVHHGDQGTSLFLVGLNYFNQLINRSMIQPIYDELGQVHACRVHDMVLDLICNFSHEAKFVNVLDGTGNSISSQSNVRRLSLQNKMEDHQAKPLTNIMSMSRVRSITIFPPAVSIMPSLSMFEVLRVLDLSNCDLGKSSSLQLNLKGVGHLIHLRYLDLQGTQISELPTEIGNLQFLEVLDLDNNYELDELPSTLFKLRRLIYLNVMLYKVVPTPGVLQNMTSIEVLRGVLVSLNIIAQELGNLTRLRELKICFKDGNLDSYKLFVKSLGNLHHIESLSISYNSKETSFELMDLLGERWVPPVHLREFVSWMPSQLSALRGWIKRDPSHLSNLSELILWPVKEVQQEDVEIIGGLLSLRRLWIKSTHQTQRLLVIRADGFRCMMDFELNCGSAAQIMFEPGALPRAEVLVFSLGVRVAQEDGNCGFDLGLQGNLLSLRHDVFVRIYCGGARVGEAKEAEAAVRHALEAHPNHPPIDIEMTPYIAEGARDDDLCEENXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXARSLVKIDWSGSEYTILGANHYEEPNTGAAAQFPGTMTVDDGRSPYIVRKLRNSSGKRFYVFTGHPQQPIVWNPHEEIEIQFNRKFLIAVLTEFEADSQVFNHFARRQHR
```

For the custom templatse, we first made the following changes in PyMOL for each structure (Zdock-1.pdb, Zdock-5.pdb and Zdock-13.pdb) to use it as a template. Make sure that the file name of templates are composed of 4 letters (alphabets and numbers) in lower cases.
```
#properly change AvrSr50 into chain A
alter (chain B), ID = ID + 15305
alter (chain B), resi = int(resi) + 1156
alter chain B, chain="A"

#remove loops
select loops, ss L
remove loops

#save the structure
save zk1t.pdb
```
For the repesentative docking models, Zdock-1.pdb, Zdock-5.pdb and Zdock-13, the modified templates, zk1t.pdb, zk5t.pdb and zk15.pdb, are available in [files](https://github.com/s-kyungyong/Sr50-AvrSr50/tree/main/files).  

This step effectively places AvrSr50 into chain A, as for Sr50, as if Sr50 and AvrSr50 were a single monomeric unit. However, they are not a 'continuous' unit, because we introduce an index gap. Sr50 spans between 1 to 956, and AvrSr50 between 1157 to the end. This modifications force ColabFold to stick with the given structural template, introduce small structural adjustments in backbones, and re-model the flexible loops.  

This can be useful for multiple reasons. [1] if the docking models had steric clashes, AlphaFold will fix it. [2] to model AvrSr50 with AlphaFold, the template structure of AvrSr50 QCMJC has to be provided. When the template is given, however, without any other homologous sequences guiding the modeling, AlphaFold tends to duplicate the template structure to model AvrSr50. By removing the loop structures, AlphaFold will consider some flexibility within this region to model multimeric structures. [3] the 200 X's inserted between Sr50 and AvrSr50 allows separation of Sr50 and AvrSr50 into individual units, not 'joined' structures as if they were actually a monomer. 

If, as in most cases, Sr50 and AvrSr50 are separated into two chains (A and B), AlphFold use individual chains from the template to model individual proteins but will not consider their conformation to model the structure. 

Once the prediction is finished, the side chain needs to be relaxed. The modification below needs to be made first in PyMOL.
```
#reassign AvrSr50 into chain B
select effector, (resi 1156-1266)
alter (effector), chain="B"
alter (chain B), resi=int(resi)-1156
save Zdock1_c8745_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.reindexed.pdb
```
Once the model is modified, it can be relaxed with [amber](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/beta/relax_amber.ipynb). The three final structures were our Model I (Zdock1) and Alternative models I (Zodock13) and II (Zdock5). These files are available in [Zenodo](). 

## 2. Refinement of structural hypotheses

We used [ColabDock](https://github.com/JeffSHF/ColabDock) to refine our structural hypotheses and derive Model II and III. [A cuscomized Colab notebook]() was used, which mostly comes from this [this notebook](https://colab.research.google.com/github/JeffSHF/ColabDock/blob/dev/ColabDock.ipynb). This notebook will ask for an input structure. *ColabDock.input.pdb* in [files](https://github.com/s-kyungyong/Sr50-AvrSr50/tree/main/files) can be uploaded. This PDB file contains the partial NB-ARC domain (428-520), the full LRR domain of Sr50 (521-956) and the mature protein of AvrSr50 with its chain ID assigned to B. These are created from *Sr50.pdb* and *AvrSr50.pdb* in [files](https://github.com/s-kyungyong/Sr50-AvrSr50/tree/main/files). The partial NB-ARC domain acts as constaints to prevent the effector to be inserted to deep. ColabDock uses a lot of memory, so A40 will be necessary to run this job. 

### Model II and III 

Model II were generated with the following restraints for Sr50 and AvrSr50. The position for AvrSr50 is based on *the mature protein* (position within the full-length protein - 22).
```
#pairwise restraints:
K711  Q99
K711  N102
K824  E95
R904  E95

#input parameters
res_thres: 8.0

#advanced settings
The weights of each chain:
Chain A: 0.9
Chain B: 0.9
```

### Model III
```
#pairwise restraints:
D643 R106E
K824 E95
K824 D97
R904 Q68
R904 E95
R904 D97

#input parameters
res_thres: 12.0

#advanced settings
The weights of each chain:
Chain A: 0.95
Chain B: 0.95
```

The structure needs to be relaxed after the prediction with ColabDock, using the Amber relaxation outlined above. Then, re-index the structure in PyMOL so that the residue index is consistent. 
```
alter (chain A), resi=int(resi) + 427
alter (chain B), resi=int(resi) - 579
```

Unfortunately, there is no 'right' parameters for ColabDock. We observed that ColabDock has no algorithmic power to predict the correct topology of Sr50 and AvrSr50. It is simply a tool to get you the pose of Sr50 and AvrSr50 that makes sense to you and your experimental data. It allows you to have more involved conversations with your structural hypotheses. You can try different parameters and choose what makes the most sense to you. 

## 3. Final structural hypothesis (Model IV)

The last structural hypothesis is a simple refiement of Model III using the methedology outlined in section 1. Change the chain ID for AvrSr50 to 'A' with a 200 index gap and submit this to ColabFold as a template. Make sure to relax the side chains with Amber relaxation. Default parameters were fine. For this model, we used alphafold2_multimer_v3. 

```
num_relax: 0
template_mode: custom
model_type: alphafold2_multimer_v3
```


## 4. AlphaFold-deriven models of Sr50 and AvrSr50

Some AlphaFold models showed great similarity to our Model IV. Those models were predicted with ColabFold with the following parameters. For the template, we used the structure of AvrSr50 QCMJC (PDB:7MQQ) with lopp regions removed. This structure can be found in [files](https://github.com/s-kyungyong/Sr50-AvrSr50/tree/main/files) and is named as 'nolp.pdb'.
```
num_relax: 1
template_mode: custom
model_type: alphafold2_multimer_v3
num_recycles: 24
pair_mode: unpaired
```


## 5. Evolutionary analyses

The Shannon entropy can be found in [this repository](https://github.com/krasileva-group/Sr33-Sr50_analysis) and is associated with [the previous study](https://apsjournals.apsnet.org/doi/full/10.1094/MPMI-07-22-0154-R?rfr_dat=cr_pub++0pubmed&url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org).

Happy modeling! 


