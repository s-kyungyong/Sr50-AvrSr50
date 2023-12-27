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

select helix, ss H
remove helix
select loop, ss L+T
remove loop
save avrb.pdb


For ColabDock, we also have [a cuscomized Colab notebook](), which comes from this [this notebook](https://colab.research.google.com/github/JeffSHF/ColabDock/blob/dev/ColabDock.ipynb). As AlphaLink2, the restraints can be manually defined as variables. This notebook will ask for an input structure. ColabDock.input.pdb in [files](https://github.com/s-kyungyong/Sr50-AvrSr50/tree/main/files) can be uploaded. This PDB file contains the LRR domain of Sr50 (521-956) and the mature protein of AvrSr50 with its chain ID assigned to B. We failed to run the notebook with the free T4 GPUs, and high memory GPUs will be needed for this job. Furthermore, ColabDock occasionally disort the input structures to meet the docking restraints. In this case, we remodel the pose by submitting the best model to ColabFold as a template, as outlined above. 
PyMOL>alter (chain A), resi=int(resi)+427


max_recyling_iters: 3

### Model II

Model II were generated with the following restraints for Sr50 and AvrSr50. The position for AvrSr50 is based on the mature protein.
```
K711  Q99
K711  N102
K824  E95
R904  E93
R904  E95
```


### Model III
```
++K824  D97
K711  Q99
K711  N102
K824  E95
R904  E93
R904  E95
```

### Model III
```
++R904  Q68
++K824  D97
K711  Q99
K711  N102
K824  E95
R904  Q68
R904  E93
R904  E95
```

### Model IV

```
++D641  R106
++D643  R106
K711  Q99
K711  N102
K824  E95
K824  D97
R904  Q68
R904  E93
R904  E95
```


Model II 0.372 



NB-ARC: Sr50 198-519 

KKVSIVGFGGLGKTTLARAVYDKIKGDFDCRAFVPVGQNPDMK
KVLRDILIDLGNPHSDLAILDDKQLVKKLHDFLENKRYLVIIDDIWDEMLWEGINFAFSN
RNNLGSRLITTTRNFDVSKSCCLSADDSIYKMKPLSTDDSRRLFHKRIFPDAGGCPSEFQ
QVSEDILKKCGGVPLAIITIASALASGQHVKPKHEWDILLQSLGSGVTKDNSLVEMRRIL
SFSYYNLPSHLKTCLLYLCIYPEDSTIGRDRLIWKWVAEGFVHHGDQGTSLFLVGLNYFN
QLINRSMIQPIYDELGQVHACRVHDMVLDLICNFSHEAK

mafft 7.312 fasttree 2.1.10
makeblastdb -in MLA.family.fasta -out MLA.family -dbtype 'prot'
blastp -query NBARC.fasta -db MLA.family -outfmt "6 std qlen slen" -max_hsps 1 -evalue 1e-5 -out MLA.family.against.NBARC.out
python select_sequences.py NBARC MLA.family.NBARC.only.fasta

mafft --maxiterate 1000 --globalpair --thread 56 MLA.family.NBARC.only.fasta > MLA.family.NBARC.only.msa.fasta
fasttree < MLA.family.NBARC.only.msa.fasta > MLA.family.NBARC.only.msa.nwk
root tree with Sr35


Bootstrap 0.986 42 leafs 
Sr50
0110fc85-6611-4798-8ce1-0bcd4b44830f
8b286388-55e3-41e4-8e27-0e8992ceee28
5970b758-56f6-43fb-be45-f5320835e5c7
02eb821f-594c-425f-94d7-2cc27c175112
f9987463-4fbd-4709-b6e0-05267a68370a
SECCEUnv1G0554160.1
SECCEUnv1G0554180.1
SECCEUnv1G0564510.1
SECCEUnv1G0562260.1
2da9b49d-6acc-474c-89af-ec1bf6a25868
76ffe333-cbf8-46ef-a6e8-ee5df402b9db
aebdd00c-8992-49ea-aefc-4e3bacde0641
TraesCS1B02G034500.1
0d99df4d-6b61-4ced-989b-b680ec51b2cd
10ce4a09-292d-47d1-8f6e-985d870d1a9c
bcf5aca8-cbb5-432f-8707-0206b9b5f7f6
d4926f3c-5e61-4b4c-bc57-884eb21b91f0
67370616-7575-4a71-add3-27074d15f819
4b373bb9-095b-4169-8475-6867c47e18d8
c206c1d4-fd71-4582-9bc4-6c68f79cd580
c008c820-bb9e-4897-a9f6-d5474770fe95
TraesCS1A02G028700.1
703de04f-7289-433a-96e8-a9315543db49
42aae0d9-00cd-41e3-86cb-90c055465ad5
404f3e17-f9a1-4246-90b3-b7e58a9a3c3b
43049980-d551-417b-9f05-62ed6b23ab91
01caa5dc-a2be-4bb8-ad92-7417a2e51d37
af2be739-5228-47c5-9fe3-e8ec4752f8ee
abc18fe6-3c4d-470e-b479-3e256e486a6a
8c7c74c3-2151-45ea-bebc-bde4ce29517b
87c7c421-66d8-4564-abd3-aef5d577ad54
693c55a2-02c9-42cf-b15d-729275d093ee
TraesCS1D02G028200.1
b1aade35-e5ab-49ce-ad65-65acee3c3b26
366144f0-0dcb-4df9-926e-119b150fc1c2
dd83b274-ec7c-42ea-a5ef-3dc35bf3b1b2
ca107e53-bf2b-4395-b504-59da4bd5e7ac
0effb35c-07d8-4bb3-a6fd-63b799a4d4c1
bfe83166-d286-44b3-ad6a-88c693421a68
5e4d3f1f-0a7e-403d-90d7-e5e793f27dbf
627c9d39-8f70-4d6a-b428-e52e7b542294


 python select_sequences.py LRR Sr50.family.LRR.only.fasta Sr50_clade.list
 cd-hit -c 0.98 -i Sr50.family.LRR.only.fasta -o Sr50.family.LRR.only.cd-hit.c0.98
mafft --maxiterate 1000 --globalpair --thread 56 Sr50.family.LRR.only.cd-hit.c0.98 > Sr50.family.LRR.only.cd-hit.c0.98.msa.fasta
 
 Remove TraesCS1A02G028700.1 - > almost all positiosn are gaps
