# -*- coding: utf-8 -*-

from Bio.PDB import PDBParser, NeighborSearch
import warnings
from Bio import BiopythonWarning
import sys
import os
import numpy as np
from prody import parsePDB, confProDy, calcRMSD

warnings.simplefilter('ignore', BiopythonWarning)
#confProDy(verbosity="none")

def get_backbone_atoms(PDBs):
 
    ca_atoms = {}
    for p in PDBs:
        structure = parsePDB(p)
        selection = structure.select("(chid B) and (name CA)")
        ca_atoms[p] = selection    
            
    return ca_atoms

def clusterize(PDBs, rmsd_cutoff):

    clusters_found = 0
    clusters = {clusters_found: [PDBs[0]]}

    # Read all structures backbone atoms
    backbone_atoms = get_backbone_atoms(PDBs)

    for j in PDBs[1:]:
        in_cluster = False
        for cluster_id in list(clusters.keys()):
            # For each cluster representative
            representative_id = clusters[cluster_id][0]
            
            rmsd = calcRMSD(backbone_atoms[representative_id], backbone_atoms[j]).round(4)

            if rmsd <= rmsd_cutoff:
                clusters[cluster_id].append(j)
                in_cluster = True
                break

        if not in_cluster:
            clusters_found += 1
            clusters[clusters_found] = [j]
    return clusters


def identify_hvLRR(entropy_list, target_LRR_list):
        
    initial_positions = []  
    final_positions   = []
    with open(target_LRR_list) as f:
        for line in f:
            start, end = map(int, line.split()[1].split("-"))
            initial_positions.extend(range(start, end + 1))
    
    with open(entropy_list, 'r') as f:
        for line in f:
            if float(line.split()[-1]) >= 0.347:
                pos = int(line.split()[0])
                
                if pos in initial_positions: final_positions.append(pos)
    
    return final_positions
                

def check_requirement(pdb_file, cutoff, cutoff2, hv_target_LRR):

  # Parse the PDB file
  parser = PDBParser()
  structure = parser.get_structure("PDB", pdb_file)

  # Get the chains A and B from the structure
  receptor_chain = structure[0]['A'] # Receptor
  effector_chain = structure[0]['B'] # Effector

  # Create a NeighborSearch object for the atoms in receptor chain
  ns = NeighborSearch(list(effector_chain.get_atoms()))


  # Set some parameters
  e494 = False # Check if E494 touches the effector
  touch = False # Check if effector touches CC and NB-ARC
  hv_positions =  list() # Collect hvLRR targeting effector
  q99 = False # Check if Q99 touches LRR
  
  # Iterate over the residues in chain A
  for receptor_residue in receptor_chain:
      
      # Get residue number
      receptor_pos = receptor_residue.get_id()[1]
      receptor_atoms = list(receptor_residue.get_atoms())
      receptor_aminoacid = receptor_residue.resname

      # Iterate over the atoms in the residue
      # Calculate the distance from the backbone of AvrSr35 (alpha carbon & beta carbon)
      for receptor_atom in receptor_atoms:
                if receptor_atom.name == 'CB' or (receptor_atom.name == 'CA' and receptor_aminoacid == 'GLY'):
                      try:
                          for atom in ns.search(receptor_atom.coord, level="A", radius=cutoff):
                            effector_aminoacid = atom.get_parent().resname
                            effector_pos = atom.get_parent().get_id()[1]
                            if atom.name == 'CB' or (atom.name == 'CA' and effector_aminoacid == 'GLY'):
                                # Checi if 494E contacts the effector
                                if receptor_pos == 494:
                                    e494 = True
                                
                                if ( 1 <= receptor_pos <= 520) and not ( 492 <= receptor_pos <= 499):
                                    touch = True
                                elif receptor_pos > 520:
                                    if effector_pos == 99:
                                        q99 = True
                                    if receptor_pos in hv_target_LRR:
                                        hv_positions.append(receptor_pos)
                                    
                                
                                    
                      except ValueError:
                          pass
                      
                      if receptor_pos == 494:
                          try:
                              for atom in ns.search(receptor_atom.coord, level="A", radius=cutoff2):
                                effector_aminoacid = atom.get_parent().resname
                                effector_pos = atom.get_parent().get_id()[1]
                                if atom.name == 'CB' or (atom.name == 'CA' and effector_aminoacid == 'GLY'):
                                    # Checi if 494E contacts the effector
                                    e494 = True
                                
                                         
                          except ValueError:
                              pass
                      

  return e494, touch, list(set(hv_positions)), q99

# Inputs 
input_dir = '../03.Structural_analysis/Sr50/' # folder that contains the two inputs below
entropy = f'{input_dir}/Tomborski_2022_Sr50_group_entropy.txt'
target_LRR = f'{input_dir}/Sr50.LRR.targets.txt'

# Distance cutoffs
cutoff = 8.0
cutoff2= 12.0

# rmsd cutoffs
rmsd_cutoff = 3.0

#Collect hv LRR positions 
hv_target_LRR_list = identify_hvLRR(entropy, target_LRR)

# # For all PDBs check the requirment
pdb_dir = './All' # folder that contains ZDOCK, HDOCK and ClusPro models
PDBs = [ f'{pdb_dir}/{p}' for p in os.listdir(pdb_dir) if p.endswith('pdb')]

models = [] 
for i, pdb in enumerate(PDBs):
    e494, touch, hv_contact, q99  = check_requirement(pdb, cutoff, cutoff2, hv_target_LRR_list)
    print(f'{i} {pdb.split("/")[-1]}: E494 - {e494}; Q99 - {q99}; Touch CC/NB - {touch}; hvLRR - { len(hv_contact) }')
    if e494 == True and q99 == True and touch == False and len(hv_contact) >= 12:
        models.append(pdb)

# Cluster the plausible models
clusters= clusterize(models, rmsd_cutoff)

for clust in clusters:
    members = [ m.replace(pdb_dir +"/", '') for m in clusters[clust]]
    prefix =  list(set([ p.split("-")[0] for p in members]))
    
    if len(prefix) > 1:
        print(f'{clust}: {", ".join(members)}')
