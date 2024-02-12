import sys
myList = [[r'../'],
          sys.path]
myList = [aPath for aPart in myList for aPath in aPart]
sys.path = myList
import sknetwork
import numpy as np
import scipy as sp
import MDAnalysis as mdan
from TrajNetworkAnalysis import *

traj_first = 0
traj_step  = 1
traj_last  = -1

sys_name  = "lk_8"


traj_path = f'lk_8_example_traj.dcd'
topo_path = f'lk_8_init.data'


RDF_FirstMinima_AcidicTracts = {
    'rpL5':  {'A0':13.0, 'A1':13.0, 'A2':12.5},
}



Peptides_Charged_Beads = {
 'rpL5': ['44', '45', '46', '47', '50', '56', '57', '58'],
 }

N130_Tract_Charged_Beads = {
 'A0': ['6', '7', '10', '12'],
 'A1': ['22', '24', '25', '27'],
 'A2': ['33', '34', '35', '37', '39', '40', '41', '42', '43']
 }
 
for a_tract in ['A0', 'A1', 'A2']:
    analysis_func = n130_network_analysis_unweighted_adjacency_matrix_between_atom_selections

    adj_mats = analysis_func(top_file_path=topo_path, trj_file_path=traj_path,
                             at_sel_1=N130_Tract_Charged_Beads[a_tract],
                             at_sel_2=Peptides_Charged_Beads[sys_name],
                             trj_first=traj_first, trj_last=traj_last, trj_step=traj_step,
                             bond_cutoff=RDF_FirstMinima_AcidicTracts[sys_name][a_tract])
                                                
    # Removing self-edges
    for i, _ in enumerate(adj_mats):
        adj_mats[i] -= np.diag(np.diag(adj_mats[i]))
                                                
    np.savez_compressed(f"{sys_name}_charged_{a_tract}", a=np.sum(adj_mats, axis=-1))

                                            
