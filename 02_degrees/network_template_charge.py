import sys
#Including the module from the directory one level above.
myList = [[r'../'],
          sys.path]
myList = [aPath for aPart in myList for aPath in aPart]
sys.path = myList
import sknetwork
import numpy as np
import scipy as sp
import MDAnalysis as mdan
from TrajNetworkAnalysis import *

traj_first = -1001
traj_step  = 1
traj_last  = -1

sys_name  = "##SYSNAME##"
rep_num   = "##REPSNUM##"
traj_dir  = "##RUNSDIR##"
topo_dir  = "##TOPODIR##"


traj_path = f'{traj_dir}/{sys_name}/Run_{rep_num}/{sys_name}_nvt.dcd'
topo_path = f'{topo_dir}/{sys_name}_init.data'


RDF_FirstMinima_AcidicTracts = {
    'lk_8':  {'A0':13.0, 'A1':13.0, 'A2':12.5},
    'lk_16': {'A0':14.0, 'A1':13.0, 'A2':13.0},
    'L_10':  {'A0':13.0, 'A1':13.0, 'A2':12.5},
    'L_16':  {'A0':13.0, 'A1':13.0, 'A2':12.5},
    'L_20':  {'A0':13.0, 'A1':13.0, 'A2':12.0},
    'no_A0': {'A0':9.00, 'A1':13.0, 'A2':12.0},
    'no_A1': {'A0':13.0, 'A1':9.00, 'A2':12.0},
    'no_A2': {'A0':13.0, 'A1':12.5, 'A2':9.0},
}



Peptides_Charged_Beads = {
 'lk_8': ['44', '45', '46', '47', '50', '56', '57', '58'],
 'L_10': ['46', '47', '48', '49', '60', '61', '62'],
 'L_16': ['46', '47', '48', '49', '66', '67', '68'],
 'L_20': ['46', '47', '48', '49', '70', '71', '72'],
 'no_A0': ['44', '45', '46', '47', '50', '56', '57', '58'],
 'no_A1': ['44', '45', '46', '47', '50', '56', '57', '58'],
 'no_A2': ['44', '45', '46', '47', '50', '56', '57', '58'],
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
    
    # Calculate the degrees by taking a row-sum of the adjacency matrices, and saving them
                        
    np.savez_compressed(f"{sys_name}_charged_{a_tract}", a=np.sum(adj_mats, axis=-1))

                                            
