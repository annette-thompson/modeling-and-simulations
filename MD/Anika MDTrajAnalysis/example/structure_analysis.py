#!/ usr / bin / env python
import mdtraj as md
import numpy as np
import argparse
import sys
import os.path
import warnings
import pandas as pd

#Silence MDTraj Warnings
warnings.filterwarnings("ignore")

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of RMSD for BB Atoms in proteins and/or ligands')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-e', required=False, default=0, type=int, help= 'If input trajectory is equilibrated input equilibration time otherwise input 0(ns)')
parser.add_argument('-m', required=False, type=int, default = 0, help= 'Supply the number of missing terminal residues(default 0)')
parser.add_argument('-ref', required=False, type=str, default = 'none', help='Reference structure for RMSD')
parser.add_argument('-mref', required=False, type=int, default = 0, help= 'Supply the number of missing terminal residues for the reference structure(default 0)')
parser.add_argument('-ref2', required=False, type=str, default = 'none', help='Reference structure for RMSD if two proteins present')
parser.add_argument('-aa', required=False, default=0, type=int, help= 'Number of atoms in receptor')
parser.add_argument('-aa_ref', required=False, default=0, type=int, help= 'Number of atoms in protein 1 for reference')
parser.add_argument('-l', required=False, default = 'none', type = str, help= 'Ligand name if ligand analysis should be performed')
parser.add_argument('-lref', required=False, type=str, default = 'none', help='Reference structure for Ligand COM RMSD')
parser.add_argument('-s', required=False, default = 'none', type = str, help= 'File name for residue ranges for computed RMSD (sect_name ref_res_initial ref_res_final traj_res_initial traj_res_final)')
parser.add_argument('-rn', required=False, default = 'self', help='Reference name for RMSD')
parser.add_argument('-ln', required=False, default = 'self', help='Reference name for Ligand RMSD')

#Import Arguments
args = parser.parse_args()
File_traj = args.t
File_gro = args.g
eq_time = args.e
miss_res = args.m
ref = args.ref
ref2 = args.ref2
miss_ref = args.mref
aa_atom = args.aa
aa_atom_ref = args.aa_ref
lig = args.l
lig_ref = args.lref
lig_name = args.ln
rmsd_sect = args.s
ref_name = args.rn


#Source custom functions
current_directory = os.path.dirname(os.path.realpath(__file__))
prefix = current_directory.rsplit('/',1)[0]
sys.path.insert(1, prefix + '/Traj_process/')
import load_data 

sys.path.insert(1, prefix + '/protein_inter/')
import process_traj

sys.path.insert(1, prefix + '/ligand_analysis/')
import lig_motion

#Load Trajectory
if ref_name == 'self':
    rm_uncorr = False
else:
    rm_uncorr = True
    
traj = load_data.mdtraj_load(File_traj, File_gro, True, rm_uncorr)
top = traj.topology

#Set protein offset based on missing residues
offset = 1 + miss_res
offset_ref = 1 + miss_ref

#Add label to files if the trajectory is not equilibrated
if eq_time != 0:
    name_add = ''
else:
    name_add = '_unequil'

#Check if output should go to seperate directory
if not os.path.exists('./rmsd/'):
    os.mkdir('rmsd')

#Load reference
if ref != 'none':
    if aa_atom != 0:
        if aa_atom_ref == 0:
            aa_atom_ref = aa_atom
        print(aa_atom_ref)
        ref_bb = load_data.load_ref(ref, 'backbone and index <= ' + str(aa_atom_ref))
        traj_bb = traj.atom_slice(top.select('backbone and index <= ' + str(aa_atom))) #Backbond atoms of protein 1 only
        ref2_bb = load_data.load_ref(ref2, 'backbone and index > ' + str(aa_atom_ref))
        traj2_bb = traj.atom_slice(top.select('backbone and index > ' + str(aa_atom))) #Backbond atoms of protein 1 only
    elif aa_atom == 0 and aa_atom_ref != 0:
        ref_bb = load_data.load_ref(ref, 'backbone and index <= ' + str(aa_atom_ref))
        traj_bb = traj.atom_slice(top.select('backbone')) #Backbond atoms of protein only
    else:
        ref_bb = load_data.load_ref(ref, 'backbone')
        traj_bb = traj.atom_slice(top.select('backbone')) #Backbond atoms of protein only

#Compute full BB RMSD
if ref != 'none' and ref_name == 'self':
    #Calculate RMSF from reference structure
    rmsf_data = md.rmsf(traj_bb, ref_bb, parallel=True, precentered=False)

    #Calculate RMSD for full protein relative to reference structure
    rmsd_full_uncorr, t_full = process_traj.compute_rmsd(traj_bb, ref_bb)
    
    if aa_atom  != 0:
        #Calculate RMSF from reference structure
        rmsf_data2 = md.rmsf(traj2_bb, ref2_bb, parallel=True, precentered=False)

        #Calculate RMSD for full protein relative to reference structure
        rmsd_full_uncorr2 = process_traj.compute_rmsd(traj2_bb, ref2_bb, t_full)
    
    #Save RMSD and RMSF values to CSV
    if aa_atom != 0 and ref2 != 'none':
        df = pd.DataFrame({'RMSD': rmsd_full_uncorr, 'Prot2 RMSD': rmsd_full_uncorr2})
        if len(rmsf_data) > len(rmsf_data2):
            rmsf_data2 = np.pad(rmsf_data2, (0,len(rmsf_data) - len(rmsf_data2)), 'constant')
        else:
            rmsf_data =  np.pad(rmsf_data, (0,len(rmsf_data2) - len(rmsf_data)), 'constant')
        df2 = pd.DataFrame({'RMSF': rmsf_data2, 'Prot2 RMSF': rmsf_data2})
    else:
        df = pd.DataFrame({'RMSD': rmsd_full_uncorr}) 
        df2 = pd.DataFrame({'RMSF': rmsf_data})
    df.to_csv('rmsd/rmsd_ref_' + str(ref_name) + name_add + '.csv')
    df2.to_csv('rmsf_ref_' + str(ref_name) + name_add + '.csv')
    
    np.savetxt('uncorrelated_frames' + name_add + '.txt', t_full) #Save indices for uncorrelated frames to file

    print('Full BB RMSD Computed')

    #Delete unneeded arrays to save memory
    del rmsf_data; del rmsd_full_uncorr
else:
    print('Full BB RMSD Skipped')

#Compute BB RMSD by sections
if ref != 'none' and rmsd_sect != 'none':
    #Initiate Data Frame
    df = pd.DataFrame()
    #Loop through sections from input file
    sections = open(rmsd_sect, 'r').readlines()
    for i in sections:
        #Split line by spaces
        line = i.split(' ')
        if len(line) == 5:
            [name, ref_res_initial, ref_res_final, traj_res_initial, traj_res_final] = line
        elif len(line) == 3:
            [name, traj_res_initial, traj_res_final] = line
            ref_res_initial = traj_res_initial
            ref_res_final = traj_res_final
        else:
            raise Exception('Error in input file! Can only accept 3 or 5 elements per line.')
            exit()
            
        #Seperate sections
        ref_res = [int(ref_res_initial) - offset_ref, int(ref_res_final) - offset_ref]
        traj_res = [int(traj_res_initial) - offset, int(traj_res_final) - offset]

        #Determine if section in protein 1 or 2
        if traj_res[0] < traj_bb.n_residues:
             ref_sect = ref_bb.atom_slice(ref_bb.topology.select(str(ref_res[0]) + ' <= resid and resid <= ' + str(ref_res[1]))) #Limit trajectory to the section of choice
             traj_sect = traj_bb.atom_slice(traj_bb.topology.select(str(traj_res[0]) + ' <= resid and resid <= ' + str(traj_res[1]))) #Limit trajectory to the section of choice
        else:
             ref_sect = ref2_bb.atom_slice(ref2_bb.topology.select(str(ref_res[0]-traj_bb.n_residues) + ' <= resid and resid <= ' + str(ref_res[1]-traj_bb.n_residues))) #Limit trajectory to the section of choice
             traj_sect = traj2_bb.atom_slice(traj2_bb.topology.select(str(traj_res[0]-traj_bb.n_residues) + ' <= resid and resid <= ' + str(traj_res[1]-traj_bb.n_residues))) #Limit trajectory to the section of choice
        
        if traj_sect.n_atoms != ref_sect.n_atoms:
             print(f'Skipping Section {name} because different atoms selected')
        else:
            #Compute RMSD for section of interest
            if ref_name == 'self':
                rmsd_sect_uncorr = process_traj.compute_rmsd(traj_sect, ref_sect, t_full)
            else:
                rmsd_sect_uncorr = process_traj.compute_rmsd(traj_sect, ref_sect, 'N/A')
        
            #Save RMSD to file 
            df['RMSD' + name] =  rmsd_sect_uncorr
    df.to_csv('rmsd/rmsd_sections_ref_' + str(ref_name) + '.csv')

    print('Section BB RMSD Completed')
else:
    print('Section BB RMSD Skipped')

#Compute Ligand COM RMSD
if lig != 'none':
    #Load reference
    if aa_atom != 0:
        ref_ns = load_data.load_ref(lig_ref, '(backbone and index <= ' + str(aa_atom_ref) + ') or resname ' + lig)
        traj_ns = traj.atom_slice(traj.topology.select('(backbone and index <= ' + str(aa_atom) + ') or resname ' + lig))
    elif aa_atom_ref != 0:
        ref_ns = load_data.load_ref(lig_ref, '(backbone and index <= ' + str(aa_atom_ref) + ') or resname ' + lig)
        traj_ns = traj.atom_slice(traj.topology.select('backbone or resname ' + lig))
    else:
        ref_ns = load_data.load_ref(lig_ref, 'backbone or resname ' + lig)
        traj_ns = traj.atom_slice(traj.topology.select('backbone or resname ' + lig))

    if ref_ns.n_atoms != traj_ns.n_atoms:
        raise Exception(f'Reference and Trajectory sections do not match for ligand alignment')

    #Remove uncorrelated frames
    traj_uncorr = load_data.remove_uncorr('uncorrelated_frames.txt', traj_ns)
    
    #Comput Ligand COM displacment and RMSD
    displacment, lig_rmsd = lig_motion.com_rmsd(ref_ns, traj_uncorr, lig)
    
    #Output COM RMSD
    output = open('lig_com_rmsd_' + lig_name + '.txt', 'w')
    output.write(str(lig_rmsd))
    
    #Output displacment for bootstrapping
    np.savetxt('lig_com_dis_' + lig_name + '.txt', displacment)

    print('Ligand COM RMSD Completed')
    #Delete unneeded arays
    del displacment

else:
    print('Ligand COM RMSD Skipped')

#Compute Ligand Heavy atom RMSD
if lig != 'none':
    #Load reference
    ref_bb = load_data.load_ref(lig_ref, 'backbone or resname ' + lig)
    
    #Seperate Ligand heavy atoms
    ref_sect = ref_bb.atom_slice(ref_bb.topology.select('resname ' + lig)) #Limit trajectory to the section of choice
    traj_sect = traj_uncorr.atom_slice(traj_uncorr.topology.select('resname ' + lig)) #Limit trajectory to the section of choice

    #Compute RMSD for section of interest
    rmsd_sect_uncorr = md.rmsd(traj_sect, ref_sect, parallel=True, precentered=False)

    #Save RMSD to file
    np.savetxt('rmsd_lig_heavy_ref_' + lig_name + '.txt', rmsd_sect_uncorr)

    print('Ligand Heavy Atom RMSD Completed')

else:
    print('Ligand Heavy Atom RMSD Skipped')
print('-------------------------------------------------------------')
