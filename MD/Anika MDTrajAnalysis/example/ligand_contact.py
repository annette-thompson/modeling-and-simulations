#!/ usr / bin / env python
import math
import mdtraj as md
import numpy as np
import argparse
import sys
import os.path
from itertools import product

#Declare arguments
parser = argparse.ArgumentParser(description = 'Determination of Ligand Residue Contacts')
parser.add_argument('-t', required=True, help='File name for input trajectory')
parser.add_argument('-g', required=True, help= 'File name for input topology (gro format)')
parser.add_argument('-m', required=False, type=int, default = 0, help= 'Supply the number of missing terminal residues(default 0)')
parser.add_argument('-l', required=True, type = int, help= 'Ligand residue number')
parser.add_argument('-ln', required=False, default = 'LIG', type = str, help= 'Ligand name in GRO file')
parser.add_argument('-s', required=False, default = 'none', help= 'File containing sections of interest(txt)')

#Import Arguments
args = parser.parse_args()
File_traj = args.t
File_gro = args.g
miss_res = args.m
lig = args.l
lig_name = args.ln
sect = args.s
if sect.split('.')[-1] != 'txt' and sect != 'none': #Add default file extension if not in input
    sect = sect + '.txt'

#Source custom functions
current_directory = os.path.dirname(os.path.realpath(__file__))
prefix = current_directory.rsplit('/',1)[0]
sys.path.insert(1, prefix + '/Traj_process/')
import load_data 

#Load Trajectory
traj_ns = load_data.mdtraj_load(File_traj, File_gro, True, True)
top = traj_ns.topology

#Determine indices of protein and ligand residues in topology
lig_res = load_data.lig_check(lig, miss_res, traj_ns, lig_name)
prot_res = np.linspace(0, lig_res-1, num = lig_res)

#Determine distance between all protein residues and ligand
res_pairs = list(product([lig_res], prot_res))
[dist, pairs] = md.compute_contacts(traj_ns, contacts=res_pairs, scheme='closest-heavy', ignore_nonprotein = False, periodic=True, soft_min = False)

#Output % contact for each residue
if not os.path.exists('lig_inter'):
    os.mkdir('lig_inter')
output_per = open('lig_inter/residue_lig_per.txt', 'w')
output_high_contact = open('lig_inter/residue_lig_high_contact.txt', 'w')

#Compute % contact
per_contact = np.zeros(len(prot_res))
for i in range(0, len(prot_res)):
    dist_i = dist[:,i]
    
    contact = 0 #conter for protein--ligand contact
    for j in dist_i:
        if j < 0.4:
            contact += 1
    per_contact[i] = 100 * (contact/len(dist_i))
    output_per.write(str(prot_res[i]+1+miss_res) + ': ' + str(per_contact[i])+ '\n')
    if per_contact[i] > 75:
        output_high_contact.write(str(prot_res[i]+1+miss_res) + '\n')

if sect != 'none':
    #Compute the number of simultaneous contacts per section
    sections = open(sect, 'r').readlines()#Load sections of interest file

    #Loop through sections of interest
    for i in range(len(sections)):
        name1, name2, sect1, sect2 = load_data.read_sections(sections, i, miss_res, 1, 2)
        res_pairs = list(product([lig_res], sect1))

        [dist, pairs] = md.compute_contacts(traj_ns, contacts=res_pairs, scheme='closest-heavy', ignore_nonprotein = False, periodic=True, soft_min = False)
    
        #Determine simultaenrous contacts
        frames, pairs = np.shape(dist)
        simul_contact = np.zeros(frames)
        for t in range(frames):
            for r in range(pairs):
                if dist[t][r] < 0.4:
                    simul_contact[t] += 1
        np.savetxt('lig_inter/' + name1 + '_lig_contacts.txt', simul_contact)

print('Ligand Interaction Analysis Complete')
