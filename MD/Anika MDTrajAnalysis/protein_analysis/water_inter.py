def traj_sect(traj, prot_res, offset):
    import mdtraj as md
    res = prot_res - offset
    top = traj.topology
    traj_sect = top.select('resid ' + str(res)) #Select only atoms in the given section
    return traj_sect

def water_neighbor(traj, res, offset, res_not_water, threshold=0.4, atom_or_res='res'):
    import mdtraj as md
    import numpy as np

    #Seperate protein residues of interest
    if atom_or_res == 'res':
        res_ind = traj_sect(traj, res, offset)
    else:
        res_ind = res

    #Compute neighboring atoms for all residues of interest
    neighbors = md.compute_neighbors(traj, threshold, res_ind, haystack_indices=None, periodic=False)
    
    #Determine which neighbors are water molecules
    top = traj.topology
    traj_solv = top.select('water')

    water_neighbors = []
    for i in range(len(neighbors)):
        water_neighbors_i = []
        neighbors_i = neighbors[i]
        for j in neighbors_i:
            if j in traj_solv:
                water_neighbors_i.append(j)
        #Add array of water neighbors to full list
        water_neighbors.append(water_neighbors_i)
    #Delete neighbors array as no longer needed
    del neighbors

    #Determine residue from atom number
    water_res_neighbors = []
    traj_OW = top.select('water and name O')
    traj_HW1 = top.select('water and name H1')
    traj_HW2 = top.select('water and name H2')
    traj_not_water = top.select('protein or resname LIG')
    for i in range(len(water_neighbors)):
        water_neighbors_i = water_neighbors[i]
        water_res_neighbors_i = []
        for j in water_neighbors_i:
            if j in traj_OW:
                x = 1
            elif j in traj_HW1:
                x = 2
            elif j in traj_HW2:
                x = 3
            res = np.round((j - traj_not_water[-1] - x)/3) + res_not_water + 1
            if res not in water_res_neighbors_i: #remove duplicate residues
                water_res_neighbors_i.append(res)
        water_res_neighbors.append(water_res_neighbors_i)
    return water_res_neighbors

def seperate():
    if lig_res_water_neighbor != False:
        water_res_lig_prot = []
        for i in range(len(water_res_neighbors)): #Loop through frames
            prot_res_water = water_res_neighbors[i]
            lig_res_water = lig_res_water_neighbor[i]
            check = 0
            for j in prot_res_water:
                if j in lig_res_water:
                    check = 1
                    break
            water_res_lig_prot.append(check)

    if lig_res_water_neighbor == False:
        return water_res, mean_neighbors, water_res_neighbors
    else:
        return water_res, mean_neighbors, water_res_lig_prot

def prot_water_dist(traj, prot_res, offset, water_indices):
    from itertools import product
    import mdtraj as md
    import numpy as np

    #Seperate protein residues of interest
    prot_ind = traj_sect(traj, prot_res, offset)

    #Determine all solvent protein residue pairs
    prot_solv = list(product(prot_ind, water_indices))
    
    #Don't compute distances if no water molecules selected
    if not water_indices:
        dist = 0
    else:
        #Determine distance b/w all protein residues and all solvent molecules
        dist = md.compute_distances(traj, prot_solv, periodic=True)

    return dist
