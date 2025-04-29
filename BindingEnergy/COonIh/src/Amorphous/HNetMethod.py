# NOTE
# This code is used to randomize the hydrogen positions with NeighborList function and Hydrogen bond network algorithm
# date: 2025.4.29
from ase.spacegroup import crystal
from ase.io import write
from ase.build.tools import rotation_matrix
from ase.visualize import view
import os
import numpy as np

def buildIh():
    """Build a unit cell of Ih ice.
    """
    # Define the unit cell parameters
    cell_parameters = [7.50, 7.50, 7.06, 90, 90, 120]
    
    # Define the basis atoms and their positions
    basis = [
        [0.334231, 0.334231, 0.555157],
        [0.667414, 0.667414, 0.430407],
        [0.336017, 0.336017, 0.696031],
        [0.460401, 0.460401, 0.511393],
        [0.792978, 0.669243, 0.478506]
    ]
    
    # Create the crystal structure
    Ihunitcell = crystal(
        symbols=['O', 'O', 'H', 'H', 'H'],
        basis=basis,
        spacegroup=185,
        cellpar=cell_parameters
    )
    
    return Ihunitcell

def random_Ih_hydrogen_position(crystal):
    """Randomize the positions of hydrogen atoms in the Ih ice structure.
    Args:
        crystal: The crystal structure of Ih ice.
    Returns:
        The modified crystal structure with randomized hydrogen positions. 
    """
    # Make a copy of the crystal to modify
    modified_atoms = crystal.copy()
    
    # Get all oxygen atom indices
    oxygen_indices = [i for i, atom in enumerate(modified_atoms) if atom.symbol == 'O']
    
    # Fixed O-H bond length
    oh_bond_length = 0.96  # Angstrom
    
    # Find O-O connections (potential hydrogen bond sites)
    from ase.neighborlist import NeighborList
    cutoff = 3.0  # O-O distance cutoff
    nl = NeighborList([cutoff/2] * len(modified_atoms), self_interaction=False, bothways=True)
    nl.update(modified_atoms)
    
    o_o_connections = {}
    o_o_pairs = set()
    
    for i in oxygen_indices:
        indices, offsets = nl.get_neighbors(i)
        # Filter only oxygen neighbors
        o_neighbors = [j for j in indices if modified_atoms[j].symbol == 'O']
        o_o_connections[i] = []
        
        for j, offset in zip([j for j in indices if modified_atoms[j].symbol == 'O'], 
                            [offsets[list(indices).index(j)] for j in indices if modified_atoms[j].symbol == 'O']):
            if (i, j) not in o_o_pairs and (j, i) not in o_o_pairs:
                o_o_pairs.add((i, j))
            
            # Store the connection with offset
            o_position = modified_atoms.positions[i]
            neighbor_position = modified_atoms.positions[j] + np.dot(offset, modified_atoms.get_cell())
            o_o_vector = neighbor_position - o_position
            
            # Normalize and scale to place hydrogen
            h_position_vector = o_o_vector / np.linalg.norm(o_o_vector) * oh_bond_length
            o_o_connections[i].append((j, h_position_vector, offset))
    
    # Delete all hydrogen atoms
    hydrogen_indices = [i for i, atom in enumerate(modified_atoms) if atom.symbol == 'H']
    del modified_atoms[[i for i in hydrogen_indices]]
    
    # Create new hydrogen atoms with proper positions
    # For each oxygen atom, need to place 2 hydrogen atoms
    # For each O-O connection, need to ensure only one hydrogen atom
    
    # Step 1: Randomly assign hydrogen atoms to O-O connections
    from random import shuffle
    
    # Convert pairs to list for shuffling
    o_o_list = list(o_o_pairs)
    shuffle(o_o_list)
    
    # Map to track assigned hydrogens
    o_h_count = {i: 0 for i in oxygen_indices}
    assigned_o_o = set()
    
    # First pass: place hydrogens on O-O connections
    from ase import Atom
    
    for (o1, o2) in o_o_list:
        # Skip if both oxygens already have 2 hydrogens
        if o_h_count[o1] >= 2 and o_h_count[o2] >= 2:
            continue
            
        # Randomly choose which oxygen gets the hydrogen
        if o_h_count[o1] < 2 and (o_h_count[o2] >= 2 or np.random.random() < 0.5):
            # Place hydrogen near o1
            for o_idx, h_vec, offset in o_o_connections[o1]:
                if o_idx == o2:
                    # Calculate final position
                    h_pos = modified_atoms.positions[o1] + h_vec
                    modified_atoms.append(Atom('H', h_pos))
                    o_h_count[o1] += 1
                    assigned_o_o.add((o1, o2))
                    break
        elif o_h_count[o2] < 2:
            # Place hydrogen near o2
            for o_idx, h_vec, offset in o_o_connections[o2]:
                if o_idx == o1:
                    h_pos = modified_atoms.positions[o2] + h_vec
                    modified_atoms.append(Atom('H', h_pos))
                    o_h_count[o2] += 1
                    assigned_o_o.add((o2, o1))
                    break
    
    # Second pass: add remaining hydrogens to ensure each O has 2 hydrogens
    for o_idx in oxygen_indices:
        while o_h_count[o_idx] < 2:
            # Find available O-O connections
            available_connections = []
            for neighbor, h_vec, offset in o_o_connections[o_idx]:
                if (o_idx, neighbor) not in assigned_o_o and (neighbor, o_idx) not in assigned_o_o:
                    available_connections.append((neighbor, h_vec))
            
            # If no available connections, create a random direction
            if not available_connections:
                # Generate a random unit vector
                random_vec = np.random.normal(0, 1, 3)
                random_vec = random_vec / np.linalg.norm(random_vec) * oh_bond_length
                h_pos = modified_atoms.positions[o_idx] + random_vec
            else:
                # Choose a random available connection
                neighbor, h_vec = available_connections[np.random.randint(0, len(available_connections))]
                h_pos = modified_atoms.positions[o_idx] + h_vec
                assigned_o_o.add((o_idx, neighbor))
            
            # Add hydrogen
            modified_atoms.append(Atom('H', h_pos))
            o_h_count[o_idx] += 1
    
    # Optional: Sort atoms to keep O and H atoms together
    modified_atoms = sort_atoms_by_element(modified_atoms)
    
    return modified_atoms

def sort_atoms_by_element(atoms):
    """Sort atoms in the structure by element (O first, then H)"""
    # Create a new Atoms object with sorted atoms
    from ase import Atoms
    
    # Get positions and symbols
    pos_o = [atom.position for atom in atoms if atom.symbol == 'O']
    pos_h = [atom.position for atom in atoms if atom.symbol == 'H']
    
    # Create new Atoms object with correct order
    sorted_atoms = Atoms(
        symbols=['O'] * len(pos_o) + ['H'] * len(pos_h),
        positions=pos_o + pos_h,
        cell=atoms.get_cell(),
        pbc=atoms.get_pbc()
    )
    
    return sorted_atoms

if __name__ == "__main__":
    # Build the Ih structure
    Ihunitcell = buildIh()
    
    # Output directory
    out = "/home/lihuimin/projects/BEofSpecialMoleculeOnIh/CO/out"
    
    # Save initial structure
    write(os.path.join(out, 'Ih_initial.cif'), Ihunitcell)
    
    # Randomize hydrogen positions
    randomized_Ih = random_Ih_hydrogen_position(Ihunitcell)
    # Save randomized structure
    write(os.path.join(out, 'Ih_randomized.cif'), randomized_Ih)
    
    # Visualize
    view(Ihunitcell)
    view(randomized_Ih)