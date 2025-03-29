from ase.build import molecule
from ase.io import write,read
from ase import Atoms
import os
import numpy as np
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
from source.BuildIhSpecial import buildIhSpecial

# version 1 def
# def find_absorption_site(slab, site_type='ontop', atom_index=0):
#     '''Find the absorption site of CO on Ih 001 surface.
    
#     Parameters:
#     ----------
#     slab: Atoms
#         The Ih(001) slab.
#     site_type: str
#         The type of absorption site, 'ontop', 'bridge', or 'hollow'.
#     atom_index: int 
#         The atoms which will be used as the absorption site.
#     Returns:
#     ----------
#     position: ndarray
#         The (x, y, z) coordinates of the absorption site.
#     '''
#     # Get the layer oxygen atoms
#     z_coords = slab.get_positions()[:, 2]
#     top_layer_z = np.max(z_coords)
#     print(f"Top layer z-coordinate: {top_layer_z}")
    
#     # Find surface oxygen atoms (within 0.5 Å of top layer)
#     surface_atoms = [atom.index for atom in slab if atom.symbol == 'O' and atom.position[2] > top_layer_z - 0.5]
    
#     # Find all atoms near the surface
#     surface_atoms_all = [atom.index for atom in slab if atom.position[2] > top_layer_z - 0.5]
    
#     print(f"Found {len(surface_atoms)} surface oxygen atoms and {len(surface_atoms_all)} total surface atoms")
    
#     # Check if we found any surface atoms
#     if not surface_atoms:
#         print("No oxygen atoms found in the surface layer. Looking for any atoms near the surface.")
#         # Fallback - use any atoms near the surface
#         if surface_atoms_all:
#             surface_atoms = surface_atoms_all
#             print(f"Using {len(surface_atoms)} atoms near the surface.")
#         else:
#             # Last resort - just use the highest atom
#             highest_atom = np.argmax(z_coords)
#             print(f"No surface atoms found. Using highest atom (index {highest_atom}, z={z_coords[highest_atom]})")
#             surface_atoms = [highest_atom]
    
#     if site_type == 'ontop':
#         # Use the atom_index parameter to select a specific atom
#         if atom_index < len(surface_atoms):
#             selected_atom = surface_atoms[atom_index]
#         else:
#             print(f"Warning: Requested atom index {atom_index} exceeds available atoms. Using first atom.")
#             selected_atom = surface_atoms[0]
            
#         site_pos = slab.get_positions()[selected_atom].copy()
#         site_pos[2] += 2.3  # Move the CO molecule 2.0 Å above the surface
#         return site_pos
    
#     elif site_type == 'bridge':
#         # Check if we have enough atoms for a bridge site
#         if len(surface_atoms) < 2:
#             print("Not enough surface atoms for a bridge site. Falling back to ontop site.")
#             site_pos = slab.get_positions()[surface_atoms[0]].copy()
#             site_pos[2] += 2
#             return site_pos
        
#         # Use the midpoint of the first two surface atoms
#         pos1 = slab.get_positions()[surface_atoms[0]]
#         pos2 = slab.get_positions()[surface_atoms[1]]
#         site_pos = 0.5 * (pos1 + pos2)
#         site_pos[2] += 1.0
#         return site_pos
    
#     elif site_type == 'hollow':
#         # Check if we have enough atoms for a hollow site
#         if len(surface_atoms) < 3:
#             print("Not enough surface atoms for a hollow site. Falling back to available site.")
#             if len(surface_atoms) == 2:
#                 # Fall back to bridge site
#                 pos1 = slab.get_positions()[surface_atoms[0]]
#                 pos2 = slab.get_positions()[surface_atoms[1]]
#                 site_pos = 0.5 * (pos1 + pos2)
#                 site_pos[2] += 1.0
#             else:
#                 # Fall back to ontop site
#                 site_pos = slab.get_positions()[surface_atoms[0]].copy()
#                 site_pos[2] += 2
#             return site_pos
        
#         # Find three nearby surface atoms
#         pos1 = slab.get_positions()[surface_atoms[0]]
#         pos2 = slab.get_positions()[surface_atoms[1]]
#         pos3 = slab.get_positions()[surface_atoms[2]]
#         site_pos = (pos1 + pos2 + pos3) / 3
#         site_pos[2] += 1.5
#         return site_pos
    
#     else:
#         raise ValueError(f'Invalid site_type: {site_type}')
    
def find_absorption_site(slab, site_type='ontop', atom_indices=None,distance=2.8):
    '''Find the absorption site of CO on Ih 001 surface.
    
    Parameters:
    slab: Atoms
        The Ih(001) slab.
    site_type: str
        The type of absorption site, 'ontop', 'bridge', or 'hollow'.
    atom_indices: int, list, or tuple
        For 'ontop': int - index of the surface atom to use
        For 'bridge': list/tuple of 2 indices - the two atoms to form the bridge
        For 'hollow': list/tuple of 3 indices - the three atoms to form the hollow site
        If None, default atoms will be used based on site_type
    distance: int
        Initial distance from molecule to absorted surface

    Returns:
    position: ndarray
        The (x, y, z) coordinates of the absorption site.
    '''
    # Get the layer oxygen atoms
    z_coords = slab.get_positions()[:, 2]
    top_layer_z = np.max(z_coords)
    print(f"Top layer z-coordinate: {top_layer_z}")
    
    # Find surface oxygen atoms (within 0.5 Å of top layer)
    surface_atoms = [atom.index for atom in slab if atom.symbol == 'O' and atom.position[2] > top_layer_z - 0.5]
    
    # Find all atoms near the surface
    surface_atoms_all = [atom.index for atom in slab if atom.position[2] > top_layer_z - 0.5]
    
    print(f"Found {len(surface_atoms)} surface oxygen atoms and {len(surface_atoms_all)} total surface atoms")
    
    # Check if we found any surface atoms
    if not surface_atoms:
        print("No oxygen atoms found in the surface layer. Looking for any atoms near the surface.")
        # Fallback - use any atoms near the surface
        if surface_atoms_all:
            surface_atoms = surface_atoms_all
            print(f"Using {len(surface_atoms)} atoms near the surface.")
        else:
            # Last resort - just use the highest atom
            highest_atom = np.argmax(z_coords)
            print(f"No surface atoms found. Using highest atom (index {highest_atom}, z={z_coords[highest_atom]})")
            surface_atoms = [highest_atom]
    
    if site_type == 'ontop':
        if atom_indices is not None:
            # Use the specified atom index
            if isinstance(atom_indices, (list, tuple)):
                # If a list was provided, use the first element
                atom_index = atom_indices[0]
            else:
                # Otherwise assume it's a single index
                atom_index = atom_indices
                
            # Verify the index is valid
            if atom_index < len(surface_atoms):
                selected_atom = surface_atoms[atom_index]
            else:
                print(f"Warning: Requested atom index {atom_index} exceeds available atoms. Using first atom.")
                selected_atom = surface_atoms[0]
        else:
            # Default to first surface atom
            selected_atom = surface_atoms[0]
            
        site_pos = slab.get_positions()[selected_atom].copy()
        site_pos[2] += distance  # Move the CO molecule distance Å above the surface
        return site_pos
    
    elif site_type == 'bridge':
        # Determine which atoms to use for the bridge
        if atom_indices is not None and isinstance(atom_indices, (list, tuple)) and len(atom_indices) >= 2:
            # Use the specified atom indices
            idx1 = atom_indices[0]
            idx2 = atom_indices[1]
            
            # Verify indices are valid
            if idx1 < len(surface_atoms) and idx2 < len(surface_atoms):
                atom1 = surface_atoms[idx1]
                atom2 = surface_atoms[idx2]
            else:
                print(f"Warning: One or more requested atom indices exceed available atoms. Using first two atoms.")
                if len(surface_atoms) < 2:
                    print("Not enough surface atoms for a bridge site. Falling back to ontop site.")
                    atom1 = surface_atoms[0]
                    site_pos = slab.get_positions()[atom1].copy()
                    site_pos[2] += 2.3
                    return site_pos
                atom1 = surface_atoms[0]
                atom2 = surface_atoms[1]
        else:
            # Default to first two surface atoms
            if len(surface_atoms) < 2:
                print("Not enough surface atoms for a bridge site. Falling back to ontop site.")
                atom1 = surface_atoms[0]
                site_pos = slab.get_positions()[atom1].copy()
                site_pos[2] += 2.3
                return site_pos
            atom1 = surface_atoms[0]
            atom2 = surface_atoms[1]
        
        # Calculate bridge position
        pos1 = slab.get_positions()[atom1]
        pos2 = slab.get_positions()[atom2]
        site_pos = 0.5 * (pos1 + pos2)
        site_pos[2] += distance  # Adjust height above surface
        return site_pos
    
    elif site_type == 'hollow':
        # Determine which atoms to use for the hollow site
        if atom_indices is not None and isinstance(atom_indices, (list, tuple)) and len(atom_indices) >= 3:
            # Use the specified atom indices
            idx1 = atom_indices[0]
            idx2 = atom_indices[1]
            idx3 = atom_indices[2]
            
            # Verify indices are valid
            if idx1 < len(surface_atoms) and idx2 < len(surface_atoms) and idx3 < len(surface_atoms):
                atom1 = surface_atoms[idx1]
                atom2 = surface_atoms[idx2]
                atom3 = surface_atoms[idx3]
            else:
                print(f"Warning: One or more requested atom indices exceed available atoms. Using default atoms.")
                # Fall back to default behavior
                if len(surface_atoms) < 3:
                    print("Not enough surface atoms for a hollow site. Falling back to available site.")
                    if len(surface_atoms) == 2:
                        # Fall back to bridge site
                        atom1 = surface_atoms[0]
                        atom2 = surface_atoms[1]
                        pos1 = slab.get_positions()[atom1]
                        pos2 = slab.get_positions()[atom2]
                        site_pos = 0.5 * (pos1 + pos2)
                        site_pos[2] += 2.0
                        return site_pos
                    else:
                        # Fall back to ontop site
                        atom1 = surface_atoms[0]
                        site_pos = slab.get_positions()[atom1].copy()
                        site_pos[2] += 2.3
                        return site_pos
                atom1 = surface_atoms[0]
                atom2 = surface_atoms[1]
                atom3 = surface_atoms[2]
        else:
            # Default to first three surface atoms
            if len(surface_atoms) < 3:
                print("Not enough surface atoms for a hollow site. Falling back to available site.")
                if len(surface_atoms) == 2:
                    # Fall back to bridge site
                    atom1 = surface_atoms[0]
                    atom2 = surface_atoms[1]
                    pos1 = slab.get_positions()[atom1]
                    pos2 = slab.get_positions()[atom2]
                    site_pos = 0.5 * (pos1 + pos2)
                    site_pos[2] += 2.0
                    return site_pos
                else:
                    # Fall back to ontop site
                    atom1 = surface_atoms[0]
                    site_pos = slab.get_positions()[atom1].copy()
                    site_pos[2] += 2.3
                    return site_pos
            atom1 = surface_atoms[0]
            atom2 = surface_atoms[1]
            atom3 = surface_atoms[2]
        
        # Calculate hollow site position (centroid of three atoms)
        pos1 = slab.get_positions()[atom1]
        pos2 = slab.get_positions()[atom2]
        pos3 = slab.get_positions()[atom3]
        site_pos = (pos1 + pos2 + pos3) / 3
        site_pos[2] += distance  # Adjust height above surface
        return site_pos
    
    else:
        raise ValueError(f'Invalid site_type: {site_type}')
def place_molecule_on_slab(slab,molecule,position):
    '''Place a molecule at a specific position on the slab
    
    Parameters:
    ----------
    slab: Atoms
        The surface slab
    molecule: Atoms
        The molecule to place on the slab
    position: array-like
        The position for the molecule's first atom
    
    Returns:
    ----------
    Combined: Atoms object
        The combined system with molecule on slab
    '''
    #Create a copy of the molecule
    molecule = molecule.copy()
    #Calculate the displacement vector
    disp = position - molecule.positions[0]
    #Move the molecule
    molecule.translate(disp)
    #Combine the slab and the molecule
    combined = slab + molecule
    # Write the combined structure
    out = "/home/lihuimin/projects/BEofSpecialMoleculeOnIh/CO/out"
    write(os.path.join(out,'combined.cif'), combined)
    return combined
