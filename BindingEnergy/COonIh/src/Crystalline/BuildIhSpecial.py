# Version: 2025.3.28
# from ase.spacegroup import crystal
# from ase.build import surface, cut
# from ase.io import write, read
# from ase.build import molecule
# import os

# def buildIhSpecial(indices=(0, 0, 1), plane=(1, 1, 1)):
#     '''functions: Build ice Ih slab with specified Miller indices.
#     The dataset comes from https://next-gen.materialsproject.org/materials/mp-703459?formula=H2O.
#     After finishing building the slab, the structure is written to a CIF file in your out directory.
    
#     Parameters:
#     ----------
#     indices : tuple
#         Miller indices defining the surface orientation, default (0,0,1)
#     plane : tuple
#         Repetition of the surface unit cell in x, y, z directions, default (1,1,1)
    
#     Returns:
#     ----------
#     Atoms object
#         The constructed ice Ih slab with specified orientation
#     '''
#     # Basis positions
#     basis = [
#         [0.334231, 0.334231, 0.555157],
#         [0.667414, 0.667414, 0.430407],
#         [0.336017, 0.336017, 0.696031],
#         [0.460401, 0.460401, 0.511393],
#         [0.792978, 0.669243, 0.478506],
#     ]
#     # Lattice parameters
#     a = 7.50
#     b = 7.50
#     c = 7.06
#     alpha = 90
#     beta = 90
#     gamma = 120
#     # Lattice
#     cellpar = [a, b, c, alpha, beta, gamma]
#     # Create the crystal
#     Ih = crystal(symbols=['O', 'O', 'H', 'H', 'H'],
#                 basis=basis, spacegroup=185, cellpar=cellpar)
#     Ih.translate([0, 0, 1])
#     Ih.wrap()
    
#     # Extract Miller indices to create surface name
#     h, k, l = indices
#     surface_name = f"Ih{h}{k}{l}"
    
#     # Build the slab with specified Miller indices
#     n_bilayers = 1
#     layers = n_bilayers//2 + n_bilayers%2
#     slab = surface(Ih, indices, layers=layers, vacuum=10)
    
#     if n_bilayers % 2 == 1:
#         # remove extra layers
#         extra_indices = [atom.index for atom in slab if atom.z < 7]
#         del slab[extra_indices]
    
#     # Center the slab
#     slab.center(vacuum=10, axis=2)
    
#     orthogonal = False
#     if orthogonal:
#         slab = cut(slab, (2, 1, 0), (0, 1, 0), (0, 0, 1))
#         slab.rotate(30, '-z', rotate_cell=True)
#         slab.translate([1, 0, 0])
#         slab.wrap()
#     slab = slab.repeat(plane)
#     # Write the structure with appropriate name based on Miller indices
#     out_dir = "/home/lihuimin/projects/BEofSpecialMoleculeOnIh/CO/out"
#     file_name = f"{surface_name}.cif"
#     write(os.path.join(out_dir, file_name), slab)
    
#     return slab



from ase.spacegroup import crystal
from ase.build import surface, cut
from ase.io import write, read
from ase.build import molecule
import os
import numpy as np

def buildIhSpecial(indices=(0, 0, 1), plane=(1, 1, 1)):
    '''functions: Build ice Ih slab with specified Miller indices.
    The dataset comes from https://next-gen.materialsproject.org/materials/mp-703459?formula=H2O.
    After finishing building the slab, the structure is written to a CIF file in your out directory.
    
    Parameters:
    ----------
    indices : tuple
        Miller indices defining the surface orientation, default (0,0,1)
    plane : tuple
        Repetition of the surface unit cell in x, y, z directions, default (1,1,1)
    
    Returns:
    ----------
    Atoms object
        The constructed ice Ih slab with specified orientation
    '''
    # Basis positions
    basis = [
        [0.334231, 0.334231, 0.555157],
        [0.667414, 0.667414, 0.430407],
        [0.336017, 0.336017, 0.696031],
        [0.460401, 0.460401, 0.511393],
        [0.792978, 0.669243, 0.478506],
    ]
    # Lattice parameters
    a = 7.50
    b = 7.50
    c = 7.06
    alpha = 90
    beta = 90
    gamma = 120
    # Lattice
    cellpar = [a, b, c, alpha, beta, gamma]
    # Create the crystal
    Ih = crystal(symbols=['O', 'O', 'H', 'H', 'H'],
                basis=basis, spacegroup=185, cellpar=cellpar)
    Ih.translate([0, 0, 1])
    Ih.wrap()
    
    # Extract Miller indices to create surface name
    h, k, l = indices
    surface_name = f"Ih{h}{k}{l}"
    
    # Build the slab with specified Miller indices - use more layers for better surfaces
    n_layers = 1  # Use 1 complete bilayers
    slab = surface(Ih, indices, layers=n_layers, vacuum=10.0) # Increased vacuum to 10 Å
    
    # Repeat the slab in x,y directions only - avoids z-stacking issues
    if len(plane) == 3:
        slab = slab.repeat((plane[0], plane[1], plane[2]))
    else:
        slab = slab.repeat(plane)
    
    # Remove bottom layers - use percentage instead of hard-coded value
    # Find z positions and use a cutoff based on overall cell height
    z_pos = slab.positions[:,2]
    z_min, z_max = np.min(z_pos), np.max(z_pos)
    cell_height = slab.cell[2,2]
    
    # Keep top 60% of the slab (removes ~40% from bottom)
    z_cutoff = z_min + (z_max - z_min) * 0.4
    
    # Remove atoms below the cutoff
    indices_to_remove = [atom.index for atom in slab if atom.position[2] < z_cutoff]
    if indices_to_remove:
        slab = slab[list(set(range(len(slab))) - set(indices_to_remove))]
    
    # Center the slab with more vacuum (15 Å)
    slab.center(vacuum=15, axis=2)
    
    # Make the cell orthogonal if requested
    orthogonal = False
    if orthogonal:
        slab = cut(slab, (2, 1, 0), (0, 1, 0), (0, 0, 1))
        slab.rotate(30, '-z', rotate_cell=True)
        slab.translate([1, 0, 0])
        slab.wrap()
    
    # Write the structure with appropriate name based on Miller indices
    out_dir = "/home/lihuimin/projects/BEofSpecialMoleculeOnIh/CO/out"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    file_name = f"{surface_name}_{plane[0]}x{plane[1]}.cif"
    write(os.path.join(out_dir, file_name), slab)
    
    # Print details for verification
    num_O = len([atom for atom in slab if atom.symbol == 'O'])
    num_H = len([atom for atom in slab if atom.symbol == 'H'])
    print(f"Created {surface_name} surface with dimensions {plane}")
    print(f"Contains {num_O} oxygen atoms and {num_H} hydrogen atoms")
    
    return slab