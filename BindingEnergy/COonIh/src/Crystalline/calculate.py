from ase.optimize import BFGS
from tblite.ase import TBLite
from ase.io import write,read
import numpy as np
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))

def optimazation(object, calculator='GFN2-xTB',elec_temp=30):
    '''
    Optimize the structure by BFGS algorithm.
    Save the optimized structure to a file on the 'out' directory.
    '''
    # Set up the calculator with modified parameters for better SCF convergence
    if calculator == 'GFN2-xTB':
        calc = TBLite(
            method='GFN2-xTB',
            # Add SCF convergence parameters
            accuracy=1.0,  # Increase accuracy
            max_iterations=500,  # Increase max SCF iterations
            electronic_temperature=elec_temp,  # Adjust electronic temperature
            mixer=0.1  # Use a smaller mixing parameter for more stable convergence
        )
    elif calculator == 'GFN1-xTB':
        calc = TBLite(
            method='GFN1-xTB',
            accuracy=1.0,
            max_iterations=500,
            electronic_temperature=elec_temp,
            mixer=0.1
        )
    else:
        raise ValueError(f"Invalid calculator: {calculator}")
    
    # Attach calculator to the object
    object.calc = calc
    
    # Optimize the structure with gentler parameters
    opt = BFGS(object, logfile='optimization.log')
    opt.run(fmax=0.05, steps=50)  # Increase force tolerance slightly
    
    # Write the optimized structure to a file
    out = "/home/lihuimin/projects/BEofSpecialMoleculeOnIh/CO/out"
    if object == 'Ih001':
        write(os.path.join(out,'Ih001_opt.cif'), object)
    elif object == 'complexsystem':
        write(os.path.join(out,'complexsystem_opt.cif'), object)
        print("Opt file had been written in "+out+"complexsystem_opt.cif")
    else:
        write(os.path.join(out,'combined_opt.cif'), object)
        print("Opt file had been written in "+out+"combined_opt.cif")
    return object
# version1
#  def calculate(atoms,slab,object,calculator):
#     '''Calculate the binding energy of System by TBlite.
    
#     Parameters:
#     ----------
#     atoms: Goal atoms.
#     slab: be used as absorbent surface.
#     object: Atoms or Surface.
#         The object to be calculated.
#     Returns:
#     ----------
#     BE: float
#         The binding energy of atoms abosorbed on slab.
#     '''
#     # Set up the calculator
#     if calculator == 'GFN2-xTB':
#         calc = TBLite(method='GFN2-xTB')
#     elif calculator == 'GFN1-xTB':
#         calc = TBLite(method='GFN1-xTB')
#     else:
#         raise ValueError(f"Invalid calculator: {calculator}")
#     atoms.calc = calc
#     slab.calc = calc
#     object.calc = calc
#     # Calculate the binding energy
#     atomsenergy = atoms.get_potential_energy()
#     slabenergy = slab.get_potential_energy()
#     objectenergy = object.get_potential_energy()
#     BE = atomsenergy + slabenergy - objectenergy
#     return BE
def calculate(molecule, slab, complex_system, calculator='GFN2-xTB', electronic_temperature=30):
    """
    Calculate binding energy between molecule and slab.
    
    Parameters:
    ----------
    molecule : Atoms
        The isolated molecule.
    slab : Atoms
        The isolated slab.
    complex_system : Atoms
        The combined molecule-slab system.
    calculator : str
        The calculator to use ('GFN1-xTB', 'GFN2-xTB').
    electronic_temperature : float
        The electronic temperature in K.
        
    Returns:
    ----------
    float
        The binding energy in eV.
    """
    # Create calculators with consistent parameters
    if calculator == 'GFN1-xTB':
        calc = TBLite(
            method="GFN1-xTB", 
            electronic_temperature=electronic_temperature,
            max_iterations=1500,
        )
    elif calculator == 'GFN2-xTB':
        calc = TBLite(
            method="GFN2-xTB", 
            electronic_temperature=electronic_temperature,
            max_iterations=1500,
        )
    else:
        raise ValueError(f'Invalid calculator: {calculator}')
    molecule.calc = calc
    slab.calc = calc
    complex_system.calc = calc
    try:
        # Calculate energies
        # Existing energy calculation code
        # ...
        atomsenergy = molecule.get_potential_energy()
        slabenergy = slab.get_potential_energy()
        objectenergy = complex_system.get_potential_energy()
        BE = atomsenergy + slabenergy - objectenergy
        return BE
    except Exception as e:
        raise ValueError(f'Error in energy calculation: {str(e)}')

