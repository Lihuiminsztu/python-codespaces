from ase.spacegroup import crystal
from ase.io import write
import os
def buildIh():
    '''functions:Build Ih unit cell.
    The dataset comes from https://next-gen.materialsproject.org/materials/mp-703459?formula=H2O.
    After finishing building the Ih unit cell, the structure is written to Ihunitcell.cif in your out directory.
    '''
    Ihunitcell = crystal(
    symbols=['O','O','H','H','H'],
    basis=[
        [0.334231, 0.334231, 0.555157],
        [0.667414, 0.667414, 0.430407],
        [0.336017, 0.336017, 0.696031],
        [0.460401, 0.460401, 0.511393],
        [0.792978, 0.669243, 0.478506]
    ],
    spacegroup=185,
    cellpar=[7.50, 7.50, 7.06, 90, 90, 120]
    )

    out = "/home/lihuimin/projects/BEofSpecialMoleculeOnIh/CO/out"
    write(os.path.join(out,'Ihunitcell.cif'), Ihunitcell)
    return Ihunitcell