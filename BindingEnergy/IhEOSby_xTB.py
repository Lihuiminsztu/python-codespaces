from tblite.ase import TBLite  
from ase.spacegroup import crystal
from ase.eos import EquationOfState
import numpy as np
import matplotlib.pyplot as plt

a = 7.50
b = 7.50
c = 7.06
atoms = crystal(
    symbols=['O', 'O', 'H', 'H', 'H'],
    basis=[
        (0.334231, 0.334231, 0.555157),  # O1
        (0.667414, 0.667414, 0.430407),  # O2
        (0.336017, 0.336017, 0.696031),  # H1
        (0.460401, 0.460401, 0.511393),  # H2
        (0.792978, 0.669243, 0.478506)   # H3
    ],
    spacegroup=185,  # P6_3cm
    cellpar=[a, b, c, 90, 90, 120]
)

# 设置xTB计算器


#拟合不同体积能量值以得到平衡状态下体积、能量和体积弹性模量
scale_factors = np.linspace(0.90, 1.10, 11)
volumes = []    
energies = []
print("Calculating EOS")
for scale in scale_factors:
    atoms_scaled = atoms.copy()
    atoms_scaled.set_cell([a*scale, b*scale, c*scale], scale_atoms=True)
    calculator = TBLite(
        method="GFN2-xTB",  
    )
    atoms_scaled.calc = calculator
    print(scale)
    print("calculating能量、体积")
    energy = atoms_scaled.get_potential_energy()
    volume = atoms_scaled.get_volume()
    volumes.append(volume)
    energies.append(energy)
print("Fitting EOS")
eos = EquationOfState(volumes, energies, eos="birchmurnaghan")
v0, e0, B = eos.fit()
eos.plot(filename="eos2.png")