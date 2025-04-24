# CO Adsorption on Ice Ih Project Setup Guide

## Repository Structure and Setup

Welcome to the CO adsorption on Ice Ih project. This repository contains computational tools to study how carbon monoxide (CO) interacts with different ice Ih surfaces using tight-binding methods.

## Important Setup Notes

### File Paths Configuration

After cloning this repository, you will need to manually update several file paths in the scripts to match your local environment:

- `.cif` structure files and other outputs are saved to directories specified in the code
- The main output directory is currently set to `/home/lihuimin/projects/BEofSpecialMoleculeOnIh/CO/out/`
- All scripts use absolute paths that need to be modified based on your installation location

Look for paths like these in the code and update them:
```python
out_dir = "/home/lihuimin/projects/BEofSpecialMoleculeOnIh/CO/out"
```
## Using the Testing Module
The test_with_logging.py module provides a comprehensive way to run CO adsorption simulations with various parameters and automatically log the results.

#### Key Features:
  - Automated creation of Ice Ih surfaces with specified Miller indices
  - Support for different adsorption sites (ontop, bridge, hollow)
  - Flexible molecule orientation and positioning options
  - Record-keeping of simulation parameters and results 
  - Calculation of binding energies using GFN1-xTB and GFN2-xTB methods
###  Basic Usage

You can run a simulation from the command line:
```python 
python test_with_logging.py --miller 0,0,1 --plane 1,1,1 --temp 30 --site ontop --atom 0 --dist 2.5 --calc GFN2-xTB ....
```

### Available Options

| Option      | Description                          | Default   | Example Values       |
|-------------|--------------------------------------|-----------|----------------------|
| `--miller`  | Miller indices as "h,k,l"           | `0,0,1`   | `0,1,0`, `1,1,0`     |
| `--plane`   | Surface cell size as "nx,ny,nz"     | `1,1,1`   | `2,2,1`, `3,3,2`     |
| `--temp`    | Electronic temperature (K)          | `300.0`   | `30.0`, `500.0`      |
| `--site`    | Adsorption site type                | `ontop`   | `bridge`, `hollow`   |
| `--atom`    | Atom indices for site               | `0`       | `1`, `0,1`, `0,1,2`  |
| `--dist`    | Initial distance (Å)                | `2.5`     | `2.0`, `3.0`         |
| `--calc`    | Calculator type                     | `GFN2-xTB`| `GFN1-xTB`           |
| `--orient`  | Molecule orientation                | `vertical`| `horizontal`, `tilted`|
| `--angle`   | Rotation angle (degrees)            | `0.0`     | `45.0`, `90.0`       |
| `--face`    | Atom index facing the surface       | `0`       | `1`                  |
| `--invert`  | Invert molecular orientation        | `False`   | `--invert`           |
| `--no-log`  | Disable combined log file           | -         | -                    |

### Molecule Orientation Options

The `absorption.py` module has been enhanced with more control over how the CO molecule is positioned:

- **Vertical orientation** (default): Molecule aligned along the z-axis with specified atom facing the surface
- **Horizontal orientation**: Molecule aligned parallel to the surface with optional rotation angle
- **Tilted orientation**: Molecule tilted at specified angle from vertical position
- **Atom facing**: Control which atom in the molecule faces toward the surface (0=C, 1=O for CO)
- **Inversion**: Option to invert the molecule orientation

### Output and Logging
Each simulation produces:
- CIF files of initial and optimized structures
- A JSON record of simulation parameters and results
- Entries in a combined simulation log file

## Development Status
This project is currently in version 0.1 with the following status:
- Core functionality is implemented and working
- The Ice Ih builder module supports different surface orientations
- Binding energy calculations are functional with GFN1-xTB and GFN2-xTB methods
- Automatic logging and record-keeping is implemented
- Advanced molecule positioning and orientation controls are available

#### Upcoming Features
The following enhancements are planned for future releases:
- Support for additional adsorbate molecules
- Implementation of more advanced surface building techniques
- Integration with higher-level computational methods
- Extended analysis tools for binding energy trends
- Improved visualization capabilities

### All scripts will be regularly updated in this repository as development continues.