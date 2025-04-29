# Amorphous Ih Studies

## Overview
This files open a new direction of Building ASW(Amorphous Solid Water) work. As of now, i have updated two scripts both used to change the Ih crystal come from BuildIh function in the other directory "Crystalline", however, one of these scripts run by Tetrahedral methods and another one is try on Hydrogen bond network methods which basic on Ice Rules(Bernal-Fowler Rules).

## Main Features

I am trying to interim building object from Crystalline Ih to Amorphous Ih. I used the tetrahedral method to change the orientation of H atoms connected to bilayer oxygen atoms. In "TetrahedralMethod.py" python script, codes can come ture three main features:
- 1. Moving positions of Hydrogen atoms to special tetrahedral vertices.
- 2. Orient top layer hydrogen facing directions
- 3. Analyze atomic distribution in the system

By the way, I'm trying on randomizing H atoms facing orientation when typing "HNetMethod.py" file. Till now, codes had made sure the Bernal-Fowler ice rules are followed during hydrogen randomization. The script implements a hydrogen bond network approach that:
- Creates a neighbor list to identify O-O connections
- Randomly assigns hydrogen atoms while ensuring each oxygen forms exactly 2 covalent bonds


## Dependecies 
- Atomic Simulation Environment(ASE)
- Numpy 

## Notes

**if u try clone this code, please change the default saved filepath.**

