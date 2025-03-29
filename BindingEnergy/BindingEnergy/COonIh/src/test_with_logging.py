# version 1.0
# import sys
# import os
# import json
# import argparse
# import datetime
# from ase.visualize import view
# from ase.build import molecule
# import numpy as np

# sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
# from source.BuildIhSpecial import buildIh001
# from source.BuildIh import buildIh
# from source.absorption import find_absorption_site, place_molecule_on_slab
# from source.calculate import optimazation, calculate

# def get_next_record_number():
#     """Get the next available record number"""
#     records_dir = "/home/lihuimin/projects/BEofSpecialMoleculeOnIh/CO/out/records"
    
#     # Create directory if it doesn't exist
#     if not os.path.exists(records_dir):
#         os.makedirs(records_dir)
#         return 1
    
#     # Find all existing record files
#     existing_files = [f for f in os.listdir(records_dir) 
#                      if f.startswith("record") and f.endswith(".json")]
    
#     if not existing_files:
#         return 1
    
#     # Extract numbers from filenames
#     numbers = []
#     for filename in existing_files:
#         try:
#             # Extract number between "record" and ".json"
#             num_str = filename[6:-5]  
#             if num_str:
#                 numbers.append(int(num_str))
#         except ValueError:
#             continue
    
#     if not numbers:
#         return 1
    
#     # Return the next number in sequence
#     return max(numbers) + 1

# def save_record(params, results):
#     """Save test run record to the next available numbered file"""
#     record_num = get_next_record_number()
#     records_dir = "/home/lihuimin/projects/BEofSpecialMoleculeOnIh/CO/out/records"
    
#     # Create record data
#     record = {
#         "run_number": record_num,
#         "timestamp": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
#         "parameters": params,
#         "results": results
#     }
    
#     # Save to file
#     filename = f"record{record_num}.json"
#     filepath = os.path.join(records_dir, filename)
    
#     with open(filepath, 'w') as f:
#         json.dump(record, f, indent=2)
    
#     print(f"Test record saved to: {filepath}")
#     return filepath

# def run_simulation(temperature, site_type, distance, calculator, repeat_size=(1,1,1), record=True):
#     """
#     Run a complete adsorption simulation with the specified parameters
    
#     Parameters:
#     ----------
#     temperature : float
#         Electronic temperature in K
#     site_type : str
#         Adsorption site type ('ontop', 'bridge', 'hollow')
#     distance : float
#         Initial adsorption distance in Å
#     calculator : str
#         Calculator to use ('GFN1-xTB', 'GFN2-xTB')
#     repeat_size : tuple
#         Size of the repeated cell (nx, ny, nz)
#     record : bool
#         Whether to record results to the log file
        
#     Returns:
#     ----------
#     result_dict : dict
#         Dictionary with simulation parameters and results
#     """
#     start_time = datetime.datetime.now()
    
#     # Create parameters dictionary for recording
#     params = {
#         "date": start_time.strftime("%Y-%m-%d %H:%M:%S"),
#         "temperature": temperature,
#         "site_type": site_type,
#         "distance": distance,
#         "calculator": calculator,
#         "repeat_size": repeat_size,
#         "vacuum_height": 10.0  # Default vacuum height in Å
#     }
    
#     print(f"Starting simulation with parameters: {params}")
    
#     try:
#         # Build the slab
#         Ih1 = buildIh001()
        
#         # Apply repeat if specified
#         if repeat_size != (1,1,1):
#             Ih1 = Ih1.repeat(repeat_size)
        
#         # Build CO molecule
#         co = molecule("CO")
        
#         # Find adsorption site with specified distance
#         pos = find_absorption_site(Ih1, site_type=site_type, distance=distance)
        
#         # Place molecule on slab
#         complexsystem = place_molecule_on_slab(Ih1, co, pos)
        
#         # Create output directory for this run
#         run_id = f"{site_type}_{distance}A_{calculator}_{temperature}K"
#         out_dir = f"/home/lihuimin/projects/BEofSpecialMoleculeOnIh/CO/out/{run_id}"
#         if not os.path.exists(out_dir):
#             os.makedirs(out_dir)
        
#         # Save initial structure
#         from ase.io import write
#         write(os.path.join(out_dir, "initial.cif"), complexsystem)
        
#         # Optimize with specified parameters
#         complexsystem = optimazation(complexsystem, calculator=calculator, 
#                                    electronic_temperature=temperature,
#                                    logfile=os.path.join(out_dir, "opt.log"))
        
#         # Save optimized structure
#         write(os.path.join(out_dir, "optimized.cif"), complexsystem)
        
#         # Calculate binding energy
#         BE = calculate(co, Ih1, complexsystem, calculator=calculator, 
#                      electronic_temperature=temperature)
        
#         success = True
        
#     except Exception as e:
#         BE = None
#         error_message = str(e)
#         success = False
#         out_dir = None
    
#     end_time = datetime.datetime.now()
#     duration = (end_time - start_time).total_seconds()
    
#     # Create result dictionary
#     results = {
#         "binding_energy": None if BE is None else float(BE),
#         "binding_energy_formatted": None if BE is None else f"{BE:.2f} eV",
#         "duration_seconds": duration,
#         "duration_formatted": str(datetime.timedelta(seconds=int(duration))),
#         "output_directory": out_dir,
#         "success": success
#     }
    
#     if not success:
#         results["error"] = error_message
    
#     # Print results to console
#     if success:
#         print(f"\nSimulation Results ({run_id}):")
#         print(f"Binding energy: {BE:.2f} eV")
#     else:
#         print(f"Calculation failed: {error_message}")
    
#     print(f"Calculation time: {datetime.timedelta(seconds=int(duration))}")
    
#     if out_dir:
#         print(f"Output saved to: {out_dir}")
    
#     # Save individual record
#     record_file = save_record(params, results)
    
#     # Record to combined log file if requested
#     if record:
#         log_file = "/home/lihuimin/projects/BEofSpecialMoleculeOnIh/CO/simulation_log.json"
        
#         # Read existing log if it exists
#         if os.path.exists(log_file):
#             with open(log_file, 'r') as f:
#                 try:
#                     log_data = json.load(f)
#                 except json.JSONDecodeError:
#                     log_data = {"simulations": []}
#         else:
#             log_data = {"simulations": []}
        
#         # Create entry for combined log
#         combined_entry = {
#             "run_number": results.get("run_number", get_next_record_number() - 1),
#             "date": params["date"],
#             "parameters": params,
#             "results": results
#         }
        
#         # Add new result
#         log_data["simulations"].append(combined_entry)
        
#         # Write updated log
#         with open(log_file, 'w') as f:
#             json.dump(log_data, f, indent=2)
        
#         print(f"Results also added to combined log: {log_file}")
    
#     return results

# def main():
#     # Set up command line arguments
#     parser = argparse.ArgumentParser(description='Run CO adsorption simulation on ice Ih(001)')
    
#     parser.add_argument('--temp', type=float, default=300.0,
#                         help='Electronic temperature in K (default: 300.0)')
    
#     parser.add_argument('--site', type=str, default='ontop', 
#                         choices=['ontop', 'bridge', 'hollow'],
#                         help='Adsorption site type (default: ontop)')
    
#     parser.add_argument('--dist', type=float, default=2.5,
#                         help='Initial adsorption distance in Å (default: 2.5)')
    
#     parser.add_argument('--calc', type=str, default='GFN2-xTB',
#                         choices=['GFN1-xTB', 'GFN2-xTB'],
#                         help='Calculator to use (default: GFN2-xTB)')
    
#     parser.add_argument('--repeat', type=str, default='1,1,1',
#                         help='Cell repeat size as "nx,ny,nz" (default: 1,1,1)')
    
#     parser.add_argument('--no-log', action='store_true',
#                         help='Disable saving to the combined log file')
    
#     args = parser.parse_args()
    
#     # Parse repeat size
#     repeat_size = tuple(map(int, args.repeat.split(',')))
    
#     # Run simulation with specified parameters
#     run_simulation(
#         temperature=args.temp,
#         site_type=args.site,
#         distance=args.dist,
#         calculator=args.calc,
#         repeat_size=repeat_size,
#         record=not args.no_log
#     )

# if __name__ == "__main__":
#     main()


import sys
import os
import json
import argparse
import datetime
from ase.constraints import FixAtoms
from ase.visualize import view
from ase.build import molecule
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)))
from source.BuildIhSpecial import buildIhSpecial
from source.absorption import find_absorption_site, place_molecule_on_slab
from source.calculate import optimazation, calculate

def get_next_record_number():
    """Get the next available record number"""
    records_dir = "/home/lihuimin/projects/BEofSpecialMoleculeOnIh/CO/out/records"
    
    # Create directory if it doesn't exist
    if not os.path.exists(records_dir):
        os.makedirs(records_dir)
        return 1
    
    # Find all existing record files
    existing_files = [f for f in os.listdir(records_dir) 
                     if f.startswith("record") and f.endswith(".json")]
    
    if not existing_files:
        return 1
    
    # Extract numbers from filenames
    numbers = []
    for filename in existing_files:
        try:
            # Extract number between "record" and ".json"
            num_str = filename[6:-5]  
            if num_str:
                numbers.append(int(num_str))
        except ValueError:
            continue
    
    if not numbers:
        return 1
    
    # Return the next number in sequence
    return max(numbers) + 1

def save_record(params, results):
    """Save test run record to the next available numbered file"""
    record_num = get_next_record_number()
    records_dir = "/home/lihuimin/projects/BEofSpecialMoleculeOnIh/CO/out/records"
    
    # Create record data
    record = {
        "run_number": record_num,
        "timestamp": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "parameters": params,
        "results": results
    }
    
    # Save to file
    filename = f"record{record_num}.json"
    filepath = os.path.join(records_dir, filename)
    
    with open(filepath, 'w') as f:
        json.dump(record, f, indent=2)
    
    print(f"Test record saved to: {filepath}")
    return filepath

def run_simulation(miller_indices, plane_size, temperature, site_type, atom_indices, 
                   distance, calculator, record=True):
    """
    Run a complete adsorption simulation with the specified parameters
    
    Parameters:
    ----------
    miller_indices : tuple
        Miller indices for the ice surface (e.g., (0,0,1))
    plane_size : tuple
        Size of the repeated cell (nx, ny, nz)
    temperature : float
        Electronic temperature in K
    site_type : str
        Adsorption site type ('ontop', 'bridge', 'hollow')
    atom_indices : int or list
        Atom indices to use for adsorption site
    distance : float
        Initial adsorption distance in Å
    calculator : str
        Calculator to use ('GFN1-xTB', 'GFN2-xTB')
    record : bool
        Whether to record results to the log file
        
    Returns:
    ----------
    result_dict : dict
        Dictionary with simulation parameters and results
    """
    start_time = datetime.datetime.now()
    
    # Create parameters dictionary for recording
    params = {
        "date": start_time.strftime("%Y-%m-%d %H:%M:%S"),
        "miller_indices": list(miller_indices),
        "plane_size": list(plane_size),
        "temperature": temperature,
        "site_type": site_type,
        "atom_indices": atom_indices,
        "distance": distance,
        "calculator": calculator
    }
    
    print(f"Starting simulation with parameters: {params}")
    
    try:
        # Build the slab with specified Miller indices and plane size
        Ih1 = buildIhSpecial(miller_indices, plane=plane_size)
        
        # Build CO molecule
        co = molecule("CO")
        
        # Find adsorption site with specified parameters
        pos = find_absorption_site(Ih1, site_type=site_type, atom_indices=atom_indices, distance=distance)
        
        # Place molecule on slab
        complexsystem = place_molecule_on_slab(Ih1, co, pos)
        
        # Create output directory for this run
        indices_str = f"{miller_indices[0]}{miller_indices[1]}{miller_indices[2]}"
        plane_str = f"{plane_size[0]}x{plane_size[1]}x{plane_size[2]}"
        run_id = f"Ih{indices_str}_{plane_str}_{site_type}_d{distance}_{calculator}"
        out_dir = f"/home/lihuimin/projects/BEofSpecialMoleculeOnIh/CO/out/{run_id}"
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        
        # Save initial structure
        from ase.io import write
        write(os.path.join(out_dir, "initial.cif"), complexsystem)
        
        # Optimize with specified parameters
        complexsystem = optimazation(complexsystem, calculator=calculator, 
                                   elec_temp=temperature)
        
        # Save optimized structure
        write(os.path.join(out_dir, "optimized.cif"), complexsystem)
        
        # Calculate binding energy
        BE = calculate(co, Ih1, complexsystem, calculator=calculator,electronic_temperature=temperature)
        
        success = True
        
    except Exception as e:
        BE = None
        error_message = str(e)
        success = False
        out_dir = None
    
    end_time = datetime.datetime.now()
    duration = (end_time - start_time).total_seconds()
    
    # Create result dictionary
    results = {
        "binding_energy": None if BE is None else float(BE),
        "binding_energy_formatted": None if BE is None else f"{BE:.2f} eV",
        "duration_seconds": duration,
        "duration_formatted": str(datetime.timedelta(seconds=int(duration))),
        "output_directory": out_dir,
        "success": success
    }
    
    if not success:
        results["error"] = error_message
    
    # Print results to console
    if success:
        print(f"\nSimulation Results ({run_id}):")
        print(f"Binding energy: {BE:.2f} eV")
    else:
        print(f"Calculation failed: {error_message}")
    
    print(f"Calculation time: {datetime.timedelta(seconds=int(duration))}")
    
    if out_dir:
        print(f"Output saved to: {out_dir}")
    
    # Save individual record
    record_file = save_record(params, results)
    
    # Record to combined log file if requested
    if record:
        log_file = "/home/lihuimin/projects/BEofSpecialMoleculeOnIh/CO/simulation_log.json"
        
        # Read existing log if it exists
        if os.path.exists(log_file):
            with open(log_file, 'r') as f:
                try:
                    log_data = json.load(f)
                except json.JSONDecodeError:
                    log_data = {"simulations": []}
        else:
            log_data = {"simulations": []}
        
        # Create entry for combined log
        combined_entry = {
            "run_number": get_next_record_number() - 1,
            "date": params["date"],
            "parameters": params,
            "results": results
        }
        
        # Add new result
        log_data["simulations"].append(combined_entry)
        
        # Write updated log
        with open(log_file, 'w') as f:
            json.dump(log_data, f, indent=2)
        
        print(f"Results also added to combined log: {log_file}")
    
    return results

def main():
    # Set up command line arguments
    parser = argparse.ArgumentParser(description='Run CO adsorption simulation on ice Ih surface')
    
    parser.add_argument('--miller', type=str, default='0,0,1',
                       help='Miller indices as "h,k,l" (default: 0,0,1)')
    
    parser.add_argument('--plane', type=str, default='1,1,1',
                       help='Plane size as "nx,ny,nz" (default: 1,1,1)')
    
    parser.add_argument('--temp', type=float, default=300.0,
                        help='Electronic temperature in K (default: 300.0)')
    
    parser.add_argument('--site', type=str, default='ontop', 
                        choices=['ontop', 'bridge', 'hollow'],
                        help='Adsorption site type (default: ontop)')
    
    parser.add_argument('--atom', type=str, default='0',
                        help='Atom indices for adsorption site (default: 0)')
    
    parser.add_argument('--dist', type=float, default=2.5,
                        help='Initial adsorption distance in Å (default: 2.5)')
    
    parser.add_argument('--calc', type=str, default='GFN2-xTB',
                        choices=['GFN1-xTB', 'GFN2-xTB'],
                        help='Calculator to use (default: GFN2-xTB)')
    
    parser.add_argument('--no-log', action='store_true',
                        help='Disable saving to the combined log file')
    
    args = parser.parse_args()
    
    # Parse Miller indices
    miller_indices = tuple(map(int, args.miller.split(',')))
    
    # Parse plane size
    plane_size = tuple(map(int, args.plane.split(',')))
    
    # Parse atom indices (could be single integer or list)
    try:
        if ',' in args.atom:
            atom_indices = list(map(int, args.atom.split(',')))
        else:
            atom_indices = int(args.atom)
    except ValueError:
        print("Error: atom indices must be integers")
        return
    
    # Run simulation with specified parameters
    run_simulation(
        miller_indices=miller_indices,
        plane_size=plane_size,
        temperature=args.temp,
        site_type=args.site,
        atom_indices=atom_indices,
        distance=args.dist,
        calculator=args.calc,
        record=not args.no_log,
    )

if __name__ == "__main__":
    main()