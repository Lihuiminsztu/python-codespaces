# region 中文版本
# date: 2025-4-24

# from ase.spacegroup import crystal
# from ase.io import  write
# import os
# import numpy as np


# def calculate_tetrahedral_vertices(center, radius):
#     '''计算以center为中心的四面体顶点坐标
    
#     Args:
#         center: 四面体中心坐标
#         radius: 四面体中心到顶点的距离
    
#     Returns:
#         四个顶点的坐标数组
#     '''
#     # 定义四面体的四个顶点向量(单位向量)
#     v1 = np.array([1, 1, 1])
#     v2 = np.array([-1, -1, 1])
#     v3 = np.array([-1, 1, -1])
#     v4 = np.array([1, -1, -1])

#     # 标准化向量
#     vertices = np.array([v1, v2, v3, v4])
#     vertices = vertices / np.linalg.norm(vertices, axis=1)[:, np.newaxis] * radius
    
#     # 将顶点移动到中心位置
#     vertices += center
    
#     return vertices

# def identify_middle_layer_oxygens(atoms, z_min=0.4, z_max=0.6):
#     '''识别中间层的氧原子
    
#     Args:
#         atoms: ASE原子对象
#         z_min: 中间层z坐标最小值(分数坐标)
#         z_max: 中间层z坐标最大值(分数坐标)
    
#     Returns:
#         中间层氧原子的索引列表
#     '''
#     # 获取晶胞参数
#     cell = atoms.get_cell()
    
#     # 获取所有氧原子索引
#     oxygen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'O']
    
#     # 筛选中间层的氧原子
#     middle_layer_oxygens = []
#     for i in oxygen_indices:
#         # 获取原子的分数坐标
#         scaled_pos = atoms.get_scaled_positions()[i]
#         if z_min <= scaled_pos[2] <= z_max:
#             middle_layer_oxygens.append(i)
    
#     return middle_layer_oxygens
# def identify_top_layer_oxygens(atoms, z_min=0.8, z_max=1.0):
#     '''识别顶层的氧原子
    
#     Args:
#         atoms: ASE原子对象
#         z_min: 顶层z坐标最小值(分数坐标)
#         z_max: 顶层z坐标最大值(分数坐标)
    
#     Returns:
#         顶层氧原子的索引列表
#     '''
#     # 获取晶胞参数
#     cell = atoms.get_cell()
    
#     # 获取所有氧原子索引
#     oxygen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'O']
    
#     # 筛选顶层的氧原子
#     top_layer_oxygens = []
#     for i in oxygen_indices:
#         # 获取原子的分数坐标
#         scaled_pos = atoms.get_scaled_positions()[i]
#         if z_min <= scaled_pos[2] <= z_max:
#             top_layer_oxygens.append(i)
    
#     return top_layer_oxygens


# def find_hydrogen_connected_to_oxygen(atoms, oxygen_index, cutoff=1.1):
#     '''找出与给定氧原子相连的氢原子
    
#     Args:
#         atoms: ASE原子对象
#         oxygen_index: 氧原子的索引
#         cutoff: 氧-氢键长的截断值(Å)
    
#     Returns:
#         与该氧原子相连的氢原子索引列表
#     '''
#     # 获取氧原子位置
#     oxygen_pos = atoms.positions[oxygen_index]
    
#     # 获取所有氢原子索引
#     hydrogen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'H']
    
#     # 找出距离小于cutoff的氢原子
#     connected_hydrogens = []
#     for h_idx in hydrogen_indices:
#         h_pos = atoms.positions[h_idx]
#         distance = np.linalg.norm(h_pos - oxygen_pos)
#         if distance < cutoff:
#             connected_hydrogens.append(h_idx)
    
#     return connected_hydrogens

# def randomize_middle_layer_hydrogen(atoms):
#     '''随机化中间层氢原子的朝向
#     使用四面体顶点法重新分配与中间层氧原子相连的氢原子的位置
    
#     '''
#     # 复制原子对象
#     modified_atoms = atoms.copy()
    
#     # 识别中间层的氧原子
#     middle_layer_oxygens = identify_middle_layer_oxygens(modified_atoms)
#     print(f"找到{len(middle_layer_oxygens)}个中间层氧原子")
    
#     # 氧-氢键长
#     oh_bond_length = 0.84632  # Å
    
#     # 处理每个中间层氧原子
#     for o_idx in middle_layer_oxygens:
#         # 找出与该氧原子相连的氢原子
#         connected_hydrogens = find_hydrogen_connected_to_oxygen(modified_atoms, o_idx)
        
#         # 如果找到了相连的氢原子
#         if len(connected_hydrogens) > 0:
#             # 获取氧原子位置
#             o_pos = modified_atoms.positions[o_idx]
            
#             # 计算四面体顶点
#             tetrahedral_vertices = calculate_tetrahedral_vertices(o_pos, oh_bond_length)
            
#             # 随机选择两个顶点用于放置氢原子
#             h_indices = np.random.choice(range(4), size=min(len(connected_hydrogens), 2), replace=False)
#             h_positions = tetrahedral_vertices[h_indices]
            
#             # 更新氢原子位置
#             for i, h_idx in enumerate(connected_hydrogens):
#                 if i < len(h_positions):
#                     modified_atoms.positions[h_idx] = h_positions[i]
    
#     return modified_atoms
    
# def orient_toplayer_hydrogens(atoms):
#     '''将顶层其中一个氢原子朝向z轴负方向
#     Returns:
#         obeject 
#     '''
#     # 复制原子对象
#     modified_atoms = atoms.copy()
    
#     top_layer_oxygens = identify_top_layer_oxygens(modified_atoms)
#     print(f"找到{len(top_layer_oxygens)}个顶层氧原子")
#     # 氧-氢键长
#     oh_bond_length = 0.84632  # Å
#     # 处理其中一个顶层氢原子朝向z轴负方向
#     for o_idx in top_layer_oxygens:
#         # 找出与该氧原子相连的氢原子
#         connected_hydrogens = find_hydrogen_connected_to_oxygen(modified_atoms, o_idx)
        
#         # 如果找到了相连的氢原子
#         if len(connected_hydrogens) > 0:
#             # 获取氧原子位置
#             o_pos = modified_atoms.positions[o_idx]
            
#             # 计算新的氢原子位置
#             new_h_pos = np.array([o_pos[0], o_pos[1], o_pos[2] - 1.0053])
            
#             # 更新氢原子位置
#             for h_idx in connected_hydrogens:
#                 modified_atoms.positions[h_idx] = new_h_pos
#                 break  
#     return modified_atoms
    
# def analyze_hydrogen_distribution(atoms):
#     '''分析氢原子分布情况
#     '''
#     # 获取所有氧原子和氢原子索引
#     oxygen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'O']
#     hydrogen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'H']
    
#     # 统计每个氧原子连接的氢原子数量
#     o_h_counts = {}
#     for o_idx in oxygen_indices:
#         connected_hydrogens = find_hydrogen_connected_to_oxygen(atoms, o_idx)
#         o_h_counts[o_idx] = len(connected_hydrogens)
    
#     print("氧原子连接的氢原子统计:")
#     print(f"氧原子总数: {len(oxygen_indices)}")
#     print(f"氢原子总数: {len(hydrogen_indices)}")
    
#     # 统计连接不同数量氢原子的氧原子数
#     h_count_distribution = {}
#     for count in o_h_counts.values():
#         if count not in h_count_distribution:
#             h_count_distribution[count] = 0
#         h_count_distribution[count] += 1
    
#     for count, num_oxygens in sorted(h_count_distribution.items()):
#         print(f"连接{count}个氢原子的氧原子数: {num_oxygens}")

# if __name__ == "__main__":
#     # 设置输出目录
#     output_dir = "/home/lihuimin/projects/BEofSpecialMoleculeOnIh/CO/out"
#     # 搭建Ih unitcell
#     Ihunitcell = crystal(
#     symbols=['O','O','H','H','H'],
#     basis=[
#         [0.334231, 0.334231, 0.555157],
#         [0.667414, 0.667414, 0.430407],
#         [0.336017, 0.336017, 0.696031],
#         [0.460401, 0.460401, 0.511393],
#         [0.792978, 0.669243, 0.478506]
#     ],
#     spacegroup=185,
#     cellpar=[7.50, 7.50, 7.06, 90, 90, 120]
#     )

#     print("初始结构分析:")
#     analyze_hydrogen_distribution(Ihunitcell)
    
#     # 随机化中间层氢原子
#     randomized_Ihunitcell = randomize_middle_layer_hydrogen(Ihunitcell)
#     randomized_Ihunitcell = orient_toplayer_hydrogens(randomized_Ihunitcell)
    
#     print("\n随机化后结构分析:")
#     analyze_hydrogen_distribution(randomized_Ihunitcell)
#     # 保存结构
#     write(os.path.join(output_dir, 'Ih_middle_layer_randomized.cif'), randomized_Ihunitcell)

# endregion

# date: 2025-4-24

from ase.spacegroup import crystal
from ase.io import write
import os
import numpy as np


def calculate_tetrahedral_vertices(center, radius):
    '''Calculate the coordinates of tetrahedral vertices centered at the given point
    
    Args:
        center: Tetrahedral center coordinates
        radius: Distance from center to vertices
    
    Returns:
        Array of four vertex coordinates
    '''
    # Define four tetrahedral vertex vectors (unit vectors)
    v1 = np.array([1, 1, 1])
    v2 = np.array([-1, -1, 1])
    v3 = np.array([-1, 1, -1])
    v4 = np.array([1, -1, -1])

    # Normalize vectors
    vertices = np.array([v1, v2, v3, v4])
    vertices = vertices / np.linalg.norm(vertices, axis=1)[:, np.newaxis] * radius
    
    # Move vertices to center position
    vertices += center
    
    return vertices

def identify_middle_layer_oxygens(atoms, z_min=0.4, z_max=0.6):
    '''Identify oxygen atoms in the middle layer
    
    Args:
        atoms: ASE atoms object
        z_min: Minimum z-coordinate of middle layer (fractional)
        z_max: Maximum z-coordinate of middle layer (fractional)
    
    Returns:
        List of indices of middle layer oxygen atoms
    '''
    # Get cell parameters
    cell = atoms.get_cell()
    
    # Get all oxygen atom indices
    oxygen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'O']
    
    # Filter oxygen atoms in the middle layer
    middle_layer_oxygens = []
    for i in oxygen_indices:
        # Get fractional coordinates of the atom
        scaled_pos = atoms.get_scaled_positions()[i]
        if z_min <= scaled_pos[2] <= z_max:
            middle_layer_oxygens.append(i)
    
    return middle_layer_oxygens

def identify_top_layer_oxygens(atoms, z_min=0.8, z_max=1.0):
    '''Identify oxygen atoms in the top layer
    
    Args:
        atoms: ASE atoms object
        z_min: Minimum z-coordinate of top layer (fractional)
        z_max: Maximum z-coordinate of top layer (fractional)
    
    Returns:
        List of indices of top layer oxygen atoms
    '''
    # Get cell parameters
    cell = atoms.get_cell()
    
    # Get all oxygen atom indices
    oxygen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'O']
    
    # Filter oxygen atoms in the top layer
    top_layer_oxygens = []
    for i in oxygen_indices:
        # Get fractional coordinates of the atom
        scaled_pos = atoms.get_scaled_positions()[i]
        if z_min <= scaled_pos[2] <= z_max:
            top_layer_oxygens.append(i)
    
    return top_layer_oxygens


def find_hydrogen_connected_to_oxygen(atoms, oxygen_index, cutoff=1.1):
    '''Find hydrogen atoms connected to a given oxygen atom
    
    Args:
        atoms: ASE atoms object
        oxygen_index: Index of the oxygen atom
        cutoff: O-H bond length cutoff value (Å)
    
    Returns:
        List of hydrogen atom indices connected to the oxygen
    '''
    # Get oxygen atom position
    oxygen_pos = atoms.positions[oxygen_index]
    
    # Get all hydrogen atom indices
    hydrogen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'H']
    
    # Find hydrogen atoms with distance less than cutoff
    connected_hydrogens = []
    for h_idx in hydrogen_indices:
        h_pos = atoms.positions[h_idx]
        distance = np.linalg.norm(h_pos - oxygen_pos)
        if distance < cutoff:
            connected_hydrogens.append(h_idx)
    
    return connected_hydrogens

def randomize_middle_layer_hydrogen(atoms):
    '''Randomize the orientation of middle layer hydrogen atoms
    Reassign positions of hydrogen atoms connected to middle layer oxygen atoms
    using tetrahedral vertex method
    
    '''
    # Copy atoms object
    modified_atoms = atoms.copy()
    
    # Identify oxygen atoms in the middle layer
    middle_layer_oxygens = identify_middle_layer_oxygens(modified_atoms)
    print(f"Found {len(middle_layer_oxygens)} oxygen atoms in middle layer")
    
    # O-H bond length
    oh_bond_length = 0.84632  # Å
    
    # Process each middle layer oxygen atom
    for o_idx in middle_layer_oxygens:
        # Find hydrogen atoms connected to this oxygen
        connected_hydrogens = find_hydrogen_connected_to_oxygen(modified_atoms, o_idx)
        
        # If connected hydrogen atoms are found
        if len(connected_hydrogens) > 0:
            # Get oxygen atom position
            o_pos = modified_atoms.positions[o_idx]
            
            # Calculate tetrahedral vertices
            tetrahedral_vertices = calculate_tetrahedral_vertices(o_pos, oh_bond_length)
            
            # Randomly select two vertices for hydrogen atoms
            h_indices = np.random.choice(range(4), size=min(len(connected_hydrogens), 2), replace=False)
            h_positions = tetrahedral_vertices[h_indices]
            
            # Update hydrogen atom positions
            for i, h_idx in enumerate(connected_hydrogens):
                if i < len(h_positions):
                    modified_atoms.positions[h_idx] = h_positions[i]
    
    return modified_atoms
    
def orient_toplayer_hydrogens(atoms):
    '''Orient one hydrogen atom from top layer oxygen atoms towards negative z direction
    '''
    # Copy atoms object
    modified_atoms = atoms.copy()
    
    top_layer_oxygens = identify_top_layer_oxygens(modified_atoms)
    print(f"Found {len(top_layer_oxygens)} oxygen atoms in top layer")
    
    # O-H bond length
    oh_bond_length = 0.84632  # Å
    
    # Orient one hydrogen atom towards negative z direction
    for o_idx in top_layer_oxygens:
        # Find hydrogen atoms connected to this oxygen
        connected_hydrogens = find_hydrogen_connected_to_oxygen(modified_atoms, o_idx)
        
        # If connected hydrogen atoms are found
        if len(connected_hydrogens) > 0:
            # Get oxygen atom position
            o_pos = modified_atoms.positions[o_idx]
            
            # Calculate new hydrogen atom position
            new_h_pos = np.array([o_pos[0], o_pos[1], o_pos[2] - 1.0053])
            
            # Update hydrogen atom position
            for h_idx in connected_hydrogens:
                modified_atoms.positions[h_idx] = new_h_pos
                break  
    return modified_atoms
    
def analyze_hydrogen_distribution(atoms):
    '''Analyze hydrogen atom distribution
    '''
    # Get all oxygen and hydrogen atom indices
    oxygen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'O']
    hydrogen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'H']
    
    # Count number of hydrogen atoms connected to each oxygen
    o_h_counts = {}
    for o_idx in oxygen_indices:
        connected_hydrogens = find_hydrogen_connected_to_oxygen(atoms, o_idx)
        o_h_counts[o_idx] = len(connected_hydrogens)
    
    print("Oxygen-hydrogen connection statistics:")
    print(f"Total oxygen atoms: {len(oxygen_indices)}")
    print(f"Total hydrogen atoms: {len(hydrogen_indices)}")
    
    # Count distribution of oxygen atoms with different numbers of hydrogen
    h_count_distribution = {}
    for count in o_h_counts.values():
        if count not in h_count_distribution:
            h_count_distribution[count] = 0
        h_count_distribution[count] += 1
    
    for count, num_oxygens in sorted(h_count_distribution.items()):
        print(f"Oxygen atoms connected to {count} hydrogen atoms: {num_oxygens}")

if __name__ == "__main__":
    # Set output directory
    output_dir = "/home/lihuimin/projects/BEofSpecialMoleculeOnIh/CO/out"
    # Build Ih unit cell
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

    print("Initial structure analysis:")
    analyze_hydrogen_distribution(Ihunitcell)
    
    # Randomize middle layer hydrogen atoms
    randomized_Ihunitcell = randomize_middle_layer_hydrogen(Ihunitcell)
    randomized_Ihunitcell = orient_toplayer_hydrogens(randomized_Ihunitcell)
    
    print("\nStructure analysis after randomization:")
    analyze_hydrogen_distribution(randomized_Ihunitcell)
    # Save structure
    write(os.path.join(output_dir, 'Ih_middle_layer_randomized.cif'), randomized_Ihunitcell)
    