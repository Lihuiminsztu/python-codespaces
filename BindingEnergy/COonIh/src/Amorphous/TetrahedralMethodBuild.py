from ase.spacegroup import crystal
from ase.io import  write
import os
import numpy as np


def calculate_tetrahedral_vertices(center, radius):
    '''计算以center为中心的四面体顶点坐标
    
    Args:
        center: 四面体中心坐标
        radius: 四面体中心到顶点的距离
    
    Returns:
        四个顶点的坐标数组
    '''
    # 定义四面体的四个顶点向量(单位向量)
    v1 = np.array([1, 1, 1])
    v2 = np.array([-1, -1, 1])
    v3 = np.array([-1, 1, -1])
    v4 = np.array([1, -1, -1])

    # 标准化向量
    vertices = np.array([v1, v2, v3, v4])
    vertices = vertices / np.linalg.norm(vertices, axis=1)[:, np.newaxis] * radius
    
    # 将顶点移动到中心位置
    vertices += center
    
    return vertices

# region trying example
# def calculate_tetrahedral_vertices_v2(center, radius):
#     '''计算以center为中心的四面体顶点坐标
    
#     Args:
#         center: 四面体中心坐标
#         radius: 四面体中心到顶点的距离
    
#     Returns:
#         四个顶点的坐标数组
#     '''
#     # 定义四面体的四个顶点向量(单位向量)
#     v1 = np.array([0,0,1])
#     v2 = np.array([np.sqrt(6)/4, np.sqrt(6)/4, -1/2])
#     v3 = np.array([-np.sqrt(6)/4, np.sqrt(6)/4, -1/2])
#     v4 = np.array([-np.sqrt(6)/4, -np.sqrt(6)/4, -1/2])
#     # 标准化向量
#     vertices = np.array([v1, v2, v3, v4])
#     vertices = vertices / np.linalg.norm(vertices, axis=1)[:, np.newaxis] * radius
#     # 将顶点移动到中心位置
#     vertices += center
#     return vertices
# endregion

def calculate_tetrahedral_vertices_v2(center, radius, z_axis_oriented=True):
    '''计算以center为中心的四面体顶点坐标，适用于水分子中的氢原子位置调整
    
    Args:
        center: 四面体中心坐标（氧原子位置）
        radius: 四面体中心到顶点的距离（OH键长）
        z_axis_oriented: 是否使用z轴定向的四面体模型
    
    Returns:
        四个顶点的坐标数组
    '''
    if z_axis_oriented:
        # 定义z轴定向的四面体顶点，一个顶点沿z轴正方向，一个顶点沿z轴负方向
        v1 = np.array([0, 0, 1])  # z轴正方向
        v2 = np.array([0, 0, -1])  # z轴负方向
        # 在xy平面上的两个顶点，形成等边三角形
        v3 = np.array([np.sqrt(3)/2, 0, 0])
        v4 = np.array([-np.sqrt(3)/4, np.sqrt(3)*3/4, 0])
    else:
        # 原有的四面体顶点定义
        v1 = np.array([0, 0, 1])
        v2 = np.array([np.sqrt(6)/4, np.sqrt(6)/4, -1/2])
        v3 = np.array([-np.sqrt(6)/4, np.sqrt(6)/4, -1/2])
        v4 = np.array([-np.sqrt(6)/4, -np.sqrt(6)/4, -1/2])
    
    # 标准化向量
    vertices = np.array([v1, v2, v3, v4])
    vertices = vertices / np.linalg.norm(vertices, axis=1)[:, np.newaxis] * radius
    # 将顶点移动到中心位置
    vertices += center
    return vertices

def identify_middle_layer_oxygens(atoms, z_min=0.4, z_max=0.6):
    '''识别中间层的氧原子
    
    Args:
        atoms: ASE原子对象
        z_min: 中间层z坐标最小值(分数坐标)
        z_max: 中间层z坐标最大值(分数坐标)
    
    Returns:
        中间层氧原子的索引列表
    '''
    # 获取晶胞参数
    cell = atoms.get_cell()
    
    # 获取所有氧原子索引
    oxygen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'O']
    
    # 筛选中间层的氧原子
    middle_layer_oxygens = []
    for i in oxygen_indices:
        # 获取原子的分数坐标
        scaled_pos = atoms.get_scaled_positions()[i]
        if z_min <= scaled_pos[2] <= z_max:
            middle_layer_oxygens.append(i)
    
    return middle_layer_oxygens

def identify_top_layer_oxygens(atoms, z_min=0.8, z_max=1.0):
    '''识别顶层的氧原子
    
    Args:
        atoms: ASE原子对象
        z_min: 顶层z坐标最小值(分数坐标)
        z_max: 顶层z坐标最大值(分数坐标)
    
    Returns:
        顶层氧原子的索引列表
    '''
    # 获取晶胞参数
    cell = atoms.get_cell()
    
    # 获取所有氧原子索引
    oxygen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'O']
    
    # 筛选顶层的氧原子
    top_layer_oxygens = []
    for i in oxygen_indices:
        # 获取原子的分数坐标
        scaled_pos = atoms.get_scaled_positions()[i]
        if z_min <= scaled_pos[2] <= z_max:
            top_layer_oxygens.append(i)
    
    return top_layer_oxygens

def find_hydrogen_connected_to_oxygen(atoms, oxygen_index, cutoff=1.1):
    '''找出与给定氧原子相连的氢原子
    
    Args:
        atoms: ASE原子对象
        oxygen_index: 氧原子的索引
        cutoff: 氧-氢键长的截断值(Å)
    
    Returns:
        与该氧原子相连的氢原子索引列表
    '''
    # 获取氧原子位置
    oxygen_pos = atoms.positions[oxygen_index]
    
    # 获取所有氢原子索引
    hydrogen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'H']
    
    # 找出距离小于cutoff的氢原子
    connected_hydrogens = []
    for h_idx in hydrogen_indices:
        h_pos = atoms.positions[h_idx]
        distance = np.linalg.norm(h_pos - oxygen_pos)
        if distance < cutoff:
            connected_hydrogens.append(h_idx)
    
    return connected_hydrogens

def randomize_middle_layer_hydrogen(atoms):
    '''随机化中间层氢原子的朝向
    使用四面体顶点法重新分配与中间层氧原子相连的氢原子的位置
    
    '''
    # 复制原子对象
    modified_atoms = atoms.copy()
    
    # 识别中间层的氧原子
    middle_layer_oxygens = identify_middle_layer_oxygens(modified_atoms)
    print(f"找到{len(middle_layer_oxygens)}个中间层氧原子")
    
    # 氧-氢键长
    oh_bond_length = 0.84632  # Å
    
    # 处理每个中间层氧原子
    for o_idx in middle_layer_oxygens:
        # 找出与该氧原子相连的氢原子
        connected_hydrogens = find_hydrogen_connected_to_oxygen(modified_atoms, o_idx)
        
        # 如果找到了相连的氢原子
        if len(connected_hydrogens) > 0:
            # 获取氧原子位置
            o_pos = modified_atoms.positions[o_idx]
            
            # 计算四面体顶点
            tetrahedral_vertices = calculate_tetrahedral_vertices(o_pos, oh_bond_length)
            
            # 随机选择两个顶点用于放置氢原子
            h_indices = np.random.choice(range(4), size=min(len(connected_hydrogens), 2), replace=False)
            h_positions = tetrahedral_vertices[h_indices]
            
            # 更新氢原子位置
            for i, h_idx in enumerate(connected_hydrogens):
                if i < len(h_positions):
                    modified_atoms.positions[h_idx] = h_positions[i]
    
    return modified_atoms
    
def orient_toplayer_hydrogens(atoms):
    '''将顶层其中一个氢原子朝向z轴负方向
    Returns:
        obeject 
    '''
    # 复制原子对象
    modified_atoms = atoms.copy()
    
    top_layer_oxygens = identify_top_layer_oxygens(modified_atoms)
    print(f"找到{len(top_layer_oxygens)}个顶层氧原子")
    # 氧-氢键长
    oh_bond_length = 0.84632  # Å
    # 处理其中一个顶层氢原子朝向z轴负方向
    for o_idx in top_layer_oxygens:
        # 找出与该氧原子相连的氢原子
        connected_hydrogens = find_hydrogen_connected_to_oxygen(modified_atoms, o_idx)
        
        # 如果找到了相连的氢原子
        if len(connected_hydrogens) > 0:
            # 获取氧原子位置
            o_pos = modified_atoms.positions[o_idx]
            
            # 计算新的氢原子位置
            new_h_pos = np.array([o_pos[0], o_pos[1], o_pos[2] - 1.0053])
            
            # 更新氢原子位置
            for h_idx in connected_hydrogens:
                modified_atoms.positions[h_idx] = new_h_pos
                break  
    return modified_atoms

# region intial try. def orient_top_middle_hydrogens
# def orient_top_middle_hydrogens(atoms):
#     '''找到中间层中偏z轴正向的氧原子层的朝z轴正方向的H原子
#     ,并沿z轴转动到与氧原子相连的另一个H原子位置对称的四面体位置
#     '''
#     modified_atoms = atoms.copy()
#     # 识别O、H原子
#     hydrogen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'H']
#     oxygen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'O']
    
#     topmiddle_layer_oxygens = []
#     for i in oxygen_indices:
#         # 获取原子的分数坐标
#         scaled_pos = atoms.get_scaled_positions()[i]
#         if 0.45 <= scaled_pos[2] <= 0.55:
#            topmiddle_layer_oxygens.append(i)
#     print(f"找到{len(topmiddle_layer_oxygens)}个中间层氧原子")


    
#     # 找出距离小于cutoff的氢原子
#     oxygen_pos = atoms.positions[oxygen_indices[0]]
#     cutoff = 1.1
#     connected_hydrogens = []
#     for h_idx in hydrogen_indices:
#         h_pos = atoms.positions[h_idx]
#         distance = np.linalg.norm(h_pos - oxygen_pos)
#         if distance < cutoff:
#             connected_hydrogens.append(h_idx)
#     # 计算四面体顶点
#     tetrahedral_vertices = calculate_tetrahedral_vertices(oxygen_pos, 0.84632)
#     # 找出沿z轴正方向的氢原子
# endregion

def analyze_hydrogen_distribution(atoms):
    '''分析氢原子分布情况
    '''
    # 获取所有氧原子和氢原子索引
    oxygen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'O']
    hydrogen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'H']
    
    # 统计每个氧原子连接的氢原子数量
    o_h_counts = {}
    for o_idx in oxygen_indices:
        connected_hydrogens = find_hydrogen_connected_to_oxygen(atoms, o_idx)
        o_h_counts[o_idx] = len(connected_hydrogens)
    
    print("氧原子连接的氢原子统计:")
    print(f"氧原子总数: {len(oxygen_indices)}")
    print(f"氢原子总数: {len(hydrogen_indices)}")
    
    # 统计连接不同数量氢原子的氧原子数
    h_count_distribution = {}
    for count in o_h_counts.values():
        if count not in h_count_distribution:
            h_count_distribution[count] = 0
        h_count_distribution[count] += 1
    
    for count, num_oxygens in sorted(h_count_distribution.items()):
        print(f"连接{count}个氢原子的氧原子数: {num_oxygens}")

# region version1 def orient_top_middle_hydrogens
# def orient_top_middle_hydrogens(atoms):
#     '''找到蓝线标记的中间层水分子，调整其中朝向z轴正方向的H原子
#     将其调整到与同一水分子中另一个H原子相邻的四面体顶点位置
#     （确保不与原水分子中两个氢原子的位置重叠）
    
#     Args:
#         atoms: ASE原子对象
    
#     Returns:
#         修改后的ASE原子对象
#     '''
#     modified_atoms = atoms.copy()
#     # 识别O、H原子
#     hydrogen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'H']
#     oxygen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'O']
    
#     # 根据CIF文件找到蓝线所在中间层的氧原子（z约为0.43-0.56的氧原子）
#     middle_layer_oxygens = []
#     for i in oxygen_indices:
#         # 获取原子的分数坐标
#         scaled_pos = atoms.get_scaled_positions()[i]
#         if 0.43 <= scaled_pos[2] <= 0.56:
#             middle_layer_oxygens.append(i)
    
#     print(f"找到{len(middle_layer_oxygens)}个中间层氧原子")
    
#     # OH键长
#     oh_bond_length = 0.84632  # Å
    
#     # 记录修改的氢原子数量
#     modified_h_count = 0
    
#     # 处理每个中间层氧原子
#     for o_idx in middle_layer_oxygens:
#         # 获取氧原子位置
#         o_pos = atoms.positions[o_idx]
        
#         # 找出与该氧原子相连的氢原子
#         connected_hydrogens = find_hydrogen_connected_to_oxygen(atoms, o_idx)
        
#         # 如果找到正好两个氢原子（一个水分子）
#         if len(connected_hydrogens) == 2:
#             # 计算两个氢原子的位置
#             h1_pos = atoms.positions[connected_hydrogens[0]]
#             h2_pos = atoms.positions[connected_hydrogens[1]]
            
#             # 计算氢原子相对于氧原子的位移向量
#             h1_vec = h1_pos - o_pos
#             h2_vec = h2_pos - o_pos
            
#             # 判断哪个氢原子更靠近z轴正方向
#             z_positive_h_idx = None
#             other_h_pos = None
            
#             if h1_vec[2] > 0 and h1_vec[2] > h2_vec[2]:
#                 z_positive_h_idx = connected_hydrogens[0]
#                 other_h_pos = h2_pos
#             elif h2_vec[2] > 0 and h2_vec[2] > h1_vec[2]:
#                 z_positive_h_idx = connected_hydrogens[1]
#                 other_h_pos = h1_pos
            
#             # 如果找到了朝向z轴正方向的氢原子
#             if z_positive_h_idx is not None:
#                 # 使用蓝线圈出的四面体模型：一个顶点在z轴正方向，其余三个顶点在xy平面形成等边三角形
                
#                 # 第一个顶点在z轴正方向
#                 v1 = np.array([0, 0, 1])
                
#                 # 其余三个顶点在xy平面形成等边三角形，相隔120度
#                 v2 = np.array([np.sqrt(2)*2/3, 0, -1/3])
#                 v3 = np.array([-np.sqrt(2)/3, np.sqrt(6)/3, -1/3])
#                 v4 = np.array([-np.sqrt(2)/3, -np.sqrt(6)/3, -1/3])
                
#                 # 标准化向量
#                 vertices = np.array([v1, v2, v3, v4])
#                 vertices = vertices / np.linalg.norm(vertices, axis=1)[:, np.newaxis] * oh_bond_length
                
#                 # 将顶点移动到氧原子位置
#                 tetrahedral_vertices = vertices + o_pos
                
#                 # 获取当前两个氢原子的位置
#                 current_h_positions = np.array([
#                     atoms.positions[connected_hydrogens[0]],
#                     atoms.positions[connected_hydrogens[1]]
#                 ])
                
#                 # 先找出当前氢原子最接近哪个四面体顶点
#                 current_h_vertices = []
#                 for h_pos in current_h_positions:
#                     distances = [np.linalg.norm(h_pos - vertex) for vertex in tetrahedral_vertices]
#                     closest_vertex_idx = np.argmin(distances)
#                     current_h_vertices.append(closest_vertex_idx)
                
#                 # 找到一个未被占用的顶点
#                 available_vertices = list(set(range(4)) - set(current_h_vertices))
                
#                 # 如果没有可用顶点，则选择与当前位置最不接近的顶点
#                 if not available_vertices:
#                     # 找出与当前两个氢原子距离最远的顶点
#                     max_min_dist = -1
#                     best_vertex_idx = -1
#                     for i in range(4):
#                         min_dist = min(np.linalg.norm(tetrahedral_vertices[i] - current_h_positions[0]), 
#                                       np.linalg.norm(tetrahedral_vertices[i] - current_h_positions[1]))
#                         if min_dist > max_min_dist:
#                             max_min_dist = min_dist
#                             best_vertex_idx = i
                    
#                     best_vertex = tetrahedral_vertices[best_vertex_idx]
#                 else:
#                     # 从可用顶点中随机选择一个
#                     best_vertex_idx = np.random.choice(available_vertices)
#                     best_vertex = tetrahedral_vertices[best_vertex_idx]
                
#                 # 更新氢原子位置
#                 modified_atoms.positions[z_positive_h_idx] = best_vertex
#                 modified_h_count += 1
    
#     print(f"共修改了{modified_h_count}个朝向z轴正方向的氢原子位置")
#     return modified_atoms
# endregion

# region version2 def orient_top_middle_hydrogens
# def orient_top_middle_hydrogens(atoms):
    # '''找到蓝线标记的中间层水分子，调整其中朝向z轴正方向的H原子
    # 将其调整到与同一水分子中另一个H原子相邻的四面体顶点位置
    # （确保不与原水分子中两个氢原子的位置重叠，并且不使用v1顶点）
    
    # Args:
    #     atoms: ASE原子对象
    
    # Returns:
    #     修改后的ASE原子对象
    # '''
    # modified_atoms = atoms.copy()
    # # 识别O、H原子
    # hydrogen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'H']
    # oxygen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'O']
    
    # # 根据CIF文件找到蓝线所在中间层的氧原子（z约为0.43-0.56的氧原子）
    # middle_layer_oxygens = []
    # for i in oxygen_indices:
    #     # 获取原子的分数坐标
    #     scaled_pos = atoms.get_scaled_positions()[i]
    #     if 0.43 <= scaled_pos[2] <= 0.56:
    #         middle_layer_oxygens.append(i)
    
    # print(f"找到{len(middle_layer_oxygens)}个中间层氧原子")
    
    # # OH键长
    # oh_bond_length = 0.84632  # Å
    
    # # 记录修改的氢原子数量
    # modified_h_count = 0
    
    # # 处理每个中间层氧原子
    # for o_idx in middle_layer_oxygens:
    #     # 获取氧原子位置
    #     o_pos = atoms.positions[o_idx]
        
    #     # 找出与该氧原子相连的氢原子
    #     connected_hydrogens = find_hydrogen_connected_to_oxygen(atoms, o_idx)
        
    #     # 如果找到正好两个氢原子（一个水分子）
    #     if len(connected_hydrogens) == 2:
    #         # 计算两个氢原子的位置
    #         h1_pos = atoms.positions[connected_hydrogens[0]]
    #         h2_pos = atoms.positions[connected_hydrogens[1]]
            
    #         # 计算氢原子相对于氧原子的位移向量
    #         h1_vec = h1_pos - o_pos
    #         h2_vec = h2_pos - o_pos
            
    #         # 判断哪个氢原子更靠近z轴正方向
    #         z_positive_h_idx = None
    #         other_h_pos = None
            
    #         if h1_vec[2] > 0 and h1_vec[2] > h2_vec[2]:
    #             z_positive_h_idx = connected_hydrogens[0]
    #             other_h_pos = h2_pos
    #         elif h2_vec[2] > 0 and h2_vec[2] > h1_vec[2]:
    #             z_positive_h_idx = connected_hydrogens[1]
    #             other_h_pos = h1_pos
            
    #         # 如果找到了朝向z轴正方向的氢原子
    #         if z_positive_h_idx is not None:
    #             # 使用蓝线圈出的四面体模型：一个顶点在z轴正方向，其余三个顶点在xy平面形成等边三角形
                
    #             # 第一个顶点在z轴负方向
    #             v1 = np.array([0, 0, 1])
                
    #             # 其余三个顶点在xy平面形成等边三角形，相隔120度
    #             v2 = np.array([2*np.sqrt(2)/3, 0, -1/3])
    #             v3 = np.array([-np.sqrt(2)/3, np.sqrt(6)/3, -1/3])
    #             v4 = np.array([-np.sqrt(2)/3, -np.sqrt(6)/3, -1/3])
                
    #             # 标准化向量
    #             vertices = np.array([v1, v2, v3, v4])
    #             vertices = vertices / np.linalg.norm(vertices, axis=1)[:, np.newaxis] * oh_bond_length
                
    #             # 将顶点移动到氧原子位置
    #             tetrahedral_vertices = vertices + o_pos
                
    #             # 获取当前两个氢原子的位置
    #             current_h_positions = np.array([
    #                 atoms.positions[connected_hydrogens[0]],
    #                 atoms.positions[connected_hydrogens[1]]
    #             ])
                
    #             # 先找出当前氢原子最接近哪个四面体顶点
    #             current_h_vertices = []
    #             for h_pos in current_h_positions:
    #                 distances = [np.linalg.norm(h_pos - vertex) for vertex in tetrahedral_vertices]
    #                 closest_vertex_idx = np.argmin(distances)
    #                 current_h_vertices.append(closest_vertex_idx)
                
    #             # 找到未被占用且不是v1的顶点 (不包含索引0，因为0是v1顶点)
    #             available_vertices = list(set(range(1, 4)) - set(current_h_vertices))
                
    #             # 如果没有可用顶点（除v1外），则选择v2,v3,v4中与当前位置最不接近的顶点
    #             if not available_vertices:
    #                 # 只考虑顶点2,3,4 (索引1,2,3)
    #                 max_min_dist = -1
    #                 best_vertex_idx = -1
    #                 for i in range(1, 4):
    #                     min_dist = min(np.linalg.norm(tetrahedral_vertices[i] - current_h_positions[0]), 
    #                                   np.linalg.norm(tetrahedral_vertices[i] - current_h_positions[1]))
    #                     if min_dist > max_min_dist:
    #                         max_min_dist = min_dist
    #                         best_vertex_idx = i
                    
    #                 best_vertex = tetrahedral_vertices[best_vertex_idx]
    #             else:
    #                 # 从可用顶点中随机选择一个
    #                 best_vertex_idx = np.random.choice(available_vertices)
    #                 best_vertex = tetrahedral_vertices[best_vertex_idx]
                
    #             # 更新氢原子位置
    #             modified_atoms.positions[z_positive_h_idx] = best_vertex
    #             modified_h_count += 1
    
    # print(f"共修改了{modified_h_count}个朝向z轴正方向的氢原子位置")
    # return modified_atoms
# endregion

# region version3 def orient_top_middle_hydrogens
# def orient_top_middle_hydrogens(atoms, adjust_percentage=0.5):
#     '''找到蓝线标记的中间层水分子，调整其中朝向z轴正方向的H原子
#     将其中一定比例(默认50%)调整到与同一水分子中另一个H原子相邻的四面体顶点位置
#     （确保不与原水分子中两个氢原子的位置重叠，并且不使用v1顶点）
    
#     Args:
#         atoms: ASE原子对象
#         adjust_percentage: 需要调整的朝向z轴正方向的H原子的比例(0.0-1.0)
    
#     Returns:
#         修改后的ASE原子对象
#     '''
#     modified_atoms = atoms.copy()
#     # 识别O、H原子
#     hydrogen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'H']
#     oxygen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'O']
    
#     # 根据CIF文件找到蓝线所在中间层的氧原子（z约为0.43-0.56的氧原子）
#     middle_layer_oxygens = []
#     for i in oxygen_indices:
#         # 获取原子的分数坐标
#         scaled_pos = atoms.get_scaled_positions()[i]
#         if 0.43 <= scaled_pos[2] <= 0.56:
#             middle_layer_oxygens.append(i)
    
#     print(f"找到{len(middle_layer_oxygens)}个中间层氧原子")
    
#     # OH键长
#     oh_bond_length = 0.84632  # Å
    
#     # 找出所有符合条件的氧原子和对应的朝向z轴正方向的氢原子
#     candidate_adjustments = []
    
#     # 处理每个中间层氧原子，收集符合调整条件的氢原子
#     for o_idx in middle_layer_oxygens:
#         # 获取氧原子位置
#         o_pos = atoms.positions[o_idx]
        
#         # 找出与该氧原子相连的氢原子
#         connected_hydrogens = find_hydrogen_connected_to_oxygen(atoms, o_idx)
        
#         # 如果找到正好两个氢原子（一个水分子）
#         if len(connected_hydrogens) == 2:
#             # 计算两个氢原子的位置
#             h1_pos = atoms.positions[connected_hydrogens[0]]
#             h2_pos = atoms.positions[connected_hydrogens[1]]
            
#             # 计算氢原子相对于氧原子的位移向量
#             h1_vec = h1_pos - o_pos
#             h2_vec = h2_pos - o_pos
            
#             # 判断哪个氢原子更靠近z轴正方向
#             z_positive_h_idx = None
#             other_h_pos = None
            
#             if h1_vec[2] > 0 and h1_vec[2] > h2_vec[2]:
#                 z_positive_h_idx = connected_hydrogens[0]
#                 other_h_pos = h2_pos
#             elif h2_vec[2] > 0 and h2_vec[2] > h1_vec[2]:
#                 z_positive_h_idx = connected_hydrogens[1]
#                 other_h_pos = h1_pos
            
#             # 如果找到了朝向z轴正方向的氢原子，加入候选列表
#             if z_positive_h_idx is not None:
#                 candidate_adjustments.append((o_idx, z_positive_h_idx, other_h_pos))
    
#     # 随机选择一半的候选调整
#     num_to_adjust = int(len(candidate_adjustments) * adjust_percentage)
#     if num_to_adjust == 0 and len(candidate_adjustments) > 0:
#         num_to_adjust = 1  # 至少调整一个
    
#     # 随机选择要调整的氢原子
#     selected_adjustments = np.random.choice(
#         len(candidate_adjustments), 
#         size=min(num_to_adjust, len(candidate_adjustments)), 
#         replace=False
#     )
    
#     # 记录修改的氢原子数量
#     modified_h_count = 0
    
#     # 进行所选的调整
#     for idx in selected_adjustments:
#         o_idx, z_positive_h_idx, other_h_pos = candidate_adjustments[idx]
#         o_pos = atoms.positions[o_idx]
        
#         # 使用蓝线圈出的四面体模型：一个顶点在z轴正方向，其余三个顶点在xy平面形成等边三角形
        
#         # 第一个顶点在z轴负方向
#         v1 = np.array([0, 0, 1])
        
#         # 其余三个顶点在xy平面形成等边三角形，相隔120度
#         v2 = np.array([2*np.sqrt(2)/3, 0, -1/3])
#         v3 = np.array([-np.sqrt(2)/3, np.sqrt(6)/3, -1/3])
#         v4 = np.array([-np.sqrt(2)/3, -np.sqrt(6)/3, -1/3])
        
#         # 标准化向量
#         vertices = np.array([v1, v2, v3, v4])
#         vertices = vertices / np.linalg.norm(vertices, axis=1)[:, np.newaxis] * oh_bond_length
        
#         # 将顶点移动到氧原子位置
#         tetrahedral_vertices = vertices + o_pos
        
#         # 获取当前两个氢原子的位置
#         current_h_positions = np.array([
#             atoms.positions[z_positive_h_idx],
#             other_h_pos
#         ])
        
#         # 先找出当前氢原子最接近哪个四面体顶点
#         current_h_vertices = []
#         for h_pos in current_h_positions:
#             distances = [np.linalg.norm(h_pos - vertex) for vertex in tetrahedral_vertices]
#             closest_vertex_idx = np.argmin(distances)
#             current_h_vertices.append(closest_vertex_idx)
        
#         # 找到未被占用且不是v1的顶点 (不包含索引0，因为0是v1顶点)
#         available_vertices = list(set(range(1, 4)) - set(current_h_vertices))
        
#         # 如果没有可用顶点（除v1外），则选择v2,v3,v4中与当前位置最不接近的顶点
#         if not available_vertices:
#             # 只考虑顶点2,3,4 (索引1,2,3)
#             max_min_dist = -1
#             best_vertex_idx = -1
#             for i in range(1, 4):
#                 min_dist = min(np.linalg.norm(tetrahedral_vertices[i] - current_h_positions[0]), 
#                               np.linalg.norm(tetrahedral_vertices[i] - current_h_positions[1]))
#                 if min_dist > max_min_dist:
#                     max_min_dist = min_dist
#                     best_vertex_idx = i
            
#             best_vertex = tetrahedral_vertices[best_vertex_idx]
#         else:
#             # 从可用顶点中随机选择一个
#             best_vertex_idx = np.random.choice(available_vertices)
#             best_vertex = tetrahedral_vertices[best_vertex_idx]
        
#         # 更新氢原子位置
#         modified_atoms.positions[z_positive_h_idx] = best_vertex
#         modified_h_count += 1
    
#     print(f"总共找到{len(candidate_adjustments)}个朝向z轴正方向的氢原子，调整了其中{modified_h_count}个")
#     return modified_atoms
# endregion

# region version4 def orient_top_middle_hydrogens
# add the new feature which is orient special ratio of H atoms whose orientation is positive z direction
def orient_water_hydrogens_by_layers(atoms, adjust_percentage=0.5, layer_count=4, min_z=0.0, max_z=1.0):
    '''调整晶胞各层水分子中朝向z轴正方向的H原子
    
    将指定比例的朝向z轴正方向的氢原子调整到与同一水分子中另一个H原子相邻的四面体顶点位置
    (确保不与原水分子中两个氢原子的位置重叠，并且不使用z轴正方向顶点)
    
    Args:
        atoms: ASE原子对象
        adjust_percentage: 需要调整的朝向z轴正方向的H原子的比例(0.0-1.0)，默认为0.5
        layer_count: 将晶胞在z方向上划分的层数，默认为4
        min_z: 考虑的最小z分数坐标，默认为0.0
        max_z: 考虑的最大z分数坐标，默认为1.0
    
    Returns:
        修改后的ASE原子对象
    '''
    modified_atoms = atoms.copy()
    # 识别O、H原子
    hydrogen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'H']
    oxygen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == 'O']
    
    # 计算每层的厚度
    layer_thickness = (max_z - min_z) / layer_count
    
    # OH键长
    oh_bond_length = 0.84632  # Å
    
    # 总共找到的朝向z轴正方向的H原子数
    total_candidates = 0
    # 总共修改的H原子数
    total_modified = 0
    
    # 为每一层分别处理
    for layer in range(layer_count):
        z_min = min_z + layer * layer_thickness
        z_max = min_z + (layer + 1) * layer_thickness
        
        # 找出当前层的氧原子
        layer_oxygens = []
        for i in oxygen_indices:
            # 获取原子的分数坐标
            scaled_pos = atoms.get_scaled_positions()[i]
            if z_min <= scaled_pos[2] < z_max:
                layer_oxygens.append(i)
        
        print(f"第{layer+1}层 (z={z_min:.2f}-{z_max:.2f}) 找到{len(layer_oxygens)}个氧原子")
        
        # 找出所有符合条件的氧原子和对应的朝向z轴正方向的氢原子
        candidate_adjustments = []
        
        # 处理当前层的每个氧原子
        for o_idx in layer_oxygens:
            # 获取氧原子位置
            o_pos = atoms.positions[o_idx]
            
            # 找出与该氧原子相连的氢原子
            connected_hydrogens = find_hydrogen_connected_to_oxygen(atoms, o_idx)
            
            # 如果找到正好两个氢原子（一个水分子）
            if len(connected_hydrogens) == 2:
                # 计算两个氢原子的位置
                h1_pos = atoms.positions[connected_hydrogens[0]]
                h2_pos = atoms.positions[connected_hydrogens[1]]
                
                # 计算氢原子相对于氧原子的位移向量
                h1_vec = h1_pos - o_pos
                h2_vec = h2_pos - o_pos
                
                # 判断哪个氢原子更靠近z轴正方向
                z_positive_h_idx = None
                other_h_pos = None
                
                if h1_vec[2] > 0 and h1_vec[2] > h2_vec[2]:
                    z_positive_h_idx = connected_hydrogens[0]
                    other_h_pos = h2_pos
                elif h2_vec[2] > 0 and h2_vec[2] > h1_vec[2]:
                    z_positive_h_idx = connected_hydrogens[1]
                    other_h_pos = h1_pos
                
                # 如果找到了朝向z轴正方向的氢原子，加入候选列表
                if z_positive_h_idx is not None:
                    candidate_adjustments.append((o_idx, z_positive_h_idx, other_h_pos))
        
        total_candidates += len(candidate_adjustments)
        
        # 如果在这一层找到了候选调整
        if candidate_adjustments:
            # 随机选择指定比例的候选调整
            num_to_adjust = int(len(candidate_adjustments) * adjust_percentage)
            if num_to_adjust == 0 and len(candidate_adjustments) > 0:
                num_to_adjust = 1  # 至少调整一个
            
            # 随机选择要调整的氢原子
            selected_adjustments = np.random.choice(
                len(candidate_adjustments), 
                size=min(num_to_adjust, len(candidate_adjustments)), 
                replace=False
            )
            
            # 进行所选的调整
            modified_count = 0
            for idx in selected_adjustments:
                o_idx, z_positive_h_idx, other_h_pos = candidate_adjustments[idx]
                o_pos = atoms.positions[o_idx]
                
                # 使用四面体模型：一个顶点在z轴正方向，其余三个顶点在xy平面形成等边三角形
                
                # 第一个顶点在z轴方向
                v1 = np.array([0, 0, 1])
                
                # 其余三个顶点在xy平面形成等边三角形，相隔120度
                v2 = np.array([2*np.sqrt(2)/3, 0, -1/3])
                v3 = np.array([-np.sqrt(2)/3, np.sqrt(6)/3, -1/3])
                v4 = np.array([-np.sqrt(2)/3, -np.sqrt(6)/3, -1/3])
                
                # 标准化向量
                vertices = np.array([v1, v2, v3, v4])
                vertices = vertices / np.linalg.norm(vertices, axis=1)[:, np.newaxis] * oh_bond_length
                
                # 将顶点移动到氧原子位置
                tetrahedral_vertices = vertices + o_pos
                
                # 获取当前两个氢原子的位置
                current_h_positions = np.array([
                    atoms.positions[z_positive_h_idx],
                    other_h_pos
                ])
                
                # 先找出当前氢原子最接近哪个四面体顶点
                current_h_vertices = []
                for h_pos in current_h_positions:
                    distances = [np.linalg.norm(h_pos - vertex) for vertex in tetrahedral_vertices]
                    closest_vertex_idx = np.argmin(distances)
                    current_h_vertices.append(closest_vertex_idx)
                
                # 找到未被占用且不是v1的顶点 (不包含索引0，因为0是v1顶点)
                available_vertices = list(set(range(1, 4)) - set(current_h_vertices))
                
                # 如果没有可用顶点（除v1外），则选择v2,v3,v4中与当前位置最不接近的顶点
                if not available_vertices:
                    # 只考虑顶点2,3,4 (索引1,2,3)
                    max_min_dist = -1
                    best_vertex_idx = -1
                    for i in range(1, 4):
                        min_dist = min(np.linalg.norm(tetrahedral_vertices[i] - current_h_positions[0]), 
                                      np.linalg.norm(tetrahedral_vertices[i] - current_h_positions[1]))
                        if min_dist > max_min_dist:
                            max_min_dist = min_dist
                            best_vertex_idx = i
                    
                    best_vertex = tetrahedral_vertices[best_vertex_idx]
                else:
                    # 从可用顶点中随机选择一个
                    best_vertex_idx = np.random.choice(available_vertices)
                    best_vertex = tetrahedral_vertices[best_vertex_idx]
                
                # 更新氢原子位置
                modified_atoms.positions[z_positive_h_idx] = best_vertex
                modified_count += 1
            
            print(f"第{layer+1}层: 找到{len(candidate_adjustments)}个朝向z轴正方向的氢原子，调整了其中{modified_count}个")
            total_modified += modified_count
    
    print(f"总计: 找到{total_candidates}个朝向z轴正方向的氢原子，调整了其中{total_modified}个 (约{total_modified/total_candidates*100:.1f}%)")
    return modified_atoms

# endregion


## 运行示例
if __name__ == "__main__":
    # 设置输出目录
    output_dir = "/home/lihuimin/projects/BEofSpecialMoleculeOnIh/CO/out"
    # 搭建Ih unitcell
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
    randomized_Ihunitcell = Ihunitcell.copy()
    # scaled_positions = randomized_Ihunitcell.get_scaled_positions()
    # atoms_to_delete = [i for i, pos in enumerate(scaled_positions) if pos[2] < 0.3 or pos[2] > 0.7]
    # atoms_to_delete.sort(reverse=True)
    
    # # 删除范围外的原子
    # for idx in atoms_to_delete:
    #     del randomized_Ihunitcell[idx]
    # # 增大结构至四倍
    # randomized_Ihunitcell *= 4

    randomized_Ihunitcell = orient_water_hydrogens_by_layers(randomized_Ihunitcell, adjust_percentage=0.5, layer_count=4, min_z=0.0, max_z=1.0)
    # randomized_Ihunitcell = orient_toplayer_hydrogens(randomized_Ihunitcell)
    # print("初始结构分析:")
    # analyze_hydrogen_distribution(Ihunitcell)
    
    # # 随机化中间层氢原子
    # # randomized_Ihunitcell = randomize_middle_layer_hydrogen(Ihunitcell)
    # # randomized_Ihunitcell = orient_toplayer_hydrogens(randomized_Ihunitcell)
    # randomized_Ihunitcell = orient_toplayer_hydrogens(Ihunitcell)
    # print("\n随机化后结构分析:")
    # analyze_hydrogen_distribution(randomized_Ihunitcell)
    
    write(os.path.join(output_dir, 'Ih_middle_layer_randomized.cif'), randomized_Ihunitcell)
    