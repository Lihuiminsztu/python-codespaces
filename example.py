import platform
import numpy as np
import os

def get_cpu_info():
    """获取CPU信息"""
    print("Python版本:", platform.python_version())
    print("处理器架构:", platform.processor())
    print("系统:", platform.system())
    print("CPU信息:")
    
    try:
        # Linux下获取CPU信息
        with open('/proc/cpuinfo', 'r') as f:
            cpu_info = f.read()
            
        # 提取CPU型号信息
        for line in cpu_info.split('\n'):
            if 'model name' in line:
                print(line.split(':')[1].strip())
                break
    except:
        print("无法获取详细CPU信息")
    
    # 尝试导入ASE和tblite
    try:
        import ase
        print("\nASE版本:", ase.__version__)
    except ImportError:
        print("\nASE未安装")
    
    try:
        import tblite
        print("tblite版本:", tblite.__version__)
    except ImportError:
        print("tblite未安装")
    
    # NumPy测试
    print("\n执行NumPy计算测试:")
    arr = np.random.rand(1000, 1000)
    result = np.dot(arr, arr.T)
    print(f"NumPy数组形状: {result.shape}")

if __name__ == "__main__":
    get_cpu_info()
