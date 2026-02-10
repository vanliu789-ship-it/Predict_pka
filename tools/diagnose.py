"""
诊断工具 - 检查网络连接和API可用性
"""

import requests
import time
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

def test_network():
    """测试基本网络连接"""
    print("1. 测试网络连接...")
    try:
        response = requests.get("https://www.google.com", timeout=5)
        print("   ✅ 网络连接正常")
        return True
    except Exception as e:
        print(f"   ❌ 网络连接失败: {e}")
        return False

def test_pubchem_access():
    """测试PubChem访问"""
    print("\n2. 测试PubChem访问...")
    try:
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/property/MolecularFormula/JSON"
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        print(f"   ✅ PubChem API可访问")
        print(f"   响应时间: {response.elapsed.total_seconds():.2f}秒")
        return True
    except Exception as e:
        print(f"   ❌ PubChem访问失败: {e}")
        return False

def test_single_compound():
    """测试获取单个化合物数据"""
    print("\n3. 测试获取阿司匹林(CID: 2244)数据...")
    
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    
    # 3.1 测试基本属性
    print("   3.1 获取基本属性...")
    try:
        url = f"{base_url}/compound/cid/2244/property/MolecularFormula,MolecularWeight,CanonicalSMILES,Title/JSON"
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()
        
        if 'PropertyTable' in data:
            props = data['PropertyTable']['Properties'][0]
            print(f"   ✅ 化合物名称: {props.get('Title', 'N/A')}")
            print(f"   ✅ SMILES: {props.get('CanonicalSMILES', 'N/A')}")
            print(f"   ✅ 分子式: {props.get('MolecularFormula', 'N/A')}")
        else:
            print("   ⚠️  数据格式异常")
            return False
    except Exception as e:
        print(f"   ❌ 获取属性失败: {e}")
        return False
    
    # 3.2 测试完整记录（包含pKa）
    print("\n   3.2 获取完整记录（包含pKa数据）...")
    try:
        url = f"{base_url}/compound/cid/2244/JSON"
        response = requests.get(url, timeout=15)
        response.raise_for_status()
        data = response.json()
        
        if 'Record' not in data:
            print("   ⚠️  未找到Record字段")
            return False
        
        # 搜索pKa
        import re
        record = data['Record']
        sections = record.get('Section', [])
        
        pka_found = False
        full_text = str(sections)  # 简单搜索
        
        if 'pka' in full_text.lower() or 'dissociation' in full_text.lower():
            print("   ✅ 找到pKa相关信息")
            pka_found = True
            
            # 尝试提取数值
            matches = re.findall(r'pKa[^0-9]*([0-9]+\.?[0-9]*)', full_text, re.IGNORECASE)
            if matches:
                print(f"   📊 可能的pKa值: {matches[:3]}")
        else:
            print("   ⚠️  未找到pKa信息（这可能是正常的）")
        
        return True
        
    except Exception as e:
        print(f"   ❌ 获取完整记录失败: {e}")
        return False

def test_rdkit():
    """测试RDKit"""
    print("\n4. 测试RDKit...")
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles("CC(=O)O")
        if mol and mol.GetNumAtoms() > 0:
            print("   ✅ RDKit工作正常")
            return True
        else:
            print("   ❌ RDKit无法创建分子")
            return False
    except Exception as e:
        print(f"   ❌ RDKit导入失败: {e}")
        return False

def main():
    print("=" * 60)
    print("PubChem数据采集工具 - 诊断程序")
    print("=" * 60)
    print()
    
    results = []
    
    # 运行所有测试
    results.append(("网络连接", test_network()))
    time.sleep(0.5)
    
    results.append(("PubChem访问", test_pubchem_access()))
    time.sleep(0.5)
    
    results.append(("数据获取", test_single_compound()))
    time.sleep(0.5)
    
    results.append(("RDKit", test_rdkit()))
    
    # 总结
    print("\n" + "=" * 60)
    print("诊断总结:")
    print("=" * 60)
    
    all_pass = True
    for name, status in results:
        icon = "✅" if status else "❌"
        print(f"{icon} {name}: {'通过' if status else '失败'}")
        if not status:
            all_pass = False
    
    print("\n" + "=" * 60)
    
    if all_pass:
        print("✅ 所有测试通过！")
        print("\n可能的问题:")
        print("  1. PubChem数据库中某些化合物确实没有pKa数据")
        print("  2. pKa数据格式不标准，提取逻辑需要调整")
        print("\n建议:")
        print("  1. 尝试运行: python tools/test_single_cid.py")
        print("  2. 使用备用数据源（本地数据集）")
    else:
        print("❌ 存在问题")
        print("\n请根据上述失败项进行修复:")
        if not results[0][1]:
            print("  - 检查网络连接")
        if not results[1][1]:
            print("  - 检查防火墙/代理设置")
            print("  - 确认可以访问 https://pubchem.ncbi.nlm.nih.gov")
        if not results[2][1]:
            print("  - PubChem API可能暂时不可用")
        if not results[3][1]:
            print("  - 重新安装RDKit: conda install -c conda-forge rdkit")
    
    print("=" * 60)

if __name__ == "__main__":
    main()
