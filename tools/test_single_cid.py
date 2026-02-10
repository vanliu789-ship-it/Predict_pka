"""
测试单个CID的数据获取，用于调试
"""

import sys
import json
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from tools.pubchem_scraper_simple import SimplePubChemScraper

def test_cid(cid: int):
    """测试单个CID"""
    print("=" * 60)
    print(f"测试CID: {cid}")
    print("=" * 60)
    
    scraper = SimplePubChemScraper()
    
    print("\n正在获取数据...")
    compound = scraper.get_compound_data(cid)
    
    if compound:
        print("\n✅ 成功获取数据:")
        print("-" * 60)
        for key, value in compound.items():
            print(f"{key:20s}: {value}")
        print("-" * 60)
        return True
    else:
        print("\n❌ 未能获取数据")
        print("\n可能原因:")
        print("  1. 该化合物没有pKa数据")
        print("  2. pKa数据格式不标准")
        print("  3. 网络/API问题")
        
        # 尝试获取原始数据
        print("\n尝试获取原始数据...")
        url = f"{scraper.base_url}/compound/cid/{cid}/JSON"
        data = scraper._request(url)
        
        if data:
            print("✅ 成功获取原始数据（但没有提取到pKa）")
            
            # 保存到文件供检查
            debug_file = Path(f"data/raw/debug_cid_{cid}.json")
            debug_file.parent.mkdir(parents=True, exist_ok=True)
            with open(debug_file, 'w', encoding='utf-8') as f:
                json.dump(data, f, indent=2, ensure_ascii=False)
            print(f"原始数据已保存至: {debug_file}")
            print("你可以手动检查该文件中是否有pKa信息")
        else:
            print("❌ 无法获取原始数据")
        
        return False

def main():
    # 测试几个已知的化合物
    test_cids = [
        (2244, "阿司匹林"),
        (176, "乙酸"),
        (996, "苯酚"),
        (6115, "苯胺"),
        (5950, "丙氨酸"),
    ]
    
    print("测试多个已知化合物\n")
    
    results = []
    for cid, name in test_cids:
        print(f"\n{'='*60}")
        print(f"测试: {name} (CID: {cid})")
        print('='*60)
        
        scraper = SimplePubChemScraper()
        compound = scraper.get_compound_data(cid)
        
        if compound:
            print(f"✅ 成功: {compound['compound_name']}")
            print(f"   SMILES: {compound['smiles']}")
            print(f"   pKa: {compound['pka']}")
            results.append((name, True, compound['pka']))
        else:
            print(f"❌ 失败: 无法获取pKa数据")
            results.append((name, False, None))
        
        import time
        time.sleep(0.3)  # 避免请求过快
    
    # 总结
    print("\n" + "="*60)
    print("测试总结:")
    print("="*60)
    success_count = sum(1 for _, status, _ in results if status)
    print(f"成功: {success_count}/{len(results)}")
    print()
    
    for name, status, pka in results:
        icon = "✅" if status else "❌"
        pka_str = f"pKa={pka}" if pka else "无数据"
        print(f"{icon} {name:15s}: {pka_str}")
    
    if success_count == 0:
        print("\n⚠️  所有测试都失败了")
        print("\n可能的问题:")
        print("1. 网络连接问题（无法访问PubChem）")
        print("2. PubChem数据格式变化（需要更新提取逻辑）")
        print("3. 这些化合物的pKa数据不在预期位置")
        print("\n建议:")
        print("- 首先运行: python tools/diagnose.py")
        print("- 查看生成的 debug_cid_*.json 文件")
        print("- 考虑使用备用数据源")
    elif success_count < len(results):
        print(f"\n⚠️  部分测试失败 ({success_count}/{len(results)})")
        print("这是正常的，PubChem中不是所有化合物都有pKa数据")
    else:
        print("\n✅ 所有测试通过！工具工作正常")
        print("可以运行完整采集了")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--cid", type=int, help="测试单个CID")
    args = parser.parse_args()
    
    if args.cid:
        test_cid(args.cid)
    else:
        main()
