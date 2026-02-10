"""
快速测试脚本 - 采集少量数据验证工具可用性
使用简化版采集器，更快速可靠
"""

import sys
from pathlib import Path

# 添加项目路径
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from tools.pubchem_scraper_simple import SimplePubChemScraper

def quick_test():
    """快速测试：采集20个化合物"""
    print("=" * 60)
    print("PubChem采集工具 - 快速测试")
    print("=" * 60)
    print("\n使用简化版采集器，从已知化合物列表采集...")
    print("这将采集20个化合物用于验证工具是否正常工作")
    print("预计时间: 1-2分钟\n")
    
    # 运行采集
    try:
        # 创建临时输出目录
        temp_output = "data/raw/test"
        scraper = SimplePubChemScraper(output_dir=temp_output)
        
        # 采集20个化合物
        df = scraper.collect(target_count=20)
        
        if len(df) > 0:
            print("\n" + "=" * 60)
            print("✅ 测试成功！工具运行正常")
            print(f"\n📊 采集结果:")
            print(f"  - 成功采集: {len(df)} 个化合物")
            print(f"  - pKa范围: {df['pka'].min():.2f} ~ {df['pka'].max():.2f}")
            print(f"  - 数据文件: {scraper.output_file}")
            print("\n现在可以运行完整采集:")
            print("  python tools/pubchem_scraper_simple.py --target 1000")
            print("=" * 60)
        else:
            print("\n⚠️  未采集到数据，但工具可以正常运行")
            print("可能原因：网络问题或PubChem服务暂时不可用")
        
    except Exception as e:
        print(f"\n❌ 测试失败: {e}")
        print("\n请检查:")
        print("1. 网络连接是否正常（需访问PubChem）")
        print("2. conda环境是否正确激活: conda activate pka_predictor")
        print("3. 所有依赖是否已安装")
        import traceback
        print("\n详细错误信息:")
        traceback.print_exc()

if __name__ == "__main__":
    quick_test()
