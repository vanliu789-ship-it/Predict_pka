"""
快速测试脚本 - 采集少量数据验证工具可用性
"""

import sys
from pathlib import Path
import yaml

# 添加项目路径
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from tools.pubchem_scraper import PubChemScraper

def quick_test():
    """快速测试：采集20个化合物"""
    print("=" * 60)
    print("PubChem采集工具 - 快速测试")
    print("=" * 60)
    print("\n这将采集约20个化合物用于验证工具是否正常工作...")
    print("预计时间: 2-3分钟\n")
    
    # 加载并修改配置
    config_path = Path(__file__).parent / "config.yaml"
    with open(config_path, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    
    # 临时调整配置
    config['collection']['target_count'] = 20
    config['output']['output_file'] = 'pubchem_test.csv'
    
    # 保存临时配置
    temp_config_path = Path(__file__).parent / "config_test.yaml"
    with open(temp_config_path, 'w', encoding='utf-8') as f:
        yaml.dump(config, f)
    
    # 运行采集
    try:
        scraper = PubChemScraper(config_path=str(temp_config_path))
        scraper.run()
        
        print("\n" + "=" * 60)
        print("✅ 测试成功！工具运行正常")
        print(f"测试数据已保存至: {scraper.output_file}")
        print("\n现在可以运行完整采集:")
        print("  python tools/pubchem_scraper.py")
        print("=" * 60)
        
    except Exception as e:
        print(f"\n❌ 测试失败: {e}")
        print("\n请检查:")
        print("1. 网络连接是否正常")
        print("2. conda环境是否正确激活")
        print("3. 所有依赖是否已安装")
    
    finally:
        # 清理临时配置
        if temp_config_path.exists():
            temp_config_path.unlink()

if __name__ == "__main__":
    quick_test()
