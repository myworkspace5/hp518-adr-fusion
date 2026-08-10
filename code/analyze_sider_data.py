# -*- coding: utf-8 -*-
import csv
import sys
sys.stdout.reconfigure(encoding='utf-8')

print("=== SIDER 数据分析：抗前列腺癌药物 ADR 探索 ===")
print()

sider_file = r"C:\Users\Bazinga\.qclaw\workspace\arptpred\data\sider\meddra_all_se.tsv"

# SIDER 格式（无表头，6列）：
# 列1: DrugBank ID (CID100...)
# 列2: PubChem CID (CID000...)
# 列3: UMLS CUI (C...)
# 列4: 术语类型 (PT=Preferred Term, LLT=Lower Level Term)
# 列5: MedDRA 概念 ID (C...)
# 列6: 副作用名称

print("正在读取 SIDER 数据...")
print()

drug_side_effects = {}
drug_names = {}
side_effect_counts = {}

total_rows = 0
with open(sider_file, 'r', encoding='utf-8') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        total_rows += 1
        if len(row) >= 6:
            drug_id = row[0]  # DrugBank ID
            side_effect = row[5]  # 副作用名称
            
            # 统计每个药物的副作用
            if drug_id not in drug_side_effects:
                drug_side_effects[drug_id] = []
            drug_side_effects[drug_id].append(side_effect)
            
            # 统计副作用出现频次
            if side_effect not in side_effect_counts:
                side_effect_counts[side_effect] = 0
            side_effect_counts[side_effect] += 1

print(f"数据读取完成！")
print(f"  总记录数: {total_rows:,}")
print(f"  唯一药物数: {len(drug_side_effects):,}")
print(f"  唯一副作用数: {len(side_effect_counts):,}")
print()

# 查找抗前列腺癌相关药物
print("=== 查找抗前列腺癌相关药物 ===")
print()

target_drugs = [
    "enzalutamide", "abiraterone", "apalutamide", "darolutamide",
    "docetaxel", "cabazitaxel", "prednisone", "degarelix",
    "leuprolide", "goserelin", "triptorelin", "buserelin",
    "flutamide", "bicalutamide", "nilutamide",
    "PSMA", "lutetium", "radium",
]

# 先读取药物名称映射表
drug_name_file = r"C:\Users\Bazinga\.qclaw\workspace\arptpred\data\sider\drug_names.tsv"
drug_name_map = {}
try:
    with open(drug_name_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) >= 2:
                drug_id = row[0]
                drug_name = row[1]
                drug_name_map[drug_id] = drug_name
    print(f"药物名称映射表已加载: {len(drug_name_map)} 个药物")
except Exception as e:
    print(f"药物名称映射表加载失败: {e}")

print()

# 在 SIDER 中搜索相关药物
found_drugs = []
for drug_id, drug_name in drug_name_map.items():
    for target in target_drugs:
        if target.lower() in drug_name.lower():
            found_drugs.append((drug_id, drug_name))
            break

print(f"找到 {len(found_drugs)} 个抗前列腺癌相关药物:")
print()
for drug_id, drug_name in found_drugs[:15]:  # 显示前15个
    se_count = len(drug_side_effects.get(drug_id, []))
    print(f"  {drug_id}: {drug_name} ({se_count} 个副作用)")

print()
print("=== 分析恩扎卢胺（Enzalutamide）ADR 谱 ===")
print()

# 重点分析恩扎卢胺（作为对照，与 AR-PROTAC 对比）
enzalutamide_ids = [did for did, dname in found_drugs if "enzalutamide" in dname.lower()]

if enzalutamide_ids:
    for eid in enzalutamide_ids:
        se_list = drug_side_effects.get(eid, [])
        print(f"药物: {drug_name_map.get(eid, eid)}")
        print(f"  DrugBank ID: {eid}")
        print(f"  副作用数量: {len(se_list)}")
        print()
        print("  前20个副作用:")
        for i, se in enumerate(se_list[:20]):
            print(f"    {i+1}. {se}")
        print()
        
        # 统计高频副作用
        from collections import Counter
        se_counter = Counter(se_list)
        print("  最常见的10个副作用:")
        for se, count in se_counter.most_common(10):
            print(f"    {se}: {count} 次报告")
        print()
else:
    print("未在 SIDER 中找到恩扎卢胺，尝试搜索 bicalutamide（比卡鲁胺）作为对照")
    bicalutamide_ids = [did for did, dname in found_drugs if "bicalutamide" in dname.lower()]
    if bicalutamide_ids:
        for bid in bicalutamide_ids[:1]:
            se_list = drug_side_effects.get(bid, [])
            print(f"药物: {drug_name_map.get(bid, bid)}")
            print(f"  副作用数量: {len(se_list)}")
            print()
            print("  前15个副作用:")
            for i, se in enumerate(se_list[:15]):
                print(f"    {i+1}. {se}")
            print()

print("=== 为 AR-PROTAC 分析做准备 ===")
print()
print("✓ SIDER 数据已加载，可用于:")
print("  1. 训练 ADR 预测模型的负样本（已知药物-副作用关联）")
print("  2. 作为知识图谱的种子数据")
print("  3. 与 AR-PROTAC 预测结果对比验证")
print()
print("下一步:")
print("  1. 手动从 PubChem 获取 AR-PROTAC 结构")
print("  2. 下载 DrugBank 5.1.1 XML 数据")
print("  3. 开始构建 AR-PROTAC 毒理知识图谱")
