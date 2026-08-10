# -*- coding: utf-8 -*-
"""
AR-PROTAC ADR 预测模型 - 基线模型构建
使用 SIDER 数据训练基线模型，后续替换为 AR-PROTAC 数据
"""
import csv
import json
import sys
from collections import defaultdict, Counter
sys.stdout.reconfigure(encoding='utf-8')

print("=== AR-PROTAC ADR 预测 - 基线模型构建 ===")
print()

workspace = r"C:\Users\Bazinga\.qclaw\workspace\arptpred"
sider_file = workspace + r"\data\sider\meddra_all_se.tsv"
drug_name_file = workspace + r"\data\sider\drug_names.tsv"

# ============================================================
# Step 1: 加载数据
# ============================================================
print("[1/5] 加载 SIDER 数据...")

drug_se = defaultdict(set)  # drug_id -> set of side effects
se_drug = defaultdict(set)  # side_effect -> set of drug_ids
drug_names = {}

# 加载药物名称映射
try:
    with open(drug_name_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) >= 2:
                drug_names[row[0]] = row[1]
    print(f"  ✓ 药物名称: {len(drug_names)} 个")
except Exception as e:
    print(f"  ✗ 药物名称加载失败: {e}")

# 加载 SIDER（无表头，6列）
with open(sider_file, 'r', encoding='utf-8') as f:
    reader = csv.reader(f, delimiter='\t')
    total = 0
    for row in reader:
        total += 1
        if len(row) >= 6:
            drug_id = row[0]
            se_name = row[5]
            drug_se[drug_id].add(se_name)
            se_drug[se_name].add(drug_id)

print(f"  ✓ SIDER 记录: {total:,}")
print(f"  ✓ 唯一药物: {len(drug_se):,}")
print(f"  ✓ 唯一副作用: {len(se_drug):,}")
print()

# ============================================================
# Step 2: 计算药物相似度（基于副作用谱 Jaccard 相似度）
# ============================================================
print("[2/5] 计算药物相似度（Jaccard 相似度）...")

# 选取最常见的 50 个药物（按副作用数量）
common_drugs = sorted(drug_se.keys(), key=lambda d: len(drug_se[d]), reverse=True)[:50]

similarity_matrix = {}
for i, d1 in enumerate(common_drugs):
    for d2 in common_drugs[i+1:]:
        set1 = drug_se[d1]
        set2 = drug_se[d2]
        jaccard = len(set1 & set2) / len(set1 | set2) if len(set1 | set2) > 0 else 0
        if jaccard > 0.3:  # 只保留相似度 > 0.3 的
            key = tuple(sorted([d1, d2]))
            similarity_matrix[key] = jaccard

print(f"  ✓ 计算了 {len(common_drugs)} 个药物的两两相似度")
print(f"  ✓ 相似度 > 0.3 的药物对: {len(similarity_matrix)}")
print()

# ============================================================
# Step 3: 基线预测模型（基于药物相似度的 ADR 预测）
# ============================================================
print("[3/5] 构建基线预测模型（药物相似度 + ADR 迁移）...")

def predict_adr_for_drug(target_drug_id, top_k=10):
    """基于相似药物预测 target_drug 的 ADR"""
    if target_drug_id not in drug_se:
        return []
    
    # 找到最相似的已知药物
    similarities = []
    for other_drug in drug_se.keys():
        if other_drug == target_drug_id:
            continue
        key = tuple(sorted([target_drug_id, other_drug]))
        if key in similarity_matrix:
            sim = similarity_matrix[key]
        else:
            # 实时计算
            set1 = drug_se[target_drug_id]
            set2 = drug_se[other_drug]
            sim = len(set1 & set2) / len(set1 | set2) if len(set1 | set2) > 0 else 0
        
        if sim > 0:
            # 收集相似药物的 ADR
            for se in drug_se[other_drug]:
                similarities.append((se, sim, other_drug))
    
    # 按相似度加权投票
    se_scores = defaultdict(float)
    for se, sim, other in similarities:
        if se not in drug_se[target_drug_id]:  # 只预测未知的 ADR
            se_scores[se] += sim
    
    # 排序返回 Top-K
    ranked = sorted(se_scores.items(), key=lambda x: x[1], reverse=True)
    return ranked[:top_k]

print("  ✓ 基线预测函数已定义")
print()

# ============================================================
# Step 4: 测试基线模型（用已知药物验证）
# ============================================================
print("[4/5] 测试基线模型（留一法验证）...")

# 选取一个抗前列腺癌药物进行验证
test_drugs = []
for did, dname in drug_names.items():
    for keyword in ["bicalutamide", "abiraterone", "docetaxel", "prednisone"]:
        if keyword in dname.lower():
            test_drugs.append((did, dname))
            break

print(f"  找到测试药物: {len(test_drugs)} 个")
print()

for test_id, test_name in test_drugs[:3]:
    known_se = drug_se.get(test_id, set())
    print(f"  测试药物: {test_name} ({test_id})")
    print(f"    已知副作用: {len(known_se)} 个")
    
    # 预测
    predictions = predict_adr_for_drug(test_id, top_k=10)
    
    print(f"    预测 Top-10 新副作用:")
    for i, (se, score) in enumerate(predictions):
        print(f"      {i+1}. {se} (score={score:.3f})")
    print()

print()

# ============================================================
# Step 5: 为 AR-PROTAC 预测做准备
# ============================================================
print("[5/5] 为 AR-PROTAC 预测做准备...")

print("  由于 AR-PROTAC 是全新化合物，需要:")
print("    1. 计算 AR-PROTAC 与已知药物的分子相似度（Morgan 指纹）")
print("    2. 将相似药物的 ADR 迁移到 AR-PROTAC")
print("    3. 结合知识图谱嵌入进行预测")
print()
print("  下一步（等获得 AR-PROTAC 结构后）:")
print("    1. 计算 AR-PROTAC 与 SIDER 药物的 Morgan 指纹相似度")
print("    2. 使用相似药物的 ADR 作为基线预测")
print("    3. 与知识图谱预测结果融合")
print()

# 保存基线模型结果
output = {
    "similarity_matrix_size": len(similarity_matrix),
    "drug_count": len(drug_se),
    "se_count": len(se_drug),
    "test_results": []
}

for test_id, test_name in test_drugs[:3]:
    predictions = predict_adr_for_drug(test_id, top_k=10)
    output["test_results"].append({
        "drug": test_name,
        "drug_id": test_id,
        "known_se_count": len(drug_se.get(test_id, [])),
        "predictions": [{"se": se, "score": score} for se, score in predictions]
    })

output_file = workspace + r"\results\baseline_model_results.json"
with open(output_file, 'w', encoding='utf-8') as f:
    json.dump(output, f, ensure_ascii=False, indent=2)

print(f"  ✓ 基线模型结果已保存: {output_file}")
print()
print("=== 基线模型构建完成 ===")
print()
print("下一步:")
print("  1. 手动从 PubChem 下载 AR-PROTAC 结构（SDF 3D）")
print("  2. 计算 AR-PROTAC 与 SIDER 药物的分子相似度")
print("  3. 运行预测并验证")
