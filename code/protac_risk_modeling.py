# -*- coding: utf-8 -*-
"""
AR-PROTAC ADR 预测 - Part 3: PROTAC 特异性风险建模
分析 E3 连接酶 (CRBN/VHL) 相关脱靶毒性 + PROTAC 机制特异性 ADR
"""
import json
import os
import sys
import time
import urllib.request
import urllib.parse
import ssl
from collections import defaultdict, Counter

sys.stdout.reconfigure(encoding='utf-8')

import paths as _paths
faers_dir = _paths.FAERS_DIR
output_dir = _paths.RESULTS_DIR
os.makedirs(output_dir, exist_ok=True)

ssl_ctx = ssl.create_default_context()

print("=" * 60)
print("Part 3: PROTAC 特异性风险建模")
print("=" * 60)
print()

# ============================================================
# Step 1: 获取已知 PROTAC 药物的 FAERS 数据
# 已上市/进入临床后期的 PROTAC:
# - ARV-110 (bavdegalutamide) - AR PROTAC
# - ARV-471 (vepdegestrant) - ER PROTAC  
# 这些可能 FAERS 中报告较少，但值得尝试
# 同时获取 IMiD 药物数据（CRBN 结合剂，作为 E3 连接酶毒性参照）
# ============================================================

protac_and_e3_drugs = {
    # 已知/疑似 PROTAC (FAERS 中可能有少量报告)
    "bavdegalutamide": {"category": "PROTAC-AR", "cn_name": "Bavdegalutamide(ARV-110)", "e3_ligase": "VHL", "target": "AR"},
    "vepdegestrant":     {"category": "PROTAC-ER", "cn_name": "Vepdegestrant(ARV-471)", "e3_ligase": "VHL", "target": "ER"},
    
    # CRBN 结合 IMiD 类药物（E3 连接酶毒性参照）
    "lenalidomide":      {"category": "IMiD-CRBN", "cn_name": "来那度胺", "e3_ligase": "CRBN", "target": "IKZF1/3"},
    "pomalidomide":      {"category": "IMiD-CRBN", "cn_name": "泊马度胺", "e3_ligase": "CRBN", "target": "IKZF1/3"},
    "thalidomide":       {"category": "IMiD-CRBN", "cn_name": "沙利度胺", "e3_ligase": "CRBN", "target": "unknown"},
    
    # VHL 抑制剂（E3 连接酶另一类配体）
    # 目前无已上市的 VHL 抑制剂，但可以用分子胶作为参照
    
    # 其他分子胶/靶向降解剂
    "ixazomib":          {"category": "Proteasome-Inhibitor", "cn_name": "伊沙佐米", "e3_ligase": "N/A", "target": "Proteasome"},
}

def fetch_faes_adrs(drug_name):
    base_url = "https://api.fda.gov/drug/event.json"
    params = urllib.parse.urlencode({
        "search": f'patient.drug.medicinalproduct:"{drug_name}"',
        "count": "patient.reaction.reactionmeddrapt.exact",
        "limit": "100"
    })
    url = f"{base_url}?{params}"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=30, context=ssl_ctx) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            results = data.get("results", [])
            return {r["term"]: r["count"] for r in results if r.get("count", 0) > 0}
    except Exception as e:
        return {}

def fetch_count(drug_name):
    base_url = "https://api.fda.gov/drug/event.json"
    params = urllib.parse.urlencode({
        "search": f'patient.drug.medicinalproduct:"{drug_name}"',
        "limit": "1"
    })
    url = f"{base_url}?{params}"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=30, context=ssl_ctx) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            return data.get("meta", {}).get("results", {}).get("total", 0)
    except:
        return 0

print("[Step 1] 获取 PROTAC/E3 连接酶相关药物的 FAERS 数据")
print("-" * 40)

e3_data = {}
for drug_name, info in protac_and_e3_drugs.items():
    print(f"  {info['cn_name']} ({drug_name})...", end=" ", flush=True)
    total = fetch_count(drug_name)
    adrs = fetch_faes_adrs(drug_name)
    e3_data[drug_name] = {**info, "total_reports": total, "adr_profile": adrs, "adr_count": len(adrs)}
    print(f"✓ {total} reports, {len(adrs)} ADRs")
    time.sleep(0.4)

print()

# ============================================================
# Step 2: E3 连接酶特异性毒性谱分析
# ============================================================
print("[Step 2] E3 连接酶特异性毒性谱分析")
print("-" * 40)

# CRBN 毒性谱（从 IMiD 数据提取）
crbn_drugs = {k: v for k, v in e3_data.items() if v["e3_ligase"] == "CRBN"}
crbn_profile = Counter()
for drug_name, data in crbn_drugs.items():
    for adr, count in data["adr_profile"].items():
        crbn_profile[adr] += count

print(f"\nCRBN 结合药物 (IMiD) 毒性谱 Top 20:")
print(f"  (这些是使用 CRBN E3 配体的 PROTAC 可能继承的风险)")
print()
crbn_top = crbn_profile.most_common(20)
for i, (adr, cnt) in enumerate(crbn_top):
    # 标注哪些药物贡献了这个 ADR
    sources = [d for d, data in crbn_drugs.items() if adr in data["adr_profile"]]
    print(f"  {i+1:2d}. {adr:<45} (总频次={cnt:>5}) [{', '.join(sources)}]")

# VHL 相关 — 从 ARV-110 等获取（如果有的话）
vhl_drugs = {k: v for k, v in e3_data.items() if v["e3_ligase"] == "VHL"}
vhl_profile = Counter()
for drug_name, data in vhl_drugs.items():
    for adr, count in data["adr_profile"].items():
        vhl_profile[adr] += count

if vhl_profile:
    print(f"\nVHL 配体 PROTAC 毒性谱:")
    for adr, cnt in vhl_profile.most_common(15):
        print(f"  - {adr}: {cnt}")
else:
    print("\nVHL 配体 PROTAC 在 FAERS 中暂无足够数据（ARV-110/471 尚在临床试验阶段）")

# ============================================================
# Step 3: PROTAC 特异性风险信号提取
# ============================================================
print()
print("[Step 3: PROTAC 特异性风险信号]")
print("-" * 40)

# HP518 使用 VHL 作为 E3 连接酶配体（基于 HINOVa 的公开信息推断）
# 如果是 CRBN 则用 CRBN 谱，这里两种都分析

# 已知的 PROTAC/分子胶特异性风险（来自文献）
protac_specific_risks = {
    "CRBN-related": {
        "teratogenicity": {
            "name": "Teratogenicity / Embryo-fetal toxicity",
            "cn_name": "致畸性 / 胚胎-胎儿毒性",
            "evidence": "沙利度胺历史教训; CRBN 介导的 SALL4/IKZF1 降解",
            "severity": "Critical",
            "mechanism": "CRBN-neosubstrate recruitment → SALL4/P63 degradation",
            "monitoring": "妊娠测试、严格避孕"
        },
        "neutropenia": {
            "name": "Neutropenia",
            "cn_name": "中性粒细胞减少症",
            "evidence": "IMiD 类药物常见剂量限制性毒性; IKZF1/3 降解影响造血",
            "severity": "High",
            "mechanism": "IKZF1/Ikaros family degradation → hematopoietic disruption",
            "monitoring": "定期 CBC 监测"
        },
        "venous_thromboembolism": {
            "name": "Venous Thromboembolism (VTE/DVT/PE)",
            "cn_name": "静脉血栓栓塞",
            "evidence": "来那度胺/泊马度胺黑框警告",
            "severity": "High",
            "mechanism": "多因素: TNFα变化、内皮细胞损伤、血小板激活",
            "monitoring": "VTE 风险评估、必要时预防性抗凝"
        },
        "secondary_malignancies": {
            "name": "Secondary Malignancies (SPM)",
            "cn_name": "第二原发恶性肿瘤",
            "evidence": "IMiD 长期使用增加 SPM 风险",
            "severity": "High",
            "mechanism": "免疫监视受损 + DNA损伤累积",
            "monitoring": "长期随访、皮肤检查"
        },
    },
    "VHL-related": {
        "hypertension": {
            "name": "Hypertension",
            "cn_name": "高血压",
            "evidence": "VHL-HIF通路参与血压调节; VHL抑制剂可能影响",
            "severity": "Moderate",
            "mechanism": "HIF-1α 稳定化 → 血管重构/水钠潴留",
            "monitoring": "定期血压监测"
        },
        "renal_effects": {
            "name": "Renal dysfunction / Electrolyte imbalance",
            "cn_name": "肾功能异常 / 电解质紊乱",
            "evidence": "VHL 参与肾细胞癌发生; VHL 抑制可能影响肾小管功能",
            "severity": "Moderate",
            "mechanism": "HIF 通路改变 → 肾脏氧感知和代谢调节异常",
            "monitoring": "肌酐/BUN/electrolytes 定期检测"
        },
        "iron_metabolism": {
            "name": "Iron metabolism dysregulation / Anemia",
            "cn_name": "铁代谢异常 / 贫血",
            "evidence": "VHL-HIF-EPO轴调控红细胞生成",
            "severity": "Moderate",
            "mechanism": "HIF-2α 稳定化 → EPO 过表达或失调",
            "monitoring": "CBC + iron studies"
        },
    },
    "PROTAC-general": {
        "hook_effect": {
            "name": "Hook effect / Loss of efficacy at high dose",
            "cn_name": "钩状效应 / 高浓度下失效",
            "evidence": "PROTAC 双结合机制在高浓度时形成二元复合物而非三元复合物",
            "severity": "Safety-relevant (dose escalation risk)",
            "mechanism": "[Protac-Target] + [Protac-E3] > [Protac-Target-E3]",
            "monitoring": "PK/PD 建模指导剂量选择"
        },
        "off_target_degradation": {
            "name": "Off-target protein degradation",
            "cn_name": "脱靶蛋白降解",
            "evidence": "PROTAC 可能通过 E3 连接酶降解非预期靶点",
            "severity": "Variable (depends on off-target)",
            "mechanism": "Neosubstrate recruitment by E3 ligase",
            "monitoring": "蛋白质组学筛选 + 广谱安全性监测"
        },
        "binary_complex_toxicity": {
            "name": "Binary complex-mediated toxicity",
            "cn_name": "二元复合物介导的毒性",
            "evidence": "未形成三元复合物的游离 PROTAC-靶点 或 PROTAC-E3 复合物可能产生独立毒性",
            "severity": "Unknown",
            "mechanism": "游离 warhead 或 E3 ligand 的脱靶效应",
            "monitoring": "与单独使用 warhead/E3-ligand 的安全性对比"
        }
    }
}

for category, risks in protac_specific_risks.items():
    cat_cn = {"CRBN-related": "CRBN 相关风险", "VHL-related": "VHL 相关风险", "PROTAC-general": "PROTAC 通用风险"}[category]
    print(f"\n{'='*50}")
    print(f"📌 {cat_cn}")
    print(f"{'='*50}")
    for risk_id, risk in risks.items():
        sev_emoji = {"Critical": "🔴🔴", "High": "🔴", "Moderate": "🟡", "Unknown/Safety-relevant": "🟠"}.get(risk["severity"], "⚪")
        print(f"\n  {sev_emoji} {risk['cn_name']} ({risk['name']})")
        print(f"     严重度: {risk['severity']}")
        print(f"     机制:   {risk['mechanism']}")
        print(f"     证据:   {risk['evidence']}")
        print(f"     监测:   {risk['monitoring']}")

# ============================================================
# Step 4: HP518 PROTAC 风险评分整合
# ============================================================
print()
print("\n[Step 4: HP518 整合风险评估矩阵]")
print("-" * 60)

# 加载之前的融合预测结果
fusion_report_path = os.path.join(output_dir, "hp518_faers_fusion_report.json")
with open(fusion_report_path, 'r', encoding='utf-8') as f:
    fusion_data = json.load(f)

# 构建 HP518 完整风险矩阵
hp518_risk_matrix = []

# A. 来自融合预测的高置信度 ADR (Top 30)
for pred in fusion_data.get("top50_predictions", [])[:30]:
    hp518_risk_matrix.append({
        "source": "Fusion_Prediction",
        "rank": pred["rank"],
        "adr": pred["adr"],
        "score": pred["fusion_score"],
        "confidence": "High" if pred["fusion_score"] > 80 else ("Medium" if pred["fusion_score"] > 60 else "Low"),
        "validated": pred.get("hp518_validated", False),
        "monitoring_priority": "P1" if pred["fusion_score"] > 90 else ("P2" if pred["fusion_score"] > 70 else "P3"),
    })

# B. PROTAC 特异性风险
for category, risks in protac_specific_risks.items():
    for risk_id, risk in risks.items():
        hp518_risk_matrix.append({
            "source": f"PROTAC_{category}",
            "adr": risk["name"],
            "cn_adr": risk["cn_name"],
            "score": None,
            "confidence": risk["severity"],
            "mechanism": risk["mechanism"],
            "evidence": risk["evidence"],
            "monitoring": risk["monitoring"],
            "monitoring_priority": "P1" if risk["severity"] in ["Critical", "High"] else ("P2" if risk["severity"] == "Moderate" else "P3"),
        })

# C. CRBN/VHL 特异性 ADR（从 FAERS 提取）
print("\n  E3 连接酶相关 FAERS 证据:")
for label, profile in [
    ("CRBN(IMiD)", crbn_profile),
]:
    if not profile:
        continue
    print(f"\n  --- {label} 特异性 ADR（可能被 HP518 继承）---")
    for adr, cnt in profile.most_common(10):
        # 检查是否已在融合预测 Top 50 中
        in_top50 = any(p["adr"].lower() == adr.lower() for p in fusion_data.get("top50_predictions", []))
        status = "✓已在融合预测中" if in_top50 else "⚡新增PROTAC特异风险"
        
        hp518_risk_matrix.append({
            "source": f"E3_FAE_{label}",
            "adr": adr,
            "faers_count": cnt,
            "status": status,
            "monitoring_priority": "P2" if cnt > 500 else "P3",
        })
        print(f"    {adr:<45} count={cnt:>6} {status}")

# 按 monitoring priority 排序
priority_order = {"P1": 0, "P2": 1, "P3": 2}
hp518_risk_matrix.sort(key=lambda x: priority_order.get(x.get("monitoring_priority", "P3"), 99))

# ============================================================
# Step 5: 输出最终整合报告
# ============================================================
print()
print("\n[Step 5: HP518 最终整合风险清单]")
print("=" * 70)
print(f"{'Pri':>3} | {'Source':<22} | {'ADR':<40} | {'Score':>6} | {'Confidence':<12}")
print("-" * 95)

for item in hp518_risk_matrix:
    pri = item.get("monitoring_priority", "P3")
    src = item.get("source", "")[:22]
    adr = item.get("adr", item.get("cn_adr", ""))[:40]
    score = str(item.get("score") or item.get("faers_count") or "-")
    conf = item.get("confidence", "-")[:12]
    print(f"{pri:>3} | {src:<22} | {adr:<40} | {score:>6} | {conf:<12}")

# 保存 PROTAC 特异性分析报告
protac_report = {
    "title": "HP518 PROTAC 特异性风险建模报告",
    "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
    
    "e3_ligase_data": {
        drug_name: {
            "cn_name": data["cn_name"],
            "category": data["category"],
            "e3_ligase": data["e3_ligase"],
            "total_reports": data["total_reports"],
            "unique_adrs": data["adr_count"]
        }
        for drug_name, data in e3_data.items()
    },
    
    "crbn_toxicity_profile": [
        {"adr": adr, "total_count": count}
        for adr, count in crbn_top
    ],
    
    "protac_specific_risks": protac_specific_risks,
    
    "hp518_integrated_risk_matrix": hp518_risk_matrix,
    
    "key_recommendations": [
        "HP518 若使用 VHL 配体：重点关注高血压、肾功能、贫血等 VHL-HIF 轴相关 AE",
        "HP518 若使用 CRBN 配体：必须执行妊娠测试 + VTE 风险评估 + 定期 CBC",
        "无论哪种 E3 配体：脱靶蛋白降解是最不可预测的风险，建议做蛋白质组学筛选",
        "Hook effect 可能在高剂量下导致疗效丧失和意外毒性，需 PK/PD 建模支持",
        "建议将本报告中 P1 项目纳入 Phase 2/3 的强制性安全性监测终点"
    ]
}

protac_report_path = os.path.join(output_dir, "hp518_protac_risk_model.json")
with open(protac_report_path, 'w', encoding='utf-8') as f:
    json.dump(protac_report, f, ensure_ascii=False, indent=2)

print(f"\n✓ PROTAC 风险模型报告: {protac_report_path}")
print("\n" + "=" * 60)
print("PROTAC 特异性风险建模完成!")
print("=" * 60)
