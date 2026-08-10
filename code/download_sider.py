# -*- coding: utf-8 -*-
import urllib.request
import sys
sys.stdout.reconfigure(encoding='utf-8')

print("=== 下载 SIDER 数据 ===")
print()

# SIDER URL
sider_url = "http://sideeffects.embl.de/media/download/meddra_all_se.tsv.gz"
out_file = r"C:\Users\Bazinga\.qclaw\workspace\arptpred\data\sider\meddra_all_se.tsv.gz"

print("正在下载: " + sider_url)
print("保存到: " + out_file)
print()

try:
    req = urllib.request.Request(sider_url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=30) as response:
        data = response.read()
        with open(out_file, 'wb') as f:
            f.write(data)
    print("下载完成！")
    
    # Check file size
    import os
    file_size = os.path.getsize(out_file)
    print("文件大小: " + str(round(file_size / 1024 / 1024, 2)) + " MB")
    
    # Decompress
    print()
    print("正在解压缩...")
    import gzip
    out_tsv = out_file.replace('.gz', '')
    with gzip.open(out_file, 'rb') as f_in:
        with open(out_tsv, 'wb') as f_out:
            f_out.write(f_in.read())
    print("解压完成: " + out_tsv)
    
    tsv_size = os.path.getsize(out_tsv)
    print("TSV 文件大小: " + str(round(tsv_size / 1024 / 1024, 2)) + " MB")
    
except Exception as e:
    print("下载失败: " + str(e))
    print()
    print("请手动下载:")
    print("  URL: " + sider_url)
    print("  保存为: " + out_file)

print()
print("---")
print()

# Download drug_names.tsv
drug_url = "http://sideeffects.embl.de/media/download/drug_names.tsv"
drug_file = r"C:\Users\Bazinga\.qclaw\workspace\arptpred\data\sider\drug_names.tsv"

print("正在下载药物名称映射...")
try:
    req2 = urllib.request.Request(drug_url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req2, timeout=30) as response:
        data2 = response.read()
        with open(drug_file, 'wb') as f:
            f.write(data2)
    print("下载完成: drug_names.tsv")
except Exception as e:
    print("下载失败: " + str(e))

print()
print("=== SIDER 数据下载完成 ===")
