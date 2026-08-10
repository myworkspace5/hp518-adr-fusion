# AR-PROTAC ADR 预测 - 数据下载脚本
# 需要 PowerShell 5.1+ / PowerShell Core

Write-Host "=== AR-PROTAC ADR 预测 - 数据下载 ===" -ForegroundColor Cyan
Write-Host ""

# 设置工作目录
$workspace = "C:\Users\Bazinga\.qclaw\workspace\arptpred"
$dataDir = Join-Path $workspace "data"

# 创建目录结构
$directories = @(
    "data",
    "data\sider",
    "data\drugbank",
    "data\pubchem",
    "data\pdb",
    "data\alphafold",
    "structures",
    "models",
    "results",
    "results\offtarget",
    "results\knowledge_graph",
    "results\predictions"
)

foreach ($dir in $directories) {
    $fullPath = Join-Path $workspace $dir
    if (-not (Test-Path $fullPath)) {
        New-Item -ItemType Directory -Force -Path $fullPath | Out-Null
        Write-Host "  创建: $dir" -ForegroundColor Gray
    }
}

Write-Host ""
Write-Host "目录结构已创建: $workspace" -ForegroundColor Green
Write-Host ""

# ============================================================
# 1. 下载 SIDER 数据（不需要注册）
# ============================================================
Write-Host "--- 1. 下载 SIDER 数据 ---" -ForegroundColor Yellow

$siderUrl = "http://sideeffects.embl.de/media/download/meddra_all_se.tsv.gz"
$siderGz = Join-Path $dataDir "sider\meddra_all_se.tsv.gz"
$siderTsv = Join-Path $dataDir "sider\meddra_all_se.tsv"

try {
    Write-Host "  正在下载 SIDER 数据..." -ForegroundColor Gray
    Invoke-WebRequest -Uri $siderUrl -OutFile $siderGz -UserAgent "Mozilla/5.0"
    Write-Host "  下载完成: meddra_all_se.tsv.gz" -ForegroundColor Green
    
    # 解压缩
    Write-Host "  正在解压缩..." -ForegroundColor Gray
    $gzipStream = [System.IO.Compression.GzipStream]::new(
        [System.IO.File]::OpenRead($siderGz),
        [System.IO.Compression.CompressionMode]::Decompress
    )
    $outFile = [System.IO.File]::Create($siderTsv)
    $gzipStream.CopyTo($outFile)
    $gzipStream.Close()
    $outFile.Close()
    Write-Host "  解压完成: meddra_all_se.tsv" -ForegroundColor Green
    
    # 显示文件信息
    $fileInfo = Get-Item $siderTsv
    Write-Host "  文件大小: $([math]::Round($fileInfo.Length/1MB, 2)) MB" -ForegroundColor Cyan
}
catch {
    Write-Host "  SIDER 下载失败: $_" -ForegroundColor Red
    Write-Host "  请手动下载: http://sideeffects.embl.de/media/download/meddra_all_se.tsv.gz" -ForegroundColor Yellow
}

Write-Host ""

# ============================================================
# 2. 下载 SIDER 药物名称映射表
# ============================================================
Write-Host "--- 2. 下载 SIDER 药物名称映射 ---" -ForegroundColor Yellow

$siderDrugUrl = "http://sideeffects.embl.de/media/download/drug_names.tsv"
$siderDrugFile = Join-Path $dataDir "sider\drug_names.tsv"

try {
    Write-Host "  正在下载药物名称映射..." -ForegroundColor Gray
    Invoke-WebRequest -Uri $siderDrugUrl -OutFile $siderDrugFile -UserAgent "Mozilla/5.0"
    Write-Host "  下载完成: drug_names.tsv" -ForegroundColor Green
}
catch {
    Write-Host "  药物名称映射下载失败: $_" -ForegroundColor Red
}

Write-Host ""

# ============================================================
# 3. 下载 OFFSIDES 数据（FDA 不良事件）
# ============================================================
Write-Host "--- 3. 下载 OFFSIDES 数据（FDA 不良事件）---" -ForegroundColor Yellow

$offSidesUrl = "http://www.bi.csic.es/~lespinosa/safety/offsides.tsv.bz2"
$offSidesFile = Join-Path $dataDir "sider\offsides.tsv.bz2"

Write-Host "  注意: OFFSIDES 数据较大 (~2GB)，可选下载" -ForegroundColor Yellow
$downloadOffsides = Read-Host "  是否下载 OFFSIDES? (y/n, 默认 n)"

if ($downloadOffsides -eq "y") {
    try {
        Write-Host "  正在下载 OFFSIDES 数据..." -ForegroundColor Gray
        Invoke-WebRequest -Uri $offSidesUrl -OutFile $offSidesFile -UserAgent "Mozilla/5.0"
        Write-Host "  下载完成: offsides.tsv.bz2" -ForegroundColor Green
    }
    catch {
        Write-Host "  OFFSIDES 下载失败: $_" -ForegroundColor Red
    }
}
else {
    Write-Host "  跳过 OFFSIDES 下载" -ForegroundColor Gray
}

Write-Host ""

# ============================================================
# 4. 获取 DrugBank 数据（需要登录后手动下载）
# ============================================================
Write-Host "--- 4. DrugBank 数据 ---" -ForegroundColor Yellow
Write-Host "  您已登录 DrugBank，请手动下载以下文件：" -ForegroundColor Cyan
Write-Host "  1. 访问: https://go.drugbank.com/releases/5-1-1/downloads/all-full-database"
Write-Host "  2. 下载: all_drugbank_v5.1.1.xml.zip"
Write-Host "  3. 解压后放置到: $dataDir\drugbank\"
Write-Host ""

# ============================================================
# 5. 搜索并下载 AR-PROTAC 结构（PubChem）
# ============================================================
Write-Host "--- 5. 搜索 AR-PROTAC 化学结构 ---" -ForegroundColor Yellow

$protacs = @(
    @{Name="ARV-110"; PubChemCIDs="155885553"},  # ARV-110 的 PubChem CID
    @{Name="HP518"; PubChemCIDs="165312872"},    # HP518 的 PubChem CID（如有）
    @{Name="Bardoxolone"; PubChemCIDs="24827720"} # 对照：Bardoxolone methyl
)

foreach ($protac in $protacs) {
    $name = $protac.Name
    $cid = $protac.PubChemCIDs
    $outFile = Join-Path $dataDir "pubchem\$name`_$cid.sdf"
    
    Write-Host "  正在获取 $name (CID: $cid)..." -ForegroundColor Gray
    
    # PubChem REST API - 获取 SDF 3D 结构
    $pubchemUrl = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/$cid/record/SDF/?record_type=3d"
    
    try {
        Invoke-WebRequest -Uri $pubchemUrl -OutFile $outFile -UserAgent "Mozilla/5.0"
        Write-Host "  已保存: $name`_$cid.sdf" -ForegroundColor Green
    }
    catch {
        Write-Host "  $name 结构获取失败（CID $cid 可能不正确）" -ForegroundColor Yellow
        Write-Host "  请手动在 PubChem 搜索: https://pubchem.ncbi.nlm.nih.gov/" -ForegroundColor Gray
    }
}

Write-Host ""

# ============================================================
# 6. 下载 PDB 蛋白质结构
# ============================================================
Write-Host "--- 6. 下载 PDB 蛋白质结构 ---" -ForegroundColor Yellow

$pdbs = @(
    @{PDB="5VA1"; Desc="hERG 通道"},
    @{PDB="6MAH"; Desc="CYP3A4"},
    @{PDB="1M2X"; Desc="糖皮质激素受体"},
    @{PDB="4CI4"; Desc="Cereblon + 沙利度胺"},
    @{PDB="2AM9"; Desc="AR LBD"}
)

foreach ($pdb in $pdbs) {
    $pdbId = $pdb.PDB
    $desc = $pdb.Desc
    $outFile = Join-Path $dataDir "pdb\$pdbId.pdb"
    
    Write-Host "  正在下载 $pdbId ($desc)..." -ForegroundColor Gray
    $pdbUrl = "https://files.rcsb.org/download/$pdbId.pdb"
    
    try {
        Invoke-WebRequest -Uri $pdbUrl -OutFile $outFile -UserAgent "Mozilla/5.0"
        Write-Host "  已保存: $pdbId.pdb" -ForegroundColor Green
    }
    catch {
        Write-Host "  $pdbId 下载失败" -ForegroundColor Yellow
    }
}

Write-Host ""

# ============================================================
# 7. 下载 AlphaFold 蛋白质结构
# ============================================================
Write-Host "--- 7. 下载 AlphaFold 蛋白质结构 ---" -ForegroundColor Yellow

$alphafoldProteins = @(
    @{Uniprot="Q9UKD9"; Desc="IKZF1"},
    @{Uniprot="Q99589"; Desc="IKZF3"},
    @{Uniprot="P16638"; Desc="CD36"},
    @{Uniprot="P10275"; Desc="AR NTD"}
)

foreach ($protein in $alphafoldProteins) {
    $uniprot = $protein.Uniprot
    $desc = $protein.Desc
    $outFile = Join-Path $dataDir "alphafold\$uniprot.pdb"
    
    Write-Host "  正在下载 AlphaFold $uniprot ($desc)..." -ForegroundColor Gray
    $afUrl = "https://alphafold.ebi.ac.uk/files/AF-$uniprot-F1-model_v4.pdb"
    
    try {
        Invoke-WebRequest -Uri $afUrl -OutFile $outFile -UserAgent "Mozilla/5.0"
        Write-Host "  已保存: AF-$uniprot-F1-model_v4.pdb" -ForegroundColor Green
    }
    catch {
        Write-Host "  AlphaFold $uniprot 下载失败" -ForegroundColor Yellow
    }
}

Write-Host ""

# ============================================================
# 完成
# ============================================================
Write-Host "=============================================" -ForegroundColor Cyan
Write-Host "数据下载完成！" -ForegroundColor Green
Write-Host ""
Write-Host "下一步：" -ForegroundColor Cyan
Write-Host "1. 手动下载 DrugBank 5.1.1 XML 文件"
Write-Host "   -> 放置到: $dataDir\drugbank\"
Write-Host "2. 在 PubChem 中搜索 ARV-110, ARV-766, HP518"
Write-Host "   -> 下载 SDF 3D 结构"
Write-Host "   -> 放置到: $dataDir\pubchem\"
Write-Host "3. 运行 Python 环境设置脚本"
Write-Host ""
Write-Host "工作目录: $workspace" -ForegroundColor Yellow
Write-Host "=============================================" -ForegroundColor Cyan
