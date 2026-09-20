# sheep-worldwide-wgs

**方法重建版，非原始代码。** 这是根据 Lv、Cao、Liu 等（2022）论文
[10.1093/molbev/msab353](https://doi.org/10.1093/molbev/msab353) 的公开方法重新编写的
分析工作流，不是论文作者当年使用脚本的归档，不保证直接重现论文结果。

**Method-based reconstruction, not the original study code.** The paper is
*Whole-Genome Resequencing of Worldwide Wild and Domestic Sheep Elucidates
Genetic Diversity, Introgression, and Agronomically Important Loci*,
Molecular Biology and Evolution 39(2):msab353.

## 包含什么 / Scope

- Snakemake：FASTQ 到 SNP/indel 的 WGS 流程、群体分析输入质控/LD pruning/IBS、
  π 与 FST 原始窗口扫描。ROH 需先解决论文参数问题；LD 需显式指定三样本子集。
- Python：PBS、ROD、fd 清理及显式假设下的 Z/BH 检验、分位数筛选、稳定的 ILS 概率；
  窗口对齐、样本抽样和输入文件哈希记录。
- 旧版工具封装：PSMC、SMC++、fastsimcoal、qpDstat、genomics_general。
  默认仅输出命令计划；显式执行时需要已有的有效输入与外部软件。
- 合成测试和 GitHub Actions 配置；不包含真实基因组、原始作者代码或第三方二进制。

SV 多工具合并、fastsimcoal 模型构建、SFS/祖先等位基因推断、sNMF、XP-CLR、SweeD、
RNA/miRNA 及富集等未能可靠自动化的部分作为明确的交接步骤，**不伪装成已实现模块**。
完整范围见 [方法与实现对照](docs/paper-methods.md) 和 [限制与参数冲突](docs/limitations.md)。

## 先运行合成示例

在仓库根目录执行，Linux/macOS，Python 3.11 或 3.12：

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-dev.txt
python -m pytest -q
snakemake --snakefile workflow/Snakefile --configfile config/config.test.yaml --dry-run --cores 1
snakemake --snakefile workflow/Snakefile --configfile config/config.test.yaml --cores 1
python scripts/check_demo.py
```

示例只使用合成小表，不需要下载任何测序数据或安装 GATK。结果位于 `results/demo/`，
包括 PBS、ROD、fd、ILS 和文件哈希记录。某些平台安装 Snakemake 的 `datrie` 依赖需要
C 编译器；也可使用 [environment.yml](environment.yml) 创建编排环境。

**测试通过仅证明这些代码检查和合成示例通过，不代表 810 只绵羊的全流程已复现。**

## 使用真实数据

1. 阅读 [workflow/README.md](workflow/README.md)，复制并填写
   [config/config.example.yaml](config/config.example.yaml)。其中 `null` 是有意保留的待确认项。
2. 按样本表与群体表模板准备文件、参考基因组、索引和坐标映射。
   大型 WGS 和参考数据不由本项目自动下载。
3. 仅启用需要的模块；使用 conda 环境或自行安装相应可执行程序。

```bash
cp config/config.example.yaml config/config.local.yaml
# 编辑 config/config.local.yaml，解决 null 参数并提供真实路径
snakemake --snakefile workflow/Snakefile --configfile config/config.local.yaml --dry-run --cores 8
snakemake --snakefile workflow/Snakefile --configfile config/config.local.yaml --use-conda --cores 8
```

真实运行可能需要大量内存、存储和计算时间。先在有代表性的小规模数据上验证。
历史论文版本与环境模板中可安装的版本不同之处必须记录，环境模板不等同于历史锁定环境。

## 辅助脚本示例

```bash
# ILS：保留 log(P)，避免把数值下溢误写为数学上的零
python scripts/window_stats.py ils --tract-bp 85400 --divergence-years 11000 \
  --generation-years 3 --recombination-rate 1.5e-8 --branch-factor 2 --output results/ils.tsv

# 先按坐标严格对齐 VCFtools π 窗口，再计算 ROD
python scripts/prepare_windows.py rod --landrace landrace.windowed.pi \
  --improved improved.windowed.pi --output results/rod_input.tsv
python scripts/window_stats.py rod --input results/rod_input.tsv --output results/rod.tsv

# 旧版工具：先查看帮助和非执行的 JSON 命令计划
python scripts/legacy_tools.py --help
```

更多输入/输出细节：[统计脚本](scripts/window_stats.md)、[旧版工具](docs/legacy-tools.md)、
[数据准备](docs/data-and-inputs.md)、[第三方来源](config/external_sources.json)。

## 来源、许可与引用

论文与补充材料没有提供可直接归档的完整作者 pipeline。本项目中所有新代码都明确属于
重建；外部脚本仅提供来源链接，不改变其版权或许可证。具体见 [NOTICE](NOTICE.md)。
新编写代码采用 MIT；请同时引用原论文和实际使用的外部工具。
[CITATION.cff](CITATION.cff) 保存论文完整作者信息和重建软件条目。
