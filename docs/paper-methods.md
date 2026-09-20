# Paper-to-code scope

Reference: Lv F-H, Cao Y-H, Liu G-J, et al. (2022), *Whole-Genome
Resequencing of Worldwide Wild and Domestic Sheep Elucidates Genetic
Diversity, Introgression, and Agronomically Important Loci*, MBE 39(2):
msab353. DOI: [10.1093/molbev/msab353](https://doi.org/10.1093/molbev/msab353).

The names below indicate implementation boundaries, not an assertion of a
fully automated reproduction of the paper. See the workflow and wrapper
guides for actual inputs and supported targets.

| Analysis | Paper settings | Delivery |
| --- | --- | --- |
| Read QC / SNP / indel calling | FastQC 0.11.9; Trimmomatic 0.36; BWA 0.7.8; SAMtools 1.9; GATK 4.1.2.0; depth >10 | Snakemake WGS module, explicit user parameters and depth gate |
| SNP hard filters | QUAL 30; QD 2; MQ 40; FS 60; SOR 3; MQRankSum −12.5; ReadPosRankSum −8 | GATK filtering and PASS biallelic SNP extraction |
| Indel hard filters | QUAL 30; QD 2; FS 200; ReadPosRankSum −20 | GATK filtering and PASS indel extraction |
| Variant annotation | SnpEff 4.3t, Oar_rambouillet_v1.0 | External annotation handoff; requires database build and reference matching |
| SV calling | >15×; Manta 1.3.2, DELLY 0.8.3, smoove 0.2.6; 50 bp–1 Mb, multi-caller agreement | Documented handoff; consensus rules underspecified |
| Unrelated sample selection | KING 2.2.5; remove one in pairs with kinship >0.0884 | User-curated unrelated keep file; relatedness decision not invented |
| Population structure input | MAF ≥0.05; HWE ≥0.001; missing ≤0.1; prune 50/5/0.2 | Snakemake QC/pruning and IBS; HWE pooling caveat applies |
| NJ / PCA / sNMF | SplitsTree5 5.1.7-beta; EIGENSOFT 6.0.1; sNMF 1.2 K=1–20 | smartpca wrapper on prepared EIGENSTRAT; NJ/sNMF manual handoff |
| Diversity | π/He 1 Mb; VCFtools 0.1.17 | Per-individual π plus supplementary pooled π; He estimator is not reconstructed as `--het` (which is individual inbreeding) |
| ROH | PLINK1.90p; paper's four ROH flags | Gated workflow; unresolved window-kb flag documented |
| LD decay | 3 randomly chosen individuals; PopLDdecay3.41; maximum300kb | Sample selection helper; optional PopLDdecay workflow with supplied three-member lists |
| Pairwise/global FST | GCTA1.93.2beta on separate220350-SNP set; VCFtools100kb | 100kb diversity FST/GCTA manual handoff; selection FST module is separate |
| PSMC | 3 highest-depth per population; -N25 -t15 -r5; 4+25*2+4+6 | Safe legacy command wrapper; prepared consensus required |
| SMC++ | 1.14.0.dev0;7–10 unrelated high-depth per region | Wrapper on prepared masked SMC input |
| SNeP |1.11; no missing;MAF>0.05 | Manual legacy-tool handoff |
| fastsimcoal2.6 |1e6 simulations;25/65 cycles;100 optimizations;100 block bootstraps ×20 runs | Single-run command wrapper; original models/bootstraps not bundled |
| D-statistics | ADMIXTOOLS7.0.1 qpDstat; absZ>3 | Gated parameter-file wrapper; historical chromosome handling unresolved; recipient/model curation manual |
| fd / π / dxy |100kb,min500SNPs;20/50kb step conflict | External genomics_general wrappers with explicit step; new fd postprocessor |
| ILS |Gamma shape2;r1.5cM/Mb;generation3y | Tested helper with log survival |
| Improvement FST |50kb;10/25kb step conflict;top1% | Snakemake raw scans + explicit percentile helper |
| ROD |1−πimproved/πlandrace;10kb nonoverlap | π scans, strict window join, tested ROD helper |
| Fleece XP-CLR |v1.0;-w1 .005 200 2000 2 -p0 .95;top1% | Manual legacy CLI handoff; generic quantile helper on supplied scores; region aggregation must be specified |
| SV selection |200 permutations;p<.01 | Manual handoff; exchangeability and permutation units must be defined |
| Dwarf SweeD |Ouessant;10kb;top5% | Manual handoff; SweeD grid is not automatically a 10kb sliding window |
| Altitude PBS |CHA/EFR/HUS;99th percentile | New tested transformation of supplied aligned pairwise FST windows, not original R estimator |
| GO/KEGG |DAVID6.8;≥6genes;Bonferroni p<.05 | External handoff; save background and database version |
| RNA |FastQC,Cutadapt,TopHat,Cufflinks;IRF2BP2 chr25:7067974–7071785 | Method/inputs handoff, not automated |
| miRNA |Cutadapt18–30nt;TargetScan/miRBase/miRDeep2;sRNABenchRPM | Method/inputs handoff, not automated |
| Wet-lab validation |Sanger, histology, dual-luciferase assays | Not computationally reproduced |

Population genetic analyses require careful sample subsets. Do not reuse
MAF/HWE/LD-pruned variants for diversity, SFS, demographic inference or WGS
calling without scientific justification. Invariant/callable-site denominators
and ploidy are essential for correct downstream interpretation.
