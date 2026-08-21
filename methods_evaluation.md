# Computational methods evaluation

## Scope

This is a targeted audit of `Translational Profiling of Porewater Microbiomes with Ribo-ITP - Methods.txt` against the current workflow code, configuration, README, and the generic `coverage.ipynb` and `functional_analysis.ipynb` notebooks. It is not replacement Methods prose. Sample-specific results, logs, and notebooks for 0306 and 0313 were excluded, as were omissions in the experimental methods.

The current code is internally coherent enough for a supplementary research codebase, but several statements in the Methods describe a different estimator or filtering scope. Those cases should not be fixed by wording alone until the intended analysis is confirmed.

## Highest-priority discrepancies

1. **The 1,500-bp setting is not MEGAHIT's minimum output length.** `min_contig_size: 1500` is nf-core/mag's threshold for binning and bin QC. MEGAHIT itself retains its default minimum output contig length of 200 bp because no `--min-contig-len` option is supplied. The Methods currently attribute 1,500 bp to assembly.

2. **The Methods incorrectly say that all retained bins passed the 50%/10% criteria.** The collector copies all DAS Tool refined bins. With BUSCO as the sole bin-QC method, `gtdbtk_min_completeness: 50.0` and `gtdbtk_max_contamination: 10.0` determine eligibility for GTDB-Tk classification; they do not filter every bin copied to the supplementary output. Unless a separate filtering step is added, call these files refined bins rather than asserting that every one is a qualifying MAG.

3. **The configured binner set has six members, not five.** COMEBin is enabled in addition to CONCOCT, MaxBin2, MetaBAT2, MetaBinner, and SemiBin2.

4. **The UMI description omits a second trimming operation.** The implemented read structure is `{12-nt 5-prime UMI}{4 discarded nt}{retained insert}`. UMI-tools captures the first 12 nt as `umi_1` and removes the following 4 nt as `discard_1`. Thus a 31--56-nt Cutadapt output becomes a 15--40-nt alignment input. Saying only that a 12-nt UMI was extracted does not explain the 16-nt reduction.

5. **The null model does not use probability proportional to length times raw coverage.** It uses `gene_length * (coverage + 1)`. The same `coverage + 1` pseudocount is used when converting both observed and simulated counts to normalized gene signal. This is a result-defining mismatch.

6. **The null statistic and P-value descriptions are inaccurate.** The code calculates log2 observed/null fold change, the sample standard deviation of the 10,000 simulated category signals, separate enrichment/depletion tail probabilities, and a two-sided absolute-deviation probability. All Monte Carlo P-values use an add-one correction and each family receives a separate Benjamini--Hochberg correction. This is not a generic empirical quantile or a "null distribution of means."

## Metagenomic preprocessing, assembly, and bins

### Statements that should be corrected or expanded

- **Workflow release:** identify nf-core/mag 5.4.2 and the Apptainer profile. The release is pinned by the wrapper and is essential because nf-core defaults and tool containers change.

- **fastp:** the listed quality settings now agree with configuration: paired-end adapter detection, qualified-base threshold 15, front/tail sliding-window cutting with mean quality 15, minimum read length 15, and poly-G trimming. The configuration explicitly keeps PhiX-like reads (`keep_phix: true`), so no PhiX depletion is performed. That choice should be stated if read-cleaning steps are enumerated.

- **Assembly and backmapping:** MEGAHIT is the only enabled assembler. Assembly backmapping uses Bowtie2 `--very-sensitive-local`; `binning_map_mode: own` means each sample's own reads supply its binning coverage. The 1,500-bp threshold applies downstream to binning and bin QC, not MEGAHIT output.

- **Binning and refinement:** list all six enabled binners and name DAS Tool consistently. The configuration uses DAS Tool refinement with score threshold 0.5 and passes only refined bins into post-binning analyses. These settings are currently omitted.

- **Bin QC:** BUSCO alone is enabled (`run_busco: true`; CheckM and CheckM2 are both false), with automatic lineage selection. BUSCO's Complete and Duplicated marker percentages should be named directly. A duplicated BUSCO percentage is a redundancy proxy, not a general estimate of contamination. The configured GTDB boundary is at most 10%, whereas the Methods currently say less than 10%.

- **Taxonomic eligibility versus retention:** nf-core/mag supplies bins meeting at least 50% BUSCO completeness and at most 10% BUSCO duplication/its normalized contamination field to GTDB-Tk. The local collection step does not delete bins failing those limits. Keep these two concepts separate.

- **GTDB-Tk:** report GTDB-Tk 2.5.2 with GTDB release 226. The configuration also fixes minimum aligned amino-acid percentage at 10, minimum ANI-screen alignment fraction at 0.65, and pplacer CPUs at 1. The local handoff propagates only nonempty genus calls from GTDB-classified DAS Tool bins to their member contigs; it does not propagate higher ranks or a placeholder genus.

- **GTDB precedence:** the current code does satisfy the stated policy. A non-null `gtdb_db` forces GTDB-Tk on even if another nf-core parameter disagrees. During annotation, an existing GTDB-derived genus is never overwritten by CAT; CAT fills only contigs absent from that map.

### Rigor and interpretation

- Poly-G trimming is appropriate for data affected by Illumina two-colour no-signal G calls, but can remove real G-rich sequence in other data. Report the sequencing/QC rationale and verify pre/post-trimming fastp or MultiQC profiles.

- Keeping PhiX is reasonable only if PhiX was absent, removed upstream, or shown to be negligible. Otherwise it can enter assemblies or consume mappings. Record that no pipeline PhiX filtering was used and support the choice with QC rather than implying that it happened automatically.

- BUSCO-only assessment is reproducible and consistent with the requested analysis, but it is not a complete modern MAG-quality assessment. Report the selected BUSCO lineage and distinguish marker duplication from contamination. If publication claims depend on formal MAG quality, consider an independent contamination/chimerism check and MIMAG-style reporting as a separate, explicitly approved reanalysis; do not silently re-enable CheckM2.

- Multithreaded MEGAHIT is not guaranteed to be bitwise reproducible. nf-core/mag provides a single-CPU reproducibility option, but changing it would require reassembly. At minimum, preserve the pipeline version, executed configuration, container identifiers, and software-versions report.

## Annotation and genus assignment

### Statements that should be corrected or expanded

- **Bakta:** the wrapper runs Bakta in metagenome mode with `--meta`, preserves contig headers, and annotates each complete assembly. The environment pins Bakta 1.12.0. The configured local database is full database 6.0 dated 2025-02-24. These software/database identifiers are missing.

- **Structural-RNA definition:** only Bakta `rRNA`, `tRNA`, and `tmRNA` features are written to the exclusion BED. Generic Bakta `ncRNA` features are not included. The Ribo paragraph's explicit three-feature list is accurate; elsewhere, the broader label "ncRNA" should not be interpreted as every noncoding RNA feature.

- **Metagenomic CDS coverage:** the code uses `bedtools coverage -mean -sorted` on Bakta CDS intervals and the nf-core assembly-backmapping BAM. The resulting value is mean per-base alignment depth over each CDS, not a count of metagenomic reads assigned uniquely to that gene. It is not strand-specific.

- **KofamScan:** KofamScan 1.3.0 runs on Bakta protein sequences in `detail-tsv` mode with an E-value prefilter of `-E 1e-3`. The parser then retains only rows marked with `*`, which are the hits passing Kofam's KO-specific adaptive threshold. Multiple passing KOs per CDS are retained and deduplicated. The Methods mention the adaptive threshold but omit the E-value prefilter and multiple-hit behavior.

- **CAT:** CAT 6.0.1 is invoked in contig mode on the complete assembly with `-r 7 --sensitive`; `--diamond_mode sensitive` is not the executed syntax. CAT therefore computes calls for all contigs, after which the merge uses only calls missing from the GTDB map. The current prose is acceptable as a description of the final fallback assignments, but not as a literal description of which sequences CAT processes. CAT is optional only when both its database and taxonomy paths are null.

- **CAT parsing:** when CAT reports genus alternatives, the parser takes the genus with the highest reported score on that contig's translated line. The final table records genus and contig only; it does not preserve the method, score, or full lineage.

### Rigor and interpretation

- GTDB-first/CAT-fallback is a defensible hierarchy, but the two systems use different reference databases and classification logic. Preserve assignment source and confidence in a provenance table if taxonomic conclusions are important; the current merged two-column file loses that information and cannot expose GTDB/CAT conflicts.

- Bakta, Kofam, CAT, GTDB, and the KO-to-BRITE mapping are database-dependent. Record release dates, identifiers, and preferably checksums. The code identifies CAT's 2024-12-12 package and the local Bakta/GTDB releases, but the Kofam snapshot and study-specific BRITE mapping provenance are not encoded.

## Ribo-ITP preprocessing and alignment

### Statements that should be corrected or expanded

- **Cutadapt:** the implemented options match the stated adapter, minimum overlap 4, adapter-required filtering, length range 31--56 nt, and quality cutoff 28. The stage uses Cutadapt 5.2 and multiple worker threads.

- **UMI structure:** state the exact regex operation: `^(?P<umi_1>.{12})(?P<discard_1>.{4})`. The first 12 nt are appended to the read identifier as the UMI and bases 13--16 are discarded, leaving the rest of the adapter-trimmed read. UMI-tools 1.1.6 is pinned.

- **Matched reference:** each Ribo library is linked by `meta_id` to its matched metagenome assembly. "Aligned to contigs" is underspecified because it does not say that the per-sample matched assembly, rather than a pooled or external reference, is used.

- **Alignment options:** Bowtie2 2.5.4 runs with seed 1000, `--very-sensitive`, `--end-to-end`, `--score-min L,0,0`, `-k 2`, and `--no-unal`. The zero minimum score makes accepted end-to-end alignments exact under Bowtie2's scoring, while `-k 2` exposes whether a second alignment exists.

- **Filtering order:** the actual sequence is Cutadapt, UMI extraction/4-nt deletion, Bowtie2 alignment, removal of multimappers, removal of structural-RNA overlaps, then UMI deduplication. The current sentence order can be read as deduplication preceding the filters.

- **Multimappers:** secondary records are removed and a primary record carrying an `XS:i` alternative-alignment tag is removed. Thus only unique primary alignments are retained. "Multi-mappers were discarded" is correct but should be tied to `-k 2` and this test.

- **Structural-RNA overlap:** BEDTools removes an alignment with any positive overlap with a Bakta rRNA, tRNA, or tmRNA interval. The command does not request strand agreement, so an opposite-strand overlap is also removed.

- **Deduplication:** UMI-tools uses `--method=directional --read-length --edit-distance-threshold=0`. Deduplication is position/strand and read-length aware, but the zero edit-distance threshold collapses only exact UMI matches; despite the method name, near-neighbour UMI error correction is disabled.

### Rigor and interpretation

- Exact full-length alignment is conservative, but it preferentially loses footprints spanning assembly errors, within-population variants, or divergent organisms. This reference bias should be acknowledged, especially when comparing genera or functions with different assembly quality.

- Discarding every multimapper reduces ambiguous assignment but systematically depletes conserved, repeated, paralogous, and multi-copy genes. Report the choice as a sensitivity/specificity tradeoff rather than neutral cleaning.

- Exact-UMI deduplication avoids merging distinct near-neighbour UMIs, but it does not correct UMI sequencing errors and can therefore leave PCR duplicates. Validate whether this is suitable for the UMI complexity and sequencing error rate.

- The analysis assumes that read orientation and the 5-prime endpoint have the intended biological meaning. Computational QC should demonstrate strand specificity, post-extraction length distribution, CDS enrichment, and start/stop or frame periodicity before treating counts as translational signal.

## CDS counting and metagene analysis

### Statements that should be corrected or expanded

- **Read eligibility:** pysam excludes unmapped, secondary, supplementary, duplicate-flagged, and QC-failed alignments. This is not stated.

- **CDS assignment:** the default is largest positive reference overlap among same-strand CDSs. Each read contributes to at most one CDS. Ties remain with the first CDS in the deterministic coordinate/identifier ordering. The Methods omit the strand requirement and tie behavior.

- **Raw count output:** CDS counts are written to `output/coverage/{RIBO_ID}/{RIBO_ID}_cds_counts.tsv`; they do not overwrite the annotation master sheet. These counts include CDSs on contigs with no genus assignment.

- **Metagene selection:** the generic notebook uses 28-nt reads only, 5-prime reference coordinates, a window from -100 through +99 nt, and same-strand alignments. Start- and stop-eligible gene counts are recorded separately. A CDS boundary must be sufficiently far from a contig end to support the full window plus the read-length allowance.

- **Anchors and estimator:** positions are oriented into transcript direction. The start anchor is the CDS start; the stop anchor is the first base of the terminal codon. For each eligible CDS, each boundary profile is divided by the number of valid, same-strand 28-nt alignments whose 5-prime coordinate lies anywhere within that CDS. Genes with no such internal 5-prime end are skipped, and the normalized profiles are summed without a final division by the number of genes. The Methods' phrase "mapped 5-prime ends aligned to that gene" should specify this read-length, strand, coordinate, and zero-denominator behavior.


### Rigor and interpretation

- Largest-overlap assignment is simple and deterministic, but overlapping same-strand CDSs force an arbitrary single-gene choice, and the tie rule is not biological evidence. Quantify or at least acknowledge ambiguous overlaps if common.

- Restricting the displayed metagene to a single 28-nt class is reasonable only if that class is supported by the observed length distribution. If multiple footprint lengths carry signal, read-length-specific offsets or stratified profiles would be more informative.

- Per-gene normalization prevents highly covered genes from dominating solely through total counts, but gives low-count and high-count genes comparable total weight and can make sparse profiles noisy. Report the eligible-gene counts and consider a minimum-count sensitivity analysis if the metagene shape drives a conclusion.

## Taxonomic and functional analysis

### Statements that should be corrected or expanded

- **Genus filtering point:** normalization is calculated on the complete CDS table, including genes without a genus assignment. A classified-only copy is used for taxonomic abundance and read-length plots; descriptive BRITE summaries and the null model use the complete table. This now agrees with the Methods, but any results generated by the earlier global genus filter must be regenerated.

- **Taxonomic quantities:** the genus plots do not simply transfer labels to and count all aligned reads. The metagenome plot sums per-CDS mean DNA coverage by genus; the Ribo plot sums the per-gene metagenome-normalized Ribo signal by genus. Genera below a 0.005 proportion are pooled as `Other`. A separate read-length plot does label individual Ribo alignments by their contig genus.

- **Descriptive Ribo normalization:** for gene `g`, the code calculates Ribo RPKM from assigned counts and CDS length, then divides by mean DNA coverage. Zero or missing coverage is converted to a normalized value of zero. This should be stated because zero handling affects low-abundance genes. The quantity is a DNA-abundance-normalized relative footprint signal, not conventional translation efficiency, which would require an RNA-expression denominator.

- **BRITE aggregation:** every passing KO is expanded through the two-column KO-to-category mapping, after which duplicate gene--category pairs are removed. Thus multiple KOs or duplicate mapping rows cannot make one gene contribute more than once to the same category, while one gene may still contribute once to each of several distinct categories. `Function unknown` is removed from the descriptive plot, complete category strings define the aggregation groups, and the generic plot displays the top 15 functions and top 10 genera.

### Exact implemented null model

For each gene `g`, including genes without a genus assignment, let `L_g` be CDS length, `D_g` mean metagenomic coverage, and `c_g` assigned Ribo reads. The simulation probability is

`p_g = L_g (D_g + 1) / sum_h[L_h (D_h + 1)]`.

Each of 10,000 replicates draws `N = sum_g c_g` reads from `Multinomial(N, p)`. For both observed and simulated read vectors, gene signal is proportional to

`(c_g / L_g) / (D_g + 1)`

and is rescaled so the sum over genes is one million. Category signal is then summed over unique gene--category assignments. Categories with fewer than two assignment-summed observed reads are omitted before testing. The notebook uses the same BRITE-C mapping for the descriptive and null analyses, seed 42, and one thread by default.

The output reports:

- log2 of observed normalized category signal divided by the simulated mean, with a very small numerical epsilon;
- the sample standard deviation across replicate category signals;
- add-one upper-tail, lower-tail, and absolute-deviation two-sided Monte Carlo P-values; and
- separate Benjamini--Hochberg q-values for enrichment, depletion, and two-sided tests.

This is a Monte Carlo multinomial read-allocation model, not a permutation of observed gene labels or counts. It also uses a different pseudocount-bearing normalization from the descriptive RPKM/coverage calculation, so "calculated in the same way as above" is not correct.

### Rigor and interpretation

- The `coverage + 1` term prevents zero-coverage genes from receiving zero allocation probability and stabilizes division, but it is influential when coverage is low. Its scale and biological rationale should be justified, and sensitivity to other pseudocounts should be considered.

- Filtering categories using the observed outcome (`>=2` assignment-summed reads) changes the tested family and can favor broadly annotated categories. Define this as a prespecified exploratory filter and apply FDR only to the explicitly retained family, as the code does.

- Overlapping BRITE categories make category statistics dependent because categories remain nonexclusive, even though each gene contributes at most once to a given category. Report that nonexclusivity when interpreting category totals and FDR-adjusted results.

- The multinomial model conditions on one library's total assigned reads and treats gene allocations as independent draws from fixed probabilities. It does not model biological-replicate variation, overdispersion, uncertainty in DNA coverage, taxonomy/annotation uncertainty, or correlation among genes and categories. Its P-values support within-library exploratory enrichment/depletion, not population-level biological inference. Biological-replicate conclusions need a replicate-aware model.

- With 10,000 replicates and add-one correction, the smallest attainable Monte Carlo P-value is `1/10001`. Report the number of simulations and seed, and avoid interpreting small differences near this resolution as precise.

## Reproducibility details still missing from the Methods

The Methods should record, without reproducing the README, the analysis-defining software and data provenance:

- nf-core/mag 5.4.2, its containerized software-versions report, and the exact parameter file;
- GTDB-Tk 2.5.2 / release 226 and GTDB eligibility settings;
- Bakta 1.12.0 / full database 6.0 (2025-02-24), KofamScan 1.3.0 and Kofam database snapshot, and CAT 6.0.1 / 2024-12-12 database;
- Cutadapt 5.2, UMI-tools 1.1.6, Bowtie2 2.5.4, SAMtools 1.22, BEDTools 2.31.1, and pysam 0.23.3;
- Python, NumPy, pandas, Matplotlib, and any plotting/statistical package versions used for final figures;
- the construction, retrieval date, and checksum of the KO-to-BRITE-C mapping; and
- the exact YAML/CSV/manifest inputs or a sanitized archival copy, plus random seeds.

## Recommended disposition

Straightforward Methods corrections can be made without reanalysis: the nf-core release, true role of the 1,500-bp threshold, COMEBin, BUSCO-only QC, GTDB release and precedence, annotation versions/options, 12+4 UMI operation, matched reference, filtering order, exact deduplication settings, and detailed CDS/metagene/null parameters.

The code now follows the Methods for per-gene metagene normalization, retention of unclassified genes in functional analyses, and unique gene--category BRITE aggregation; analyses produced before those corrections must be regenerated. The remaining result-defining discrepancy is raw coverage in the Methods versus `coverage + 1` in the null probabilities and normalization. Resolve that scientific choice before changing either the code or prose, and describe the exact null statistics and Monte Carlo tests rather than silently treating them as equivalent.
