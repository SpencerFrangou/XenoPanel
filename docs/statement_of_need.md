# Statement of need

Targeted transcriptomic panels such as DriverMap measure a defined gene set by multiplexed primer-based
amplification. They are inexpensive, sensitive at low input and well suited to small or difficult samples,
which makes them attractive for xenograft and patient-derived xenograft studies in which human tumour cells
grow in a mouse host. Two analysis problems follow that a general RNA-seq pipeline does not solve. The first is
species cross-reactivity: primers designed for human transcripts can amplify mouse template, so part of a
measured value may come from the host rather than the tumour. Read-level tools that assign reads to a species
by alignment (XenofilteR, Disambiguate, BBSplit, Xenome) cannot be applied, because panel data arrive as
amplicon counts per feature rather than as reads with a genomic origin, and the degree of cross-reactivity
is a property of each primer pair that must be measured. The second problem arises once the counts are clean:
when the same cells are profiled in culture, in a tumour and again after explant, the change in a gene
between culture and tissue can be split exactly into a part the cells carry with them and a part imposed by
the tissue, and that split answers a different question from a differential-expression list.

XenoPanel addresses both in one configuration-driven pipeline. The cross-reactivity module takes a dilution
series in which target-species RNA is titrated into a constant mass of host RNA, estimates per feature the
fraction of signal that is host-derived at a chosen operating point, and writes an exclusion mask that is
applied before any comparison. The QC and differential-expression module runs the established edgeR workflow
(filterByExpr, TMM, quasi-likelihood F-tests with TREAT) over user-defined contrasts, with per-sample quality
metrics. The decomposition module applies an additive identity over three states (culture, tissue, explant)
to classify each gene's change as heritable, environmental or mixed, with optional attribution to the tissue
site. The entropy and pathway module scores transcriptome diversity, network entropy (SCENT signalling
entropy rate and CCAT), signature activity (singscore), PROGENy pathway activity and CollecTRI transcription
factor activity, so that a tumour's state can be described as a whole and not only gene by gene.

Inputs, contrasts, thresholds and modules are declared in a YAML file; nothing is specific to one experiment.
The R core writes CSV tables and a session record, and an optional Python layer draws figures from those
tables without recomputing anything. A synthetic example with planted truth and a recovery test are included,
so the pipeline can be installed, run and checked without any research data. XenoPanel was developed for a
study of breast cancer lung relapse in the Zhang laboratory at the University of Toledo and is intended for
any group analysing targeted panel data from mixed-species material.

## References

Ahdesmäki MJ, Gray SR, Johnson JH, Lai Z. Disambiguate: an open-source application for disambiguating two
species in next generation sequencing data from grafted samples. F1000Research 2016;5:2741.
doi:10.12688/f1000research.10082.2

Badia-i-Mompel P, Vélez Santiago J, Braunger J, et al. decoupleR: ensemble of computational methods to infer
biological activities from omics data. Bioinformatics Advances 2022;2:vbac016. doi:10.1093/bioadv/vbac016

Chen Y, Lun ATL, Smyth GK. From reads to genes to pathways: differential expression analysis of RNA-Seq
experiments using Rsubread and the edgeR quasi-likelihood pipeline. F1000Research 2016;5:1438.
doi:10.12688/f1000research.8987.2

Foroutan M, Bhuva DD, Lyu R, Horan K, Cursons J, Davis MJ. Single sample scoring of molecular phenotypes. BMC
Bioinformatics 2018;19:404. doi:10.1186/s12859-018-2435-4

Kluin RJC, Kemper K, Kuilman T, et al. XenofilteR: computational deconvolution of mouse and human reads in
tumor xenograft sequence data. BMC Bioinformatics 2018;19:366. doi:10.1186/s12859-018-2353-5

Schubert M, Klinger B, Klünemann M, et al. Perturbation-response genes reveal signaling footprints in cancer
gene expression. Nature Communications 2018;9:20. doi:10.1038/s41467-017-02391-6

Teschendorff AE, Enver T. Single-cell entropy for accurate estimation of differentiation potency from a
cell's transcriptome. Nature Communications 2017;8:15599. doi:10.1038/ncomms15599
