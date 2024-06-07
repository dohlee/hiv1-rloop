# hiv1-rloop

![rloop](https://dohlee-bioinfo.sgp1.digitaloceanspaces.com/img/rloop_cropped.jpg)

This repository provides analysis codes and bioinformatics pipelines for article entitled "HIV-1-induced host genomic R-loops dictate HIV-1 integration site selection."

This repository is organized into several directories:

**dripc-seq-processing**

*Related to Fig. 1, 2, 4 and Supplementary Fig. 1, 2, 4*

This directory conatins the processing pipeline for DRIPc-seq processing for R-loop identification. Sequencing reads were quality-controlled using `FastQC` v2.8, and adapters were trimmed using `Trim Galore!` v0.6.6 based on `Cutadapt` v2.8. Trimmed reads were aligned to hg38 reference genome using `bwa` v0.7.17-r1188. Read deduplication and peak calling was done using `MACS` v2.2.7.1. Since it is known that R-loops appear as both narrow and broad peaks in DRIPc-seq read alignment due to its variable length, two independent `MACS2 callpeak` runs were done for narrow and broad peak calling. Narrow and broad peaks were merged using `Bedtools` v2.26.0.


**rna-seq-processing**

*Related to Fig. 2*

This directory contains the RNA-seq data processing pipeline for the quantification of protein-coding gene expression levels. RNA-seq reads were quality-controlled and adapter-trimmed as in DRIPc-seq processing. Processed reads aligned to hg38 reference genome with GENCODE v38 gene annotation using `STAR` v2.7.3a. Gene expression quantification was done using `RSEM` v1.3.1. 

**rna-seq-processing-te**

*Related to Supplementary Fig. 2*

This directory contains the RNA-seq data processing pipeline for the quantification of transposable element expression levels using `TEtranscripts` v2.2.1. Processed reads were first aligned to hg38 reference genome with GENCODE and RepeatMasker TE annotation using `STAR` v2.7.3a. In this case, STAR options were modified as follows in order to utilize multimapping reads in downstream analyses: `--outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 --outMultimapperOrder random --runRNGseed 77 --outSAMmultNmax 1 --outFilterType BySJout --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000`. Expression levels of TEs were quantified as read counts by `TEcount` command.

**is-rloop-colocalization**

*Related to Fig. 4*

This directory contains the analysis pipeline for testing integration site vs R-loop colocalization. Enrichment of integration sites nearby R-loop peaks were tested by randomized permutation test. Randomzied R-loop peaks were generated using `bedtools shuffle` so that the number and the length distribution of the R-loop peaks were fully preserved during the randomization process. Similarly, integration sites were also randomized using `bedtools shuffle` command. Randomization was done for 100 times. Note that ENCODE blacklist regions were excluded while shuffling R-loops and integration sites to exclude inaccessible genomic regions from the analysis. For each of the observed (or randomized) integration sites, the closest observed (or randomized) R-loop peak and the corresponding genomic distance were identified using `bedtools closest` command. The distribution of the genomic distances was displayed to show the local enrichmen tof the integration sites, in terms of the increased proportion of integration sites within the 30kb window centered on R-loops compared to their randomized counterparts.

**is-seq-processing**

*Related to Fig. 4*

This directory contains the processing pipeline for HIV-1 integration site-sequencing processing. Quality-control for HIV-1 integration site-sequencing reads were done using `FastQC` v0.11.9. To discard primers and linkers specific for integration site-sequencing from reads, we used `Cutadapt` v2.8 with options `-u 49 -U 38 --minimum-length 36 --pair-filter any --action trim -q0,0 -a [LINKER] -A TGCTAGAGATTTTCCACACTGACTGGGTCTGAGGG -A GGGTCTGAGGG --no-indels --overlap 12`. This allowed the firts position of a read alignment directly represent the genomic position of HIV-1 integration. Processed reads were aligned to hg38 reference genome using `bwa` v0.7.17-r1188, and integration sites were identified by an in-house Python script. Genomic positions that are supported by more than five read alignments were considered as HIV-1 integration sites.

**vector-is-seq-processing**

*Related to Fig. 5*

This directory contains the pipeline for vector IS-seq data processing. Pipeline details are same with `is-seq-processing`.

## Citation

Park, K., Lee, D., Jeong, J., Lee, S., Kim, S., & Anh, S. Learning the histone codes with large genomic windows and three-dimensional chromatin interactions using transformer. bioRxiv (2024)

```
@article{park2024human,
  title={Human immunodeficiency virus-1 induces and targets host genomic R-loops for viral genome integration},
  author={Park, Kiwon and Lee, Dohoon and Jeong, Jiseok and Lee, Sungwon and Kim, Sun and Ahn, Kwangseog},
  journal={bioRxiv},
  pages={2024--03},
  year={2024},
  publisher={Cold Spring Harbor Laboratory}
}
```
