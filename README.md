User guide for package RNAglm
================

Installation
------------

Currently, you could install this package using this command in R.

``` r
devtools::install_github("ZhenWei10/m6ALogisticModel")
```

<hr/>
Motivation
----------

The major usages of this package:

1. Facilitate the data mining process on RNA level high-throughput data, while being able to deal with **confounding transcriptomic features** that are common in RNA genomics.

2. Automatically generate features given bioconductor annotation files, this could reduce the traditional work load for creating and screening the highly significant RNA features one by one, while unable to account the dependencies between those features.

3.  Provide a template for the insightful linear modeling on the relationship between **belonging (dummy), relative position, and length** for a given region of a transcript. For example, the previous researches suggest that exon length, relative position on 3'UTR, and stop codons are important predictors for the presence of m6A mRNA modification under different biological contexts.

We will build the following design for a given transcript region to improve the topological insight gained on that region:
Belong_Region_X + Belong_Region_X::Length_Region_X + Belong_Region_X::Position_Region_X


**Generalized linear model** is a computational technique used here; it can efficiently quantify the scientific and statistical significance for all the features, while adjust their dependencies on each others. Also, the generalized linear model can be applied on the transcriptomic data with different data forms (e.x. real valued, binary, and counts). 

Finally, the generalized linear modeling is coupled with model selection methods to reduce the potential negative effects in presence of large numbers of highly correlated features (high collinearity). Also, model selection can yield robust coefficient estimates that can generate reliable biological interpretations.

The funcions for the transcript feature annotation and the generalized linear modeling are included in this package; while at the same time, users can introduce more genomic features defined by them self using GRanges object. Effective visualization for comparations between multiple models are also implemented in this package.


Usage Demonstration --- Features Annotation
-------------------------------

``` r
  library(m6ALogisticModel)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(fitCons.UCSC.hg19)
  library(phastCons100way.UCSC.hg19)

  Additional_features_hg19 = list(
    HNRNPC_eCLIP = eCLIP_HNRNPC_gr,
    YTHDC1_TREW = YTHDC1_TREW_gr,
    YTHDF1_TREW = YTHDF1_TREW_gr,
    YTHDF2_TREW = YTHDF2_TREW_gr,
    miR_targeted_genes = miR_targeted_genes_grl,
    TargetScan = TargetScan_hg19_gr,
    Verified_miRtargets = verified_targets_gr,
    METTL3_TREW = METTL3_TREW,
    METTL14_TREW = METTL14_TREW,
    WTAP_TREW = WTAP_TREW,
    METTL16_CLIP = METTL16_CLIP,
    ALKBH5_PARCLIP = ALKBH5_PARCLIP,
    FTO_CLIP = FTO_CLIP,
    FTO_eCLIP = FTO_eCLIP
  )

  SE_CQN <- readRDS("SE_CQN_filtered.rds")

  SE_features_added <- predictors_annot(se = SE_example,
                                        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                                        bsgnm = Hsapiens,
                                        fc = fitCons.UCSC.hg19,
                                        pc = phastCons100way.UCSC.hg19,
                                        struct_hybridize = Struc_hg19,
                                        feature_lst = Additional_features_hg19,
                                        hk_genes_list = HK_hg19_eids)
```

