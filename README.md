User guide for package m6ALogisticModel
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

The major reasons for me to create this package is to:

1.  Dealing with **confounding genomic features** that are common in RNA genomics.

2.  Reduce the traditional work load for creating and screening the highly explanatory RNA features one by one, while unable to account the dependencies between those features.

3.  Provide a template for the linear modeling on the relationship between **belonging (dummy), relative position, and length** for a given feature of transcript region. The previous researches suggest that these 3 different properties of a region are important predictors for the behaviours of m6A in different biological context.

**Logistic regression modeling** seems to be a reasonable computational technique here, since it can efficiently quantify the scientific and statistical significance for all the features, while adjust their dependencies on each others. The logistic regression is coupled with bayesian model selection method on genaralized linear model for reasons of reduced model qualities and inferential power in the presence of large amount of highly correlated features.

The funcions for transcript feature annotation and the logistic modeling are included in this package; while at the same time, users can introduce more genomic features defined by them self using GRanges object. Effective visualization for comparation between multiple selected models is also implemented in this package.

<hr/>

A case study with this package.
-------------------------------

To demonstrate the application of this package, I provide a case analysis on finding the transcriptomic features that can predict the evolutionary conservation of m6A locations. The m6A locations used were reported by [RMBase2](http://rna.sysu.edu.cn/rmbase/), and the genomic conservation scores used are extracted from the bioconductor packages [phastCons100way.UCSC.hg19](http://www.bioconductor.org/packages/release/data/annotation/html/phastCons100way.UCSC.hg19.html) and [fitCons.UCSC.hg19](http://www.bioconductor.org/packages/release/data/annotation/html/fitCons.UCSC.hg19.html).

``` r
library(m6ALogisticModel)
library(fitCons.UCSC.hg19)
library(phastCons100way.UCSC.hg19)
```

The m6A location information is stored within data `SE_example` of the package `m6ALogisticMode`.

The rowRanges contains the conserved m6A sites reported by RMBase2 between hg19 and mm10.

``` r
library(SummarizedExperiment)
RMBase2_hg19_gr <- rowRanges( SE_example )
```

<hr/>

Define targets (response variables)
-----------------------------------

The PhastCons scores and fitCons scores of the m6A underlying sequences are extracted, each with bins of 1bp, 5bp, and 51bp.

``` r
library(dplyr)

Flank_sizes = c(0,2,25)

PhastCons_matrix <- lapply(Flank_sizes, 
      function(flank_size) scores(phastCons100way.UCSC.hg19,RMBase2_hg19_gr + flank_size)$scores
        ) %>% Reduce(cbind,.)

FitCons_matrix <- lapply(Flank_sizes, 
      function(flank_size) scores(fitCons.UCSC.hg19,RMBase2_hg19_gr + flank_size)$scores
        ) %>% Reduce(cbind,.)

ConsScore_matrix <- cbind(PhastCons_matrix,FitCons_matrix)

colnames(ConsScore_matrix) = paste0( rep(c("PastCons","FitCons"),
                                           each = length(Flank_sizes)),"_",
                                     rep((Flank_sizes)*2+1,2) , "bp")
```

Visualize the distribution of the scores, and choose the cut off values to convert the scores into binary variables.

``` r
Plot_df <- reshape2::melt(ConsScore_matrix)
Plot_df$X_intercept = NA
cut_offs <- c(.9,.9,.8,.7,.7,.7)
for (i in 1:(2*length(Flank_sizes))){
Plot_df$X_intercept[Plot_df$Var2 == colnames(ConsScore_matrix)[i]] <- cut_offs[i]
}
library(ggplot2)
ggplot(Plot_df) + geom_histogram(aes(x = value),fill = "grey") + facet_wrap(~Var2, nrow = 2, scales = "free_y") + theme_classic() + geom_vline(aes(xintercept = X_intercept),colour = "blue")
```

<img src="README_files/figure-markdown_github/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

``` r
#Convert the matrix with only entries of 0, 1, and NA using the cut_off defined above
for(i in 1:6) {
ConsScore_matrix[,i] <- as.numeric(ConsScore_matrix[,i] > cut_offs[i])
}

#Calculate proportions of positive instances in Target
for(i in 1:6) {
 cat( paste0(colnames(ConsScore_matrix)[i],": ", round( mean( ConsScore_matrix[,i],na.rm = T) ,3) , "\n"))
}
```

    ## PastCons_1bp: 0.692
    ## PastCons_5bp: 0.569
    ## PastCons_51bp: 0.513
    ## FitCons_1bp: 0.452
    ## FitCons_5bp: 0.448
    ## FitCons_51bp: 0.351

The targets still face the issue of not having balanced classes, but for now, we will just omit it and plug it into our logistic regression functions.

<hr/>

Generate features (predictors)
------------------------------

The hg19 `TxDb` and `BSgenome` are loaded for the purpose of extracting the transcriptomic features.

``` r
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
SE <- SummarizedExperiment(ConsScore_matrix)
rowRanges(SE) = RMBase2_hg19_gr

Feature_List_hg19 = list(
HNRNPC_eCLIP = eCLIP_HNRNPC_gr,
YTHDC1_TREW = YTHDC1_TREW_gr,
YTHDF1_TREW = YTHDF1_TREW_gr,
YTHDF2_TREW = YTHDF2_TREW_gr,
miR_targeted_genes = miR_targeted_genes_grl,
TargetScan = TargetScan_hg19_gr,
Verified_miRtargets = verified_targets_gr
)

SE <- m6ALogisticModel::predictors.annot(se = SE,
                                   txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                                   bsgnm = BSgenome.Hsapiens.UCSC.hg19,
                                   struct_hybridize = Struc_hg19,
                                   feature_lst = Feature_List_hg19,
                                   HK_genes_list = HK_hg19_eids)
```

    ## Feature 1 : 5'utr is generated.
    ## Feature 2 : 3'utr is generated.
    ## Feature 3 : cds is generated.
    ## Feature 4 : stop codons 201bp is generated.
    ## Feature 5 : start codons 201bp is generated.
    ## Feature 6 : m6Am (Transcription strat site 1bp with A) is generated.
    ## Feature 7 : exon is generated.
    ## Feature 8 : Start 50bp of the last exon is generated.
    ## Feature 9 : relative positioning on 5'utr is generated.
    ## Feature 10 : relative positioning on 3'utr is generated.
    ## Feature 11 : relative positioning on cds is generated.
    ## Feature 12 : relative positioning on exon is generated.
    ## Feature 13 : 3'UTR length (z scores)  is generated.
    ## Feature 14 : 5'UTR length (z scores)  is generated.
    ## Feature 15 : CDS length (z scores)  is generated.
    ## Feature 16 : Long exon (length > 400bp) is generated.
    ## Feature 17 : gene length-exons (z-score) is generated.
    ## Feature 18 : gene length full transcript (z-score) is generated.
    ## Feature 19 : Motif --- AAACA is generated.
    ## Feature 20 : Motif --- GAACA is generated.
    ## Feature 21 : Motif --- AGACA is generated.
    ## Feature 22 : Motif --- GGACA is generated.
    ## Feature 23 : Motif --- AAACT is generated.
    ## Feature 24 : Motif --- GAACT is generated.
    ## Feature 25 : Motif --- AGACT is generated.
    ## Feature 26 : Motif --- GGACT is generated.
    ## Feature 27 : Motif --- AAACC is generated.
    ## Feature 28 : Motif --- GAACC is generated.
    ## Feature 29 : Motif --- AGACC is generated.
    ## Feature 30 : Motif --- GGACC is generated.
    ## Feature 31 : RNA structure --- predicted hybridized region is generated.
    ## Feature 32 : RNA structure --- inferred loop structures between hybridized region is generated.
    ## Feature 33 : annotation feature --- HNRNPC_eCLIP is generated.
    ## Feature 34 : annotation feature --- YTHDC1_TREW is generated.
    ## Feature 35 : annotation feature --- YTHDF1_TREW is generated.
    ## Feature 36 : annotation feature --- YTHDF2_TREW is generated.
    ## Feature 37 : annotation feature --- miR_targeted_genes is generated.
    ## Feature 38 : annotation feature --- TargetScan is generated.
    ## Feature 39 : annotation feature --- Verified_miRtargets is generated.
    ## Feature 40 : snc RNA (<= 200bp) is generated.
    ## Feature 41 : lnc RNA (> 200bp) is generated.
    ## Feature 42 : Isoform number z score is generated.
    ## Feature 43 : House keeping genes is generated.
    ## Feature 44 : Gene level GC content z score is generated.
    ## Feature 45 : 101bp GC content z score is generated.

49 Transcriptomic features on hg19 are automatically attached by the command and data defined in `m6ALogisticModel`.

<hr/>

Model selection and inference
-----------------------------

-   Run the following commands, you will find the results saved under `save_dir` after some running times.
-   You could make this process much more efficient with reduced `MCMC_interations` number and reduced instance numbers in your matrix.

``` r
Group_list <- group_list_default[names(group_list_default) != "Evolution"]

m6ALogisticModel::logistic.modeling(SE,
                                    save_dir = "Conservation_scores",
                                    group_list = Group_list,
                                    MCMC_iterations = 100000)
```

<hr/>

``` r
sessionInfo()
```

    ## R version 3.4.2 (2017-09-28)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS Sierra 10.12.6
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] zh_CN.UTF-8/zh_CN.UTF-8/zh_CN.UTF-8/C/zh_CN.UTF-8/zh_CN.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_3.4.2  backports_1.1.2 magrittr_1.5    rprojroot_1.3-2
    ##  [5] tools_3.4.2     htmltools_0.3.6 yaml_2.1.16     Rcpp_0.12.14   
    ##  [9] stringi_1.1.6   rmarkdown_1.8   knitr_1.18      stringr_1.2.0  
    ## [13] digest_0.6.13   evaluate_0.10.1
