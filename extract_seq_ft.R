library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(magrittr)

SE_CQN <- readRDS("/Users/zhenwei/Documents/GitHub/m6ALogisticModel/SE_CQN_filtered.rds")

exonsBylongest_tx <- function(txdb , gene_only = F) {

exbtx <- exonsBy(txdb, by = "tx")

tx_length <- sum(width(exbtx))

Map_tx_gene <- suppressMessages( select(txdb,names(exbtx),columns = "GENEID",keytype = "TXID") )

nogene_indx <- is.na(Map_tx_gene$GENEID)

tx_length <- tx_length[!nogene_indx]

Map_tx_gene <- Map_tx_gene[!nogene_indx,]

Map_tx_gene <- Map_tx_gene[order(Map_tx_gene$GENEID),]

longest_indx_lst <- tapply(tx_length, Map_tx_gene$GENEID, function(x) x == max(x))

#test:# identical( gsub("\\..*","", names(unlist( longest_indx_lst ))), Map_tx_gene$GENEID)

Map_tx_gene$longest <- unlist( longest_indx_lst )

if(!gene_only){
  return(exbtx[c(Map_tx_gene$TXID[Map_tx_gene$longest],names(exbtx)[nogene_indx])])
} else {
  return(exbtx[Map_tx_gene$TXID[Map_tx_gene$longest]])
}

}

rrgs <- rowRanges(SE_CQN)

bin_gr = rrgs

txdb = TxDb.Hsapiens.UCSC.hg19.knownGene

Seq_on_mRNA(bin_gr,txdb,bsgnm,2)

library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)

bsgnm = Hsapiens

Seq_on_mRNA <- function(bin_gr,txdb,bsgnm,flank = 50) {

longest_tx <- exonsBylongest_tx(txdb,T)

bin_2_tx <- mapToTranscripts( bin_gr, longest_tx)

#Remove the mapping to multiple tx

bin_2_tx <- bin_2_tx[ !bin_2_tx$xHits  %in% bin_2_tx$xHits[duplicated( bin_2_tx$xHits )]]

bin_2_tx <- bin_2_tx + flank

TX_sequences <- extractTranscriptSeqs(bsgnm,longest_tx)

Seq_char <- substr(TX_sequences[as.character(seqnames(bin_2_tx))] ,start = start( bin_2_tx), stop = end(bin_2_tx))

Char_lst <- split(Seq_char, bin_2_tx$xHits)

Char_lst_reduced <- lapply(Char_lst,paste0,collapse = "")

Map_seq <- rep(NA,length(bin_gr))

names(Map_seq) <- seq_along(bin_gr)

Map_seq[names(Char_lst_reduced)] <- unlist(Char_lst_reduced)

return(Map_seq)
}


