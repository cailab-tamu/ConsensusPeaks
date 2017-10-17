library(ggplot2)
#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
library(biomaRt)
library(org.Hs.eg.db)



#########################################################################
#       Create Across Cell Type General Plots (Hist&Scatter)
#########################################################################

setwd('//a361-y/TERRA_1T_X/ahmads_workbench/blueprint_peaks/data/Universal_Consensus_Peaks/merged_w_1000/across_cell_types/')
M_ac = read.table("Cell_M_ac_counts.bed", stringsAsFactors = F)
N_ac = read.table("Cell_N_ac_counts.bed", stringsAsFactors = F)
T_ac = read.table("Cell_T_ac_counts.bed", stringsAsFactors = F)

del = M_ac[,4] == 0 | N_ac[,4] == 0 | T_ac[,4] == 0
T_ac_com = T_ac[!del,]
N_ac_com = N_ac[!del,]
M_ac_com = M_ac[!del,]

hist_M_ac_com = hist(M_ac_com[,4],breaks = 15)
hist_N_ac_com = hist(N_ac_com[,4],breaks = 15)
hist_T_ac_com = hist(T_ac_com[,4],breaks = 15)

ac_counts_com = rbind(hist_M_ac_com$counts,hist_N_ac_com$counts, hist_T_ac_com$counts)

ac_freq = ac_counts_com/nrow(M_ac)
rownames(ac_freq) = c("M", "N", "T")
colnames(ac_freq) = paste(round(hist_M_ac_com$mids/145*100), rep("%", 15))
barplot(ac_freq, beside=T,  main="H3K27ac", col=gray.colors(3), ylim=c(0,0.4) ,cex.axis = 0.8, cex.names = 0.8, xlab = "Peak Shared Between Percent Samples", ylab = "% Peaks" )
legend("topleft", legend=rownames(ac_freq), fill=gray.colors(3), cex=0.6)


T_ac_spec = T_ac[M_ac[,4] == 0 & N_ac[,4] == 0 & T_ac[,4] != 0,]
N_ac_spec = N_ac[M_ac[,4] == 0 & N_ac[,4] != 0 & T_ac[,4] == 0,]
M_ac_spec = M_ac[M_ac[,4] != 0 & N_ac[,4] == 0 & T_ac[,4] == 0,]

M_N_ac_spec = M_ac[M_ac[,4] != 0 & N_ac[,4] != 0 & T_ac[,4] == 0,]
M_T_ac_spec = M_ac[M_ac[,4] != 0 & N_ac[,4] == 0 & T_ac[,4] != 0,]
N_T_ac_spec = N_ac[M_ac[,4] == 0 & N_ac[,4] != 0 & T_ac[,4] != 0,]


ac_spec = c(M=nrow(M_ac_spec),N = nrow(N_ac_spec), T = nrow(T_ac_spec), NT = nrow(N_T_ac_spec), MT = nrow(M_T_ac_spec), MN = nrow(M_N_ac_spec) )
barplot(ac_spec/nrow(M_ac),ylim=c(0,0.45))


N_me = read.table("Cell_N_me_counts.bed",stringsAsFactors = F)
T_me = read.table("Cell_T_me_counts.bed", stringsAsFactors = F)
del = N_me[,4] == 0 | T_me[,4] == 0

T_me_com = T_me[!del, ]
N_me_com = N_me[!del, ]

hist_N_me_com = hist(N_me_com[,4],breaks = 10)
hist_T_me_com = hist(T_me_com[,4],breaks = 10)
me_counts_com = rbind(hist_N_me_com$counts, hist_T_me_com$counts)


me_freq = me_counts_com/nrow(N_me)

rownames(me_freq) = c("N", "T")
colnames(me_freq) = paste(round(hist_N_me_com$mids/110*100), rep("%", 11))
barplot(me_freq, beside=T,  main="H3K4me1", col=gray.colors(2), ylim=c(0,0.4) ,cex.axis = 0.8, cex.names = 0.8, xlab = "Peak Shared Between Percent Samples", ylab = "% Peaks" )
legend("topleft", legend=rownames(me_freq), fill=gray.colors(2), cex=0.6)


T_me_spec = T_me[del & T_me[,4] != 0,]
N_me_spec = N_me[del & N_me[,4] != 0,]

me_spec = c(N = nrow(N_me_spec), T = nrow(T_me_spec))
barplot(me_spec/nrow(N_me), ylim=c(0,0.3))




#########################################################################
# Analyze gene data
#########################################################################

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
gene_table = getBM(attributes=c('chromosome_name', 'start_position', 'end_position','hgnc_symbol','ensembl_gene_id'), mart = ensembl)


# Near-Far gene analysis
near_genes_per_peak = list()
for( chr in unique(M_ac[,1])){
  genes = gene_table[paste("chr",gene_table[,1],sep="")==chr,]
  peaks = which(M_ac[,1]==chr)
  
  for( p in peaks){
    diff = abs(genes[,2] - 0.5*(M_ac[p,3]+M_ac[p,2]))
    cg = genes[diff<100000,4]
    near_genes_per_peak[[p]] = list(cg)
  }
}


top = M_ac[,4]==147 & N_ac[,4]==147 & T_ac[,4]==147
near_top_peaks = unlist(near_genes_per_peak[top])
near_top_peaks = near_top_peaks[near_top_peaks != ""]


x = near_top_peaks

ids = bitr(x, fromType = "SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ggo = groupGO(gene = ids[,2], OrgDb = "org.Hs.eg.db", ont="CC", level=3, readable=T)
ggo@result[ggo@result[,3]>0, 1:4]


ego = enrichGO(gene=ids$ENTREZID, OrgDb= org.Hs.eg.db, ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = T)

kk <- enrichKEGG(gene         = ids$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)



# https://www.bioconductor.org/packages/devel/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html






