# qsub -I -q cpu
# hpc: /home/rintu.kutum/office/chs/dream/dream-seq2expr/predict-gene-expression-2022
# conda activate dream-seq2expr-chaos
# local: /Users/rintukutum/opt/anaconda3/envs/r41
# HPC: ssh rintu.kutum@10.1.4.95

rm(list=ls())
# kaos
library('kaos')
library('doMC')
library('dtt')
tr_data <- readr::read_table(
  'data/train_sequences.txt',col_names = FALSE
)
colnames(tr_data) <- c('seq','expr')
toATGC <- function(x){
  strsplit(tolower(x),split = '')[[1]]
}#Converts the nucleotide sequences into 
#individual character array of nucleotides 

kaos.features <- function(x,res=10){
  xx.kaos <- kaos::cgr(data = x, res = res)
  kaos.vec <- kaos::vectorize(xx.kaos)
  xx.dct <- dtt(x=xx.kaos$matrix,type='dct',variant=4)
  xx.svd <- svd(xx.dct)$d
  
  return(list(kbf = kaos.vec, kds = xx.svd))
}

loops <- seq(from=1, to=nrow(tr_data),by=10e3)
library('doMC')
registerDoMC(50)
tr_features <- foreach(i=1:length(loops))%dopar%{
  message('Running ',i)
  if(length(loops) != i){
    idx_start <- loops[i]
    idx_end <- loops[i+1] -1
  }else{
    idx_start <- loops[i]
    idx_end <- nrow(tr_data)
  }
  promoter_seqs <- tr_data$seq[idx_start:idx_end]
  promoter_nucs <- lapply(promoter_seqs,toATGC)
  promoter_features <- lapply(promoter_nucs,kaos.features)
  promoter_features
}
save(tr_features,file = 'data/tr_features.RData')
#---------

# seq.kaos.feature <- data.frame(
#   do.call('rbind',seq.kaos)
# )
# x.pca <- prcomp(t(as.matrix(seq.kaos.feature)))
# summary(x.pca)
# kaos_pca <- data.frame(
#   PC1 = x.pca$rotation[,1],
#   PC2 = x.pca$rotation[,5],
#   y = tr_data$expr[n_range]
# )
# 
# library(ggplot2)
# ggplot(kaos_pca,aes(PC1,PC2)) +
#   geom_point(aes(color=y),alpha=0.5,size=1)