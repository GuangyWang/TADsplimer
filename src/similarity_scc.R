# Title     : TODO
# Objective : TODO
# Created by: guangyu
# Created on: 2020-08-17


suppressMessages(library(hicrep))
suppressMessages(library("optparse"))

hicrep_similarity<-function(map1, map2, TAD){
  scc = c()
  for(i in 1:nrow(TAD)){
    start = TAD[i, 2]
    end = TAD[i, 3]
    hic1 = data.frame('chr', 1:(end - start+1), 2:(end - start+2), map1[start:end,start:end])
    hic2 = data.frame('chr', 1:(end - start+1), 2:(end - start+2), map2[start:end,start:end])
    processed <- prep(hic1, hic2, 1, 0, 5)
    scc.out = get.scc(processed, 1, 5)
    scc=c(scc,scc.out$scc)
  }
  TAD2 = data.frame(TAD,scc)
  return(TAD2)
}

option_list = list(
  make_option(c("-a", "--inputFile1"), type="character",
              help="Input file of Hi-C contact map", metavar="character"),
  make_option(c("-b", "--inputFile2"), type="character",
              help="Input file of Hi-C contact map", metavar="character"),
  make_option(c("-t", "--tad"), type="character",
              help="Input file of Hi-C contact map", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

map1 = read.table(opt$inputFile1)
map2 = read.table(opt$inputFile2)
TAD = read.table(opt$tad)

s = hicrep_similarity(map1, map2, TAD)
s = rbind(s, c(0,0,0,0))
write.table(s, paste0(opt$tad, '.scc.txt'), sep='\t', row.names = F, col.names = F, quote=F)
