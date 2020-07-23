set.seed(123)

count_contact <- function(l, m, slope){
  # l: length of TAD
  # m: mean of interaction T1: 3 ~ 4 T2: 1~2
  # slope: log(interaction) VS log(distance) -0.65 ~ -0.85
  # return: the simulation TAD
  TAD =  matrix(0, nrow = l, ncol = l)
  for(i in 1:l){
    for(j in 1:l){
      m1 = log(abs(i - j))*(slope) + m
      m2 = log(abs(i - j))*(0.17) + 0.15               
      r = rnorm(1, mean = m1, sd = m2^2)
      TAD[i, j] = TAD[i, j] + exp(1)^(r)
    }
  }
  return(TAD)
}

count_contact_background <- function(l, m, slope, sd){
  # l: length of TAD
  # m: mean of interaction T1: 3 ~ 4 T2: 1~2
  # slope: log(interaction) VS log(distance) -0.65 ~ -0.85
  # return: the simulation TAD
  TAD =  matrix(0, nrow = l, ncol = l)
  for(i in 1:l){
    for(j in 1:l){
      m1 = log(abs(i - j))*(slope) + m
      m2 = log(abs(i - j))*(0.17) + sd               
      r = rnorm(1, mean = m1, sd = m2^2)
      TAD[i, j] = TAD[i, j] + exp(1)^(r)
    }
  }
  return(TAD)
}

link_nearby <- function(TAD1, TAD2, m, slope, gama){
  # return linked TAD1 abd TAD2
  # m: mean of interaction (1~2)
  # slope: log(interaction) VS log(distance) -0.65 ~ -0.85
  # gama: contact probility of corner (0~1)
  l = nrow(TAD1) + nrow(TAD2)
  TAD_merge = data.frame(count_contact(l, m, slope))*gama
  TAD_merge[1:nrow(TAD1), 1:nrow(TAD1)] = (1 - gama)*TAD1 + TAD_merge[1:nrow(TAD1), 1:nrow(TAD1)]
  TAD_merge[(nrow(TAD1)+1):l, (nrow(TAD1)+1):l] = (1 - gama)*TAD2 + TAD_merge[(nrow(TAD1)+1):l, (nrow(TAD1)+1):l]
  return(TAD_merge)
}

simulation_TAD <-function(input){
  referenceTAD1 =  read.table(paste0(input,'1.txt'),header = F)
  referenceTAD2 =  read.table(paste0(input,'2.txt'),header = F)
  referenceTAD3 =  read.table(paste0(input,'3.txt'),header = F)
  referenceTAD4 =  read.table(paste0(input,'4.txt'),header = F)
  referenceTAD5 =  read.table(paste0(input,'5.txt'),header = F)
  referenceTAD6 =  read.table(paste0(input,'6.txt'),header = F)
  referenceTAD7 =  read.table(paste0(input,'7.txt'),header = F)
  referenceTAD8 =  read.table(paste0(input,'8.txt'),header = F)
  referenceTAD9 =  read.table(paste0(input,'9.txt'),header = F)
  
  reference = list(referenceTAD1,referenceTAD2,referenceTAD3,referenceTAD4,referenceTAD5,referenceTAD6,referenceTAD7,referenceTAD8,referenceTAD9)
  
  
  TAD_order = ceiling(runif(10, min = 0, max=9))   # order of TAD
  merge_order = ceiling(runif(4, min = 0, max=8))   # order of merge TAD
  gama = runif(10, min = 0.4, max=1)
  
  # count total length
  l = 0
  for(i in TAD_order){
    l = l+nrow(reference[[i]])
  }
  
  # add TAD
  TAD = data.frame(count_contact(l, 3 , -0.6))
  start = 1
  end = 0
  boundary = data.frame(start, end)
  for(i in TAD_order){
    start = end + 1
    end = start + nrow(reference[[i]]) -1
    boundary = rbind(boundary,c(start, end))
    TAD[start:end, start:end] = TAD[start:end, start:end] + reference[[i]]/9
  }
  TAD[is.na(TAD)] <- 50
  
  boundary = boundary[-1,]
  boundary2 = data.frame(boundary,gama)
  
  # add merger region
  TAD_m = TAD
  for(i in merge_order){
    x1 = TAD_m[boundary[i,1]:boundary[i,2], boundary[i,1]:boundary[i,2]] # contact 2 and 3 TADs
    x2 = TAD_m[boundary[i+1,1]:boundary[i+1,2], boundary[i+1,1]:boundary[i+1,2]]
    x =  link_nearby(x1,x2,4, -0.5, gama[i])
    TAD_m[boundary[i,1]:boundary[i+1,2], boundary[i,1]:boundary[i+1,2]] = TAD_m[boundary[i,1]:boundary[i+1,2], boundary[i,1]:boundary[i+1,2]]*1 + x*1
  }
  TAD_m[is.na(TAD_m)] <- 100
  return(list(TAD, TAD_m, boundary2, TAD_order, merge_order))
}


suppressMessages(library("optparse"))
options(warn=-1)

option_list = list(
  make_option(c("-i", "--input"), type="character", 
              help="path of reference data", metavar="character"),
  make_option(c("-o", "--output"), type="character", 
              help="path of output", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

input = opt$input
output = opt$output

simulateTAD = simulation_TAD(input)
write.table(simulateTAD[[1]], file.path(output, 'TAD_split.txt'), sep='\t', row.names = F)
write.table(simulateTAD[[2]], file.path(output, 'TAD_merge.txt'), sep='\t', row.names = F)
write.table(simulateTAD[[3]], file.path(output, 'TAD_boundary.txt'), sep='\t', row.names = F)
write.table(simulateTAD[[4]], file.path(output, 'order_of_all_TAD.txt'), sep='\t', row.names = F)
write.table(simulateTAD[[5]], file.path(output, 'order_of_merged_TAD.txt'), sep='\t', row.names = F)



