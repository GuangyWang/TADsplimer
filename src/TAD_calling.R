# Title     : TODO
# Objective : TODO
# Created by: guangyu
# Created on: 2020-08-16


edge_row <- function(l){
    m.binseg <- cpt.meanvar(l,test.stat = "Normal", method = "BinSeg")
    change = cpts(m.binseg)
    if(length(change)>2){
        change = change
    }
    return(change)
}

edge_detection <- function(y,fold){
    if(fold <= 1){
    p_lim = 0.6
    }
    else{p_lim = 0.4}
    freq = data.frame(table(unlist(y)))
    freq[,1] <- as.numeric(as.character(freq[,1]))
    freq[,2] <- as.numeric(as.character(freq[,2]))
    freq[,3] <- freq[,2]/max(freq[,2])

    ff = floor(length(freq[,3])/5000)
    point = freq[1,]
    i = 1
    for(i in 1:(ff+1)){
    f0 = freq[(1+(i-1)*5000):min((i*5000),length(freq[,3])),]
    bcp.1 <- bcp(f0[,3],mcmc = 5000)
    pre.p = data.frame(bcp.1$posterior.prob, bcp.1$posterior.mean[,1])
    colnames(pre.p) = c('Probability', 'X1')
    point = as.matrix(f0[pre.p[,1]>0.2,])
    point = point[!is.na(point[,1]),]
    }
    if(!is.null(nrow(point))){
    return(point[-1,])
    }
    else{
    return(point)
    }
}

boundary_probability <- function(points,cutoff){
  points = as.data.frame(points)
  d=dist(points[,1],method = "euclidean")
  if(nrow(points)>0){
    cluster = hclust(d,method = "complete")
    clusterCut <- cutree(cluster, h = cutoff)
    points[,4] = clusterCut
    pk1 = c()
    j=1
    for(i in unique(points[,4])){
      m = max(points[points[,4]==i,2])
      pk1[j]=points[points[,2]==m&points[,4]==i,1]
      j=j+1
    }
    return(pk1)
  }
  else{
    return(points)
  }
}

boundary_detection = function(file_path,cover_lim_up,cover_lim_down,fold) {
  if(fold<=1){
    range = 500
    range2 = 150
    range3 = 20
  }
  else{
    range = 500/fold
    range2 = 150/fold
    range3 = 20/fold
  }

  bps = data.frame(0,0,0,0,0,0)
  contact_map = file(file_path, "r")
  line = readLines(contact_map, n = 1)
  l = as.numeric(as.character(strsplit(line,'\t')[[1]]))
  end = length(l) - 50
  i = 1
  while ( i<end ) {
    line = readLines(contact_map, n = 1)
    l = as.numeric(as.character(strsplit(line,'\t')[[1]]))
    l1=l[i:min((i+range),length(l))]
    l1[l1>cover_lim_up]=cover_lim_up
    l1[l1<=cover_lim_down]=0
    l2=l[i:max((i-range),1)]
    l2[l2>cover_lim_up]=cover_lim_up
    l2[l2<=cover_lim_down]=0
    if ( length(line) == 0 ) {
      break
    }
    else{
      if(i>50){
        edge_right = edge_row(l1)
        edge_left = edge_row(l2)
        if(length(edge_right)>0&length(edge_left)>0){
          bps[i,1] = i
          bps[i,2] = i + edge_right[2]
          bps[i,3] = i - edge_left[2]
          bps[i,4] = i + edge_right[3]
          bps[i,5] = i - edge_left[3]
        }
      }
    }
    i = i+1
  }
  bps[,1] = 1:length(bps[,1])
  boundary_col = bps[(bps[,2]-bps[,1])<range2&(bps[,2]-bps[,1])>range3,2]
  boundary_col = boundary_col[!is.na(boundary_col)]
  boundary_row = bps[(bps[,1]-bps[,3])<range2&(bps[,1]-bps[,3])>range3,3]
  boundary_row = boundary_row[!is.na(boundary_row)]
  close(contact_map)
  out_list = list(bps,boundary_col,boundary_row)
  return(out_list)
}

TAD_detection<-function(boundary_left,boundary_right,band1,distance_cutoff,fold){
  prob_left = edge_detection(boundary_left,fold)
  if(FALSE){
    return(data.frame(0,0,0))
  }
  else{
    if(is.null(nrow(prob_left)) | nrow(prob_left)==0){
      return(data.frame(0,0,0))
    }
    else{
      prob_left2 = prob_left[order(prob_left[,1]),]
      edg_col = boundary_probability(prob_left2,distance_cutoff)
      prob_right = edge_detection(boundary_right,fold)
      if(!is.null(nrow(prob_right))){
        if(nrow(prob_right)>0){
          prob_right2 = prob_right[order(prob_right[,1]),]
          edg_row = boundary_probability(prob_right2,distance_cutoff)
          pnts = data.frame(0,0)
          t=1
          for(i in 1:length(edg_col)){
            for(j in 1:length(edg_row)){
              if(edg_col[i]-edg_row[j]>10&edg_col[i]-edg_row[j]<1000){
                pnts[t,1]=edg_col[i]
                pnts[t,2]=edg_row[j]
                t=t+1
              }
            }
          }

          coordinate = band1
          coordinate = na.omit(coordinate)
          pre_boundary=data.frame()

          for(t in 1:ceiling(length(coordinate[,1])/3000)){
            n1=1+3000*(t-1)
            n2=3000*t
            x = na.omit(coordinate[n1:min(n2,length(coordinate[,1])),])
            xy.kde <- kde2d(x[,1],x[,2], n=(n2-n1)/3,h=20)
            xy.im <- im(t(xy.kde$z), xcol=xy.kde$x, yrow=xy.kde$y)
            pnt_tmp = pnts[pnts[,2]<xy.im$xrange[2]&pnts[,2]>xy.im$xrange[1],]
            pnt_tmp = pnt_tmp[pnt_tmp[,1]<xy.im$yrange[2]&pnt_tmp[,1]>xy.im$yrange[1],]
            c.0 <- sum(xy.kde$z)
            pre_boundary_temp = data.frame()
            if(length(pnt_tmp[,1])>0){
              for(i in 1:length(pnt_tmp[,1])){
                z <- interp.im(xy.im, pnt_tmp[i,2], pnt_tmp[i,1])
                p <- sapply(z, function(a) sum(xy.kde$z[xy.kde$z < a])) / c.0
                pre_boundary_temp[i,1] = p
                pre_boundary_temp[i,2] = pnt_tmp[i,2]
                pre_boundary_temp[i,3] = pnt_tmp[i,1]
              }
              pre_boundary = rbind(pre_boundary,pre_boundary_temp)
            }
          }
          if(fold<=1){
            boundary = pre_boundary[pre_boundary[,1]>0.8,]
          }
          else(boundary = pre_boundary[pre_boundary[,1]>quantile(pre_boundary[,1],0.8),])
          return(boundary)
        }
        else{
          return(data.frame(0,0,0))
        }
      }
      else{
        return(data.frame(0,0,0))
      }
    }
  }
}

TAD_caller <- function(file, color_lim_up, color_lim_down, distance_cutoff = 25, fold = 1){
  boundary_list = boundary_detection(file, color_lim_up,color_lim_down, fold)
  boundary_left = boundary_list[[2]]
  boundary_right = boundary_list[[3]]
  if(length(boundary_right)==0|length(boundary_left)==0){
    return(t(matrix(c(0,0,0))))
  }
  else{
    pair_boundary1 = boundary_list[[1]][,c(3,2)]
    pair_boundary2 = boundary_list[[1]][,c(5,4)]
    boundary = t(matrix(c(1,1,1)))
    for(i in 1:max(1,floor(length(boundary_left)/2000))){
      b1 = boundary_left[(1+(i-1)*2000):min((i*2000),length(boundary_left))]
      b2 = boundary_right[(1+(i-1)*2000):min((i*2000),length(boundary_left))]
      b1 = na.omit(b1)
      b2 = na.omit(b2)
      b3 = pair_boundary1[min(b1,b2):max(b1,b2),]
      b4 = pair_boundary2[min(b1,b2):max(b1,b2),]
      b5 = rbind(as.matrix(b3),as.matrix(b3),as.matrix(b4))
      b5 = b5[order(b5[,1]),]
      boundary_temp = TAD_detection(b1,b2,b5,distance_cutoff,fold)
      boundary = rbind(boundary,as.matrix(boundary_temp))
    }
    return(boundary[-1,])
  }
}

suppressMessages(library("optparse"))
suppressMessages(require("MASS"))
suppressMessages(require("spatstat"))
suppressMessages(require("MASS"))
suppressMessages(require("spatstat"))
suppressMessages(require("changepoint"))
suppressMessages(require("bcp"))
suppressMessages(require("cluster"))

option_list = list(
  make_option(c("-f", "--inputFile"), type="character",
              help="Input file of Hi-C contact map", metavar="character"),
  make_option(c("-u", "--up"), type="numeric",
              help="Up parameter for Hi-C contact map", metavar="character"),
  make_option(c("-d", "--down"), type="numeric",
              help="Down parameter for Hi-C contact map", metavar="character"),
  make_option(c("-o", "--output"), type="character",
              help="Output directory name", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

tad = TAD_caller(opt$inputFile, opt$u, opt$d)
tad = rbind(tad, c(0, 0, 0))
write.table(tad, opt$output, , sep='\t', row.names = F, col.names = F, quote=F)
