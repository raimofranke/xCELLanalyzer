#functions for xCellAnalyze

#read_xcell function: to read the three files for each run: raw data, annotation, compound names
read_xcell <- function(exp_id, my_filepath){
  xcell_raw <- read.delim(paste0(my_filepath, "xcelldata_", exp_id, "/xcelldata_", exp_id, "_raw.txt"), skip = 2, dec = ",")
  fixed_colnames <- colnames(xcell_raw) %>% str_replace('Y..', '') %>% str_replace('\\.', '')
  colnames(xcell_raw) <- fixed_colnames
  xcell.anno <- read.delim(paste0(my_filepath, "xcelldata_", exp_id, "/xcelldata_", exp_id, "_anno.txt"))
  ma_index <- match(xcell.anno$Position, colnames(xcell_raw))
  compound_names <- as.vector(xcell.anno$Compound)
  colnames(xcell_raw) <- replace(fixed_colnames, ma_index, compound_names)
  
  return(xcell_raw)
}

#edit_df function to start from last meassurement before compound addition
edit_df <- function(exp_id, xcell_raw){
  tadd <- read_lines(paste0(my_filepath, "xcelldata_", exp_id, "/xcelldata_", exp_id, "_raw.txt"), n_max = 1)
  match_time_add <- match(tadd, xcell_raw$TimeInterval)
  xcell_raw_edited <- xcell_raw[(match_time_add-1):(match_time_add+799),]
  xcell_raw_edited <- xcell_raw_edited[,-2]
  xcell_raw_edited <- as_tibble(xcell_raw_edited)
  names(xcell_raw_edited)[1] <- "h"
  return(xcell_raw_edited)
}

#Normalization of cell index values
normalize_xcell <- function(xcell){
  x <- as.matrix(xcell[,2:97])
  norm <- as.vector(x[1,])
  xcell_norm <- sweep(x, 2, norm, FUN = "/")
  xcell_norm <- cbind(xcell[,1], xcell_norm)
  xcell_norm <- xcell_norm[-1,]
  xcell_norm$h <- xcell_norm$h - xcell_norm$h[1]
  return(xcell_norm)
}

#Median polish to identify outliers
do_median_polish <- function(exp_id, xcell_norm) {
  xcell_match <- read_tsv(paste0(my_filepath, "xcelldata_", exp_id, "/xcelldata_", exp_id, "_match.txt"))$compounds
  
  RNGvector<-vector(mode = "numeric")
  result.coleff<-vector(mode = "numeric")
  result_res_sum <- vector(mode = "numeric")
  
  pdf(paste0(my_filepath, "xcelldata", exp_id, ".pdf"), paper = "a4")
  for(i in 1:length(xcell_match)){
    
    
    ma<-grep(xcell_match[i], colnames(xcell_norm)) #matching using reg expr #colnames NOT in succesive order 1,2,3
    z<-xcell_norm[,ma]
    z<-as.data.frame(z)
    z<-z[,order(names(z))]
    mp<-medpolish(z)
    mp.residuals <- data.frame(mp$residuals)
    mp.residuals_abs <- data.frame(abs(mp$residuals))
    res_sum <- apply(mp.residuals_abs, 2, sum)
    rng<-diff(apply(mp.residuals, 2, range)) #calculate range for each column of mp.residuals
    RNGvector<-c(RNGvector,rng)
    result.coleff<-c(result.coleff,mp$col)
    result_res_sum <- c(result_res_sum, res_sum)
    
    plot(x=xcell_norm$h, y=z[,1], type="l", col="blue", main=colnames(z),cex.main=0.8, ylim=c(0,3), ylab="NCI", xlab="t [h]")
    lines(x=xcell_norm$h, y=z[,2], type="l", col="green")
    if (length(z[1,]) > 2){
      lines(x=xcell_norm$h, y=z[,3], type="l", col="red")
    }
    if (length(z[1,]) > 3){
      lines(x=xcell_norm$h, y=z[,4], type="l", col="orange")
    }
    
  }
  
  dev.off()
  
  RNG<-as.data.frame(RNGvector)
  result.coleff<-as.data.frame(result.coleff)
  results.xcell<-cbind(RNG, result.coleff)
  results.xcell$res_sum <- result_res_sum
  return(results.xcell)
}


#Calculate median curves
calculate_median_curves <- function(exp_id, xcell_norm){
  xcell_match <- read_tsv(paste0(my_filepath, "xcelldata_", exp_id, "/xcelldata_", exp_id, "_match.txt"))$compounds
  xcell_median <- matrix(ncol = length(xcell_match), nrow = 800)
  colnames(xcell_median) <- xcell_match
  rownames(xcell_median) <- xcell_norm$h
  for(i in 1:length(xcell_match)){
    name<-xcell_match[i]
    ma<-grep(xcell_match[i], colnames(xcell_norm))
    z<-xcell_norm[,ma]
    xcell_median[,i]<-apply(z, 1, median)
  }
  rownames(xcell_median) <- round(as.numeric(rownames(xcell_median)), digits = 4)
  return(xcell_median)
}

#Normalization with DMSO median TCRP
normalize_dmso <- function(exp_id, xcell_median){
  xcell_median_norm <- sweep(xcell_median, 1, xcell_median[,paste0(exp_id, "_DMSO")], FUN = "-")
  ma <- match(paste0(exp_id, "_DMSO"), colnames(xcell_median_norm))
  xcell_median_norm <- xcell_median_norm[,-ma]
  return(xcell_median_norm)
}

#Removal of DMSO from not normalized median TCRPs
remove_dmso <- function(exp_id, xcell_median){
  ma <- match(paste0(exp_id, "_DMSO"), colnames(xcell_median))
  xcell_median_notnorm <- xcell_median[,-ma]
  return(xcell_median_notnorm)
}


#score each replicate and calculate overall score

score1.function <- function(median.splines, dist.measure){
distmat <- dist(median.splines, method = dist.measure)
distmat <- as.matrix(distmat)
score2 <- c()
members <- c()
best.score <- c()
norm.score <- c()
position <- c()

for (i in 1:length(distmat[,1])){
  distmat.ordered <- distmat[order(distmat[,i]),]
  temp <- unlist(strsplit(toString(colnames(distmat)[i]), "_"))
  tempPos <- grep(paste(temp[1],"+", sep=""), rownames(distmat.ordered), perl=TRUE)
  position <- c(position, tempPos)
  members <- c(members, length(tempPos))
  best.score <- c(sum(c(1:length(tempPos))))
  score <- c(sum(tempPos))
  score2 <- c(score2, score)
  norm.score <- c(norm.score, (best.score/score))
}

# results
res <- list(rep = colnames(distmat),
            score = score2, members = members, normscore=norm.score)

return(res)
}

