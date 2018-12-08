createHTC <- function(file, bin.size, chr, dim = -1, reindex = FALSE, header=TRUE, inputIsFile = TRUE){
  cat("...... start createHTC \n")
  options("scipen"=100, "digits"=6)
  if(inputIsFile) {
    stopifnot(file.exists(file))
    chr.data = read.delim(file,header=header)
    head(chr.data)
    write.table(chr.data[1:10, ], file = "", col.names=T, row.names=F, sep="\t", quote=F)
    colnames(chr.data) = c("binsA","binsB","counts")  
  } else {
    chr.data <- file
    stopifnot(all(colnames(chr.data) == c("binsA", "binsB", "counts")))
  }
  # data from Rao et al. indicate bins by genome coordinates, they need to be turned into indeces
  if(reindex){
    chr.data$binsA = chr.data$binsA/bin.size
    chr.data$binsB = chr.data$binsB/bin.size
  }
  
  chr.data$binsA = chr.data$binsA+1
  chr.data$binsB = chr.data$binsB+1
  
  head(chr.data)
  
  chr.matrix = sparseMatrix(i=chr.data$binsA,j=chr.data$binsB,x=chr.data$counts)
  
  # resize matrix
  if(dim == -1)
    dim = max(ncol(chr.matrix),nrow(chr.matrix))
  
  if(ncol(chr.matrix) < dim){
    nadd = dim-ncol(chr.matrix)
    for(i in 1:nadd)
      chr.matrix = cbind(chr.matrix,rep(0,nrow(chr.matrix)))
  }
  
  if(nrow(chr.matrix) < dim){
    nadd = dim-nrow(chr.matrix)
    for(i in 1:nadd)
      chr.matrix = rbind(chr.matrix,rep(0,ncol(chr.matrix)))
  }
  
  cat(c(dim(chr.matrix),"\t"))
  cat("\n")
  
  stopifnot(ncol(chr.matrix) == dim & nrow(chr.matrix) == dim)
  
  # create symmetric matrix
  # uplo: The default is "U" unless ‘x’ already has a ‘uplo’ slot
  # before taking the upper as reference, just ensure that 3 column format stores the upper triangle matrix
  stopifnot(all(chr.data$binsA <= chr.data$binsB))
  chr.matrix = forceSymmetric(chr.matrix)
  
  #************ TMP
  # tmp_matrix <- as.data.frame(as.matrix(chr.matrix))
  # write.table(tmp_matrix, file = "foo_test_matrix_list.txt", quote=F, sep="\t", row.names=F, col.names=F)
  #*****************
  
  ranges = IRanges(start=seq(1,(nrow(chr.matrix)*bin.size),bin.size),width=rep(bin.size,nrow(chr.matrix)))
  chr.granges = GRanges(seqnames=Rle(c(chr),c(nrow(chr.matrix))), ranges = ranges)
  
  # give a name to each row, here just a consecutive number
  rownames(chr.matrix) = c(1:nrow(chr.matrix))
  colnames(chr.matrix) = c(1:ncol(chr.matrix))
  
  names(chr.granges) = rownames(chr.matrix)
  HTC = HTCexp(chr.matrix,chr.granges,chr.granges)
  return(HTC)
}