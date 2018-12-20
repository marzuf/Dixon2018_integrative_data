balancingSK<- function(x, max_iter=50, eps=1e-4){
    m <- dim(x)[1]

    ## Initialization    
    sum_ss <- matrix(rep(0, m), ncol=1)
    bias <- matrix(rep(1, m), ncol=1)
    old_dbias <- NULL
    ## Remove Diagonal ?
    
    for (it in 1:max_iter){
        message("it=",it," ", Sys.time())

        ## 1- calculate sum of W over all rows ++
        sum_ds <- rowSums(x, na.rm=TRUE)
        ##sum_ds <- sqrt(rowSums(x^2))

        ## 2- Calculate a vector of corrected ss reads
        ## NOT DONE
        
        ## 3- Calculate vector of bias
        dbias <- as.matrix(sum_ds, ncol=1) + sum_ss

        ## 4 - Renormalize bias by its mean valude over non-zero bins to avoid numerical instabilities
        dbias <- dbias/mean(dbias[dbias!=0])

        ## 5- Set zero values of bias to 1 to avoir 0/0 error
        dbias[dbias==0] <- 1
               
        ## 6- Divide W by bias BiBj for all (i,j) ++++
        x <- x/(dbias %*% t(dbias))

        ## 7- Multiple total vector of bias by additional biases
        ##bias <- bias * dbias


        ## ADDED MZ: for chr9 -> at iteration 172,
        ## error because there infinite in x from step171, sum_ds so then NA in dbias -> sum(abs(old_dbias - dbias)) vaut NA, ne peut pas évaluer la condition


        if(any(is.na(dbias))) {
          stopifnot(it == 172)
          x <- x_save
          break
        }

        if (!is.null(old_dbias) && sum(abs(old_dbias - dbias))<eps){
            message("Break at iteration ", it)
            break
        }
        old_dbias <- dbias 
        x_save <- x 
    }
    if (it == max_iter){
        message("Did not converged. Stop at iteration ",max_iter)
    }else{
        message("Converge in ",it," iteration")
    }
    return(x)
}


###################################
## IterativeCorNormalization
## ICE normlization
## 
##
## x = HTCexp object or HTClist object
## max_iter = maximum number of iteration to converge
## eps = threshold to converge
## spars.filter = Percentage of row and column to discard based on sparsity (default=0.02)
##
##################################

normICE <- function(x, max_iter=50, eps=1e-4, sparse.filter=0.02){

cat("... START normICE - chr9\n")

    if (inherits(x, "HTCexp")){
        stopifnot(isSymmetric(x))
        idata <- intdata(x)
        gr <- y_intervals(x)
    }else if (inherits(x, "HTClist")){
        idata <- HiTC:::getCombinedContacts(x)
        gr <- HiTC:::getCombinedIntervals(x)
    }

    if (!is.na(sparse.filter)){
        message("Start filtering  ...", Sys.time())
        spars <- apply(idata, 1, function(x){length(which(x==0))}) 
        spars.t <- quantile(spars[spars!=dim(idata)[1]], probs=(1-sparse.filter))
        idx <- which(spars>as.numeric(spars.t))
        idata[idx,] <- 0
        idata[,idx] <- 0
        message("Filter out ",length(idx)," rows and columns ...")
    }
    
    message("Start Iterative Correction ...")
    xmat <- balancingSK(idata, max_iter=max_iter, eps=eps)
    
    if (inherits(x, "HTCexp")){
        intdata(x) <- xmat
    }else if (inherits(x, "HTClist")){
        ##     gr <- dimnames2gr(xmat, pattern="\\||\\:|\\-", feat.names=c("name","chr","start", "end"))
        ##     xgi <- gr[[1]]
        ##     ygi <- gr[[2]]
        ##     rownames(xmat) <- id(ygi)
        ##     colnames(xmat) <- id(xgi)
        if (is.null(gr$xgi))
            x <- HiTC:::splitCombinedContacts(xmat, xgi=gr$ygi, ygi=gr$ygi)
        else
            x <- HiTC:::splitCombinedContacts(xmat, xgi=gr$xgi, ygi=gr$ygi)
    }
    x
}
