

balancingSK<- function(x, max_iter=50, eps=1e-4){

cat("balancingSK - debug_normICE.R\n")

    m <- dim(x)[1]

    ## Initialization    
    sum_ss <- matrix(rep(0, m), ncol=1)
    bias <- matrix(rep(1, m), ncol=1)
    old_dbias <- NULL
    ## Remove Diagonal ?
    
    for (it in 1:max_iter){
        message("it=",it," ", Sys.time())

if(it == 171) save(x, file="x_beforeStep1_171.Rdata")
if(it == 172) save(x, file="x_beforeStep1.Rdata")

# any is infinit x TRUE
# any is infinit x_171 FALSE

        ## 1- calculate sum of W over all rows ++
        sum_ds <- rowSums(x, na.rm=TRUE)
        ##sum_ds <- sqrt(rowSums(x^2))

if(it == 171) save(sum_ds, file="sum_ds_afterStep1_171.Rdata")
if(it == 172) save(sum_ds, file="sum_ds_afterStep1.Rdata")

# any is infinit sum_ds_afterStep1 TRUE
# any is infinit sum_ds_afterStep1_171 FALSE

        ## 2- Calculate a vector of corrected ss reads
        ## NOT DONE
        
        ## 3- Calculate vector of bias

if(it == 171) save(sum_ss, file="sum_ss_beforeStep3_171.Rdata")
if(it == 172) save(sum_ss, file="sum_ss_beforeStep3.Rdata")

        dbias <- as.matrix(sum_ds, ncol=1) + sum_ss

if(it == 171) save(dbias, file="dbias_afterStep3_171.Rdata")    
if(it == 172) save(dbias, file="dbias_afterStep3.Rdata")    

# after STEP3 -> any is.inifinite TRUE

        ## 4 - Renormalize bias by its mean valude over non-zero bins to avoid numerical instabilities
        cat("STEP4 - before: any(is.na(dbias)):\n")
        cat(as.character(any(is.na(dbias))), "\n")

        cat("STEP4 - before: all(dbias == 0):\n")
        cat(as.character(all(dbias == 0)), "\n")

if(it == 172) save(dbias, file="dbias_beforeStep4.Rdata")
if(it == 172) mean_dbias <- mean(dbias[dbias!=0])
if(it == 172) save(mean_dbias, file="mean_dbias_beforeStep4.Rdata")

if(it == 171) save(dbias, file="dbias_beforeStep4_171.Rdata")
if(it == 171) mean_dbias <- mean(dbias[dbias!=0])
if(it == 171) save(mean_dbias, file="mean_dbias_beforeStep4_171.Rdata")

        dbias <- dbias/mean(dbias[dbias!=0])
if(it == 172) save(dbias, file="dbias_afterStep4.Rdata")
if(it == 171) save(dbias, file="dbias_afterStep4_171.Rdata")

        cat("STEP4 - after: any(is.na(dbias)):\n")
        cat(as.character(any(is.na(dbias))), "\n")

# meam_dbias_beforeStep4 -> vaut Inf
# after STEP4 -> any is.na becomes TRUE

        ## 5- Set zero values of bias to 1 to avoir 0/0 error
        cat("STEP5 - before: any(is.na(dbias)):\n")
        cat(as.character(any(is.na(dbias))), "\n")

        dbias[dbias==0] <- 1

        cat("STEP5 - after: any(is.na(dbias)):\n")
        cat(as.character(any(is.na(dbias))), "\n")
               
        ## 6- Divide W by bias BiBj for all (i,j) ++++

if(it==171) save(x, file="step6_171_x.Rdata")
if(it==171) save(dbias, file="step6_171_dbias.Rdata")

        x <- x/(dbias %*% t(dbias))

        ## 7- Multiple total vector of bias by additional biases
        ##bias <- bias * dbias

        cat("head(old_dbias):\n")
        cat(head(old_dbias), "\n")
        cat("any(is.na(old_dbias)):\n")
        cat(as.character(any(is.na(old_dbias))), "\n")
        cat("head(dbias):\n")
        cat(head(dbias), "\n")
        cat("STEP7: any(is.na(dbias)):\n")
        cat(as.character(any(is.na(dbias))), "\n")
        cat("eps:\n")
        cat(eps, "\n")

if(it == 172) save(dbias, file="dbias_beforeStep7.Rdata")
if(it == 172) save(old_dbias, file="old_dbias_beforeStep7.Rdata")

if(it == 171) save(dbias, file="dbias_beforeStep7_171.Rdata")
if(it == 171) save(old_dbias, file="old_dbias_beforeStep7_171.Rdata")


        if (!is.null(old_dbias) && sum(abs(old_dbias - dbias))<eps){
            message("Break at iteration ", it)
            break
        }
        old_dbias <- dbias 
    }
    if (it == max_iter){
        message("Did not converged. Stop at iteration ",max_iter)
    }else{
        message("Converge in ",it," iteration")
    }
    return(x)
}





normICE <- function(x, max_iter=50, eps=1e-4, sparse.filter=0.02){

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
        ##     gr <- HiTC:::dimnames2gr(xmat, pattern="\\||\\:|\\-", feat.names=c("name","chr","start", "end"))
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

