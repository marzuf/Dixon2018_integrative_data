suppressPackageStartupMessages(library(IRanges, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))


#############################################################################################
#############################################################################################
# FOR DEBUG
# tolRad=0
# conservThresh=1
# weightValue = 1
# coverageThresh = 0
# tadfileHeader = FALSE
# chrSize=NULL
# ncpu=2
# fileAsInput = TRUE
# all_DT_files <- c("/media/electron/mnt/ed4/marie/scripts/EZH2_final_MAPQ/06_12_50kb_MAPQFILTER/KARPAS/DMSO/BEDfiles/DI/KARPAS_DMSO_chr9_50000_aggregCounts_ICE_domains.bed",
#                   "/media/electron/mnt/ed4/marie/scripts/EZH2_final_MAPQ/06_12_50kb_MAPQFILTER/LY19WT/DMSO/BEDfiles/DI/LY19WT_DMSO_chr9_50000_aggregCounts_ICE_domains.bed")

# get_consensus_TADs(all_DT_files=all_DT_files, tolRad=0, conservThresh=1, weightValue = NULL, coverageThresh = 0,
#                                tadfileHeader = FALSE, chrSize=NULL, ncpu=2, fileAsInput=TRUE, logFile="logtest.txt")


#############################################################################################
#############################################################################################
#' @param all_DT_files LIST containing all data frames (TAD lists) or VECTOR of strings (name of TAD files)
#' @param tolRad Numeric. The # of bp tolerance around start/end positions to consider as matching.
#' @param conservThresh Numeric (0-1). The ratio of input data in which should be conserved to be considered conserved.
#' @param weightValue If null, the weight of each TAD is the ratio of conservation. If a value is provided, each TAD has the same weight.
#' @param coverageThresh Numeric (0-1). If > 0, check that the TADs that have its boundary conserved is covered at least by this ratio. 
#' @param tadfileHeader Boolean. If the TAD files have an header or not (should be the same for all the files).
#' @param chrSize Numeric. Size of the chromosome (if NULL, the maximal end position is considered as end of the chromosome).
#' @param ncpu Numeric. Number of cpu available for running parallel.
#' @param fileAsInput True if all_DT_files contain the path to TAD files (if FALSE, all_DT_files is a list of dataframe objects)


get_consensus_TADs <- function(all_DT_files, tolRad, conservThresh, weightValue = NULL, coverageThresh = 0, 
                               tadfileHeader = FALSE, chrSize=NULL, ncpu=2, fileAsInput=TRUE, logFile=NULL) {
  if(!is.null(logFile)) {
    
    tableVariables <- c("nbrTADs_matching_before_conservThresh", 
                                     "nbrTADs_matching_after_conservThresh",
                                     "nbrTADs_matching_after_coverThresh",
                                     "nbrTADs_beforeWSI",
                                     "nbrTADs_beforeWSI_nodup",
                                     "nbrTADs_afterWSI")
    on.exit( {
	  outDT <- na.omit(outDT)
      # because we want to have all the numbers even if exiting before -> so the boxplots are comparable at the end !
      missingVar <- tableVariables[!tableVariables %in% as.character(outDT$variable)]
      if(length(missingVar) > 0) {
      tmpDT <- data.frame(variable = missingVar, value=rep(0, length(missingVar)), stringsAsFactors = FALSE)
      outDT <- rbind(outDT, tmpDT)
      }
      write.table(outDT, quote=F, sep="\t", row.names = F, col.names=F, file = logFile); cat(paste0("... written: ", logFile, "\n"))
      })
  }
  if(fileAsInput) {
    outDT <- data.frame(
      variable = c(paste0(rep("file", length(all_DT_files)), 1:length(all_DT_files)), "tolRad", "conservThresh", "weightValue", "coverageThresh", "tadfileHeader", "chrSize"),
      value = c(all_DT_files, tolRad, conservThresh, ifelse(is.null(weightValue), "-", weightValue), coverageThresh, tadfileHeader, ifelse(is.null(chrSize), "-", chrSize)),
      stringsAsFactors = FALSE)
    
  } else {
    if(is.null(names(all_DT_files))) {
      names(all_DT_files) <- paste0("tissue", 1:length(all_DT_files))
    }
    if(any(is.na(names(all_DT_files)))){
      names(all_DT_files)[is.na(names(all_DT_files))] <- paste0("tissue", 1:sum(is.na(names(all_DT_files))))
    }
    outDT <- data.frame(
      variable = c(paste0(rep("file", length(all_DT_files)), 1:length(all_DT_files)), "tolRad", "conservThresh", "weightValue", "coverageThresh", "tadfileHeader", "chrSize"),
      value = c(names(all_DT_files), tolRad, conservThresh, ifelse(is.null(weightValue), "-", weightValue), coverageThresh, tadfileHeader, ifelse(is.null(chrSize), "-", chrSize)),
      stringsAsFactors = FALSE)
    
  }
  
  if(!is.null(logFile)) system(paste0("mkdir -p ", dirname(logFile)))
  
  stopifnot(conservThresh >=0 & conservThresh <= 1 )
  stopifnot(coverageThresh >=0 & coverageThresh <= 1 )
  stopifnot(is.numeric(ncpu))
  cat(paste0("... settings - tolerance radius: ", tolRad, "\n"))
  cat(paste0("... settings - conservation threshold: ", conservThresh, "\n"))
  cat(paste0("... settings - coverage threshold: ", coverageThresh, "\n"))
  if(!is.null(weightValue)) {
    cat(paste0("... settings - weight value: ", weightValue, "\n"))
  }else {
    cat(paste0("... settings - weight value: ~ conservation ratio\n"))
  }
  if(!is.null(chrSize)) stopifnot(is.numeric(chrSize))

  if(fileAsInput) stopifnot(all(file.exists(unlist(all_DT_files))))
  cat("... concatenate TAD files \n")
  chromo_allDT <- foreach(i_file=seq_along(all_DT_files), .combine='rbind') %dopar% {
	if(fileAsInput)   {
  	  cat(paste0("... load TADs from: ", all_DT_files[[i_file]], "\n"))
      dt <- read.delim(all_DT_files[[i_file]], header=tadfileHeader, stringsAsFactors=F)
	} else{
	  dt <- all_DT_files[[i_file]]
	}
    stopifnot(ncol(dt) == 3)
    colnames(dt) <- c("chromo", "start", "end")
    stopifnot(length(unique(dt$chromo)) == 1)
    stopifnot(all (dt$end > dt$start) )
    if(fileAsInput){
      dt$tissue <- gsub("(^.+?_.+?)_chr.+", "\\1", basename(all_DT_files[[i_file]]))
    } else {
      dt$tissue <- ifelse(is.null(names(all_DT_files)[[i_file]]), paste0("tissue", i_file) , names(all_DT_files)[[i_file]])
    }
   dt
  }
  stopifnot(nrow(chromo_allDT) > 0)
  stopifnot(length(unique(chromo_allDT$chromo)) == 1)
  stopifnot(!any(is.na(chromo_allDT$start)))
  stopifnot(!any(is.na(chromo_allDT$end)))
  
  if(!is.null(logFile)) {
    tmpDT_0 <- aggregate(start~tissue, FUN=length, data=chromo_allDT)
    tmpDT <- data.frame(
      variable = c(paste0("nbrTADs_", tmpDT_0$tissue), "tot_input_nbrTADs"),
      value = c(tmpDT_0$start, sum(tmpDT_0$start)),
      stringsAsFactors=FALSE
          )    
    outDT <- rbind(outDT, tmpDT)
  }
  
  if(is.null(chrSize)) {
    chrSize <- max(chromo_allDT$end)
  }
  
  stopifnot(is.numeric(chrSize))
  cat(paste0("... found chromo: ", unique(chromo_allDT$chromo), "\n"))
  cat(paste0("... found chromo size: ", chrSize, "\n"))
  # order by start then end before creating IRanges objects
  chromo_allDT <- chromo_allDT[order(chromo_allDT$start, chromo_allDT$end),]
  chromo_allDT$idx <- as.character(1:nrow(chromo_allDT))
  nTissue <- length(unique(chromo_allDT$tissue))
  ################################################# get START overlap
  # create boundary regions for the start
  cat("... get START overlap\n")
  startBR <- IRanges(start = pmax(0,chromo_allDT$start-tolRad/2), end=pmin(chromo_allDT$start+tolRad/2, chrSize))
  elementMetadata(startBR)$tissue <- chromo_allDT$tissue
  elementMetadata(startBR)$idx <- chromo_allDT$idx
  
  startOverlapDT <- as.data.frame(findOverlaps(startBR, startBR))
  colnames(startOverlapDT) <- c("query", "matching")
  startOverlapDT$queryIdx <- elementMetadata(startBR)$idx[startOverlapDT$query]
  startOverlapDT$matchingIdx <- elementMetadata(startBR)$idx[startOverlapDT$matching]
  startOverlapDT$queryTissue <- elementMetadata(startBR)$tissue[startOverlapDT$query]
  startOverlapDT$matchingTissue <- elementMetadata(startBR)$tissue[startOverlapDT$matching]
  
  ################################################# get END overlap
  # create boundary regions for the end
  cat("... get END overlap\n")
  endBR <- IRanges(start = pmax(0,chromo_allDT$end-tolRad/2), end=pmin(chromo_allDT$end+tolRad/2))
  elementMetadata(endBR)$tissue <- chromo_allDT$tissue
  elementMetadata(endBR)$idx <- chromo_allDT$idx
  
  endOverlapDT <- as.data.frame(findOverlaps(endBR, endBR))
  colnames(endOverlapDT) <- c("query", "matching")
  endOverlapDT$queryIdx <- elementMetadata(endBR)$idx[endOverlapDT$query]
  endOverlapDT$matchingIdx <- elementMetadata(endBR)$idx[endOverlapDT$matching]
  endOverlapDT$queryTissue <- elementMetadata(endBR)$tissue[endOverlapDT$query]
  endOverlapDT$matchingTissue <- elementMetadata(endBR)$tissue[endOverlapDT$matching]
  
  stopifnot(all(elementMetadata(endBR)$idx == elementMetadata(startBR)$idx))
  stopifnot(all(elementMetadata(endBR)$idx == chromo_allDT$idx))
  
  # they should all appear once for match with itself
  stopifnot(all(unique(startOverlapDT$queryIdx) == unique(endOverlapDT$queryIdx)))
  
  cat("... filter the candidates to ensure that matching ends come after matching starts\n")
  # for each end position, for each match with a tissue, check in order to retain only the rows
  # -> that have a corresponding start match in the same tissue
  # -> that the position start match comes before the end match
  toKeep <- foreach(i_row = 1:nrow(endOverlapDT), .combine='c') %dopar% {
    # find where the start of the query domain matches in the current tissue
    corresp <- which(startOverlapDT$queryIdx == endOverlapDT$queryIdx[i_row] &
          startOverlapDT$matchingTissue == endOverlapDT$matchingTissue[i_row])
    if(length(corresp) == 0) return(FALSE)
    start_matches <- chromo_allDT$start[chromo_allDT$idx %in% startOverlapDT$matchingIdx[corresp]]
    end_matches <- chromo_allDT$end[chromo_allDT$idx %in% endOverlapDT$matchingIdx[i_row]]
    return(max(end_matches) > min(start_matches))
  }
  
  endOverlapDT <- endOverlapDT[toKeep,]
  
  if(nrow(endOverlapDT) == 0) {
    output_TAD_DT <- data.frame()
    cat(paste0("... number of consensus TADs found: ", nrow(output_TAD_DT), " (total input: ", nrow(chromo_allDT), " from ", nTissue, " datasets)\n"))
    return(output_TAD_DT)
  }

# ADDED 09.10.2018 for check
#cat(paste0("CHECK A: nrow(endOverlapDT) = ", nrow(endOverlapDT), "\n"))
  
  ################################################# get matching
  # ratio of matching tissues needed only if: 1) score required; or 2) matching threshold
  # and then also check the coverage threshold
  # (coverage threshold needed only if conservThreshold)
  if(is.null(weightValue) | conservThresh > 0) {
    cat("... get ratio of matching boundary regions\n")
    # named list: the queryIdx are the entries of the list, and each holds a vector with the corresponding matching tissue
    tissueMatchStart <- setNames(foreach(idx=as.character(unique(startOverlapDT$queryIdx))) %dopar% 
                                   { unique(startOverlapDT$matchingTissue[as.character(startOverlapDT$queryIdx) == idx]) },
                        as.character(unique(startOverlapDT$queryIdx)))
    tissueMatchEnd <- setNames(foreach(idx=as.character(unique(endOverlapDT$queryIdx))) %dopar%
                                 { unique(endOverlapDT$matchingTissue[as.character(endOverlapDT$queryIdx) == idx]) },
                        as.character(unique(endOverlapDT$queryIdx)))
    # take only the queryIdx that possibly match in both Start and End (also, the endOverlapDT has been filtered previously, not the Start)
    # (take only queryIdx that have both start and end match)
    intersectQuery <- intersect(names(tissueMatchStart), names(tissueMatchEnd))
    
    if(length(intersectQuery) == 0){
      output_TAD_DT <- data.frame()
      cat(paste0("... number of consensus TADs found: ", nrow(output_TAD_DT), " (total input: ", nrow(chromo_allDT), " from ", nTissue, " datasets)\n"))
      return(output_TAD_DT)
    }
    
    if(!is.null(logFile)) {
      tmpDT <- data.frame(variable="nbrTADs_matching_before_conservThresh", value = length(intersectQuery), stringsAsFactors = FALSE)
      outDT <- rbind(outDT, tmpDT)
    }
    
    ########## retrieve the TADs that have BOTH start and end matching
    # (should count as match if a given tissue has BOTH start and end that match)
    # for each domain (query), retrieve the tissues that have both start and end matching
    
    #### AFTER TAKING THE INTERCEPT, SHOULD ENSURE THAT THE MATCHING ENDS COME AFTER THE MATCHING START
    # for each of the queryIdx that has match in both start and end, select the tissue that match to both its start and end
    matchingCandidates <- setNames(foreach(x=intersectQuery) %dopar% {intersect(tissueMatchStart[[x]], tissueMatchEnd[[x]])}, intersectQuery)
    # remnove the queryIdx that have no matching candidates
    # remove empty ones (no intersect)
    matchingCandidates <- Filter(function(x) length(x) > 0, matchingCandidates)
    # get the number ratio of tissues matching both start and end 
    # conservTissueMatching: names = idx of the domain; value = ratio of tissue conservation
    conservTissueMatching <- setNames(foreach(i = matchingCandidates, .combine='c') %dopar% {length(i)/nTissue}, names(matchingCandidates))
    # remove the queryIdx that do not pass the conservation threshold
    idxConservedTADs <- names(conservTissueMatching)[conservTissueMatching >= conservThresh ]
    
    if(!is.null(logFile)) {
      tmpDT <- data.frame(variable="nbrTADs_matching_after_conservThresh", value = length(idxConservedTADs), stringsAsFactors = FALSE)
      outDT <- rbind(outDT, tmpDT)
    }
    
    if(length(idxConservedTADs) == 0){
      output_TAD_DT <- data.frame()
      cat(paste0("... number of consensus TADs found: ", nrow(output_TAD_DT), " (total input: ", nrow(chromo_allDT), " from ", nTissue, " datasets)\n"))
      return(output_TAD_DT)
    }
    
    ######### if wanted, check that the domains that match start and end also cover a certain %age of the domain
    ######### do this after filtering a first time for the ratio
    if(coverageThresh > 0) {
      cat("... get ratio of matching domain coverage\n")
      # FOR ALL QUERIES
      # -> retrieve their corresponding domain size
      # -> check each matching candidate to cover enough of the domain
      # => iterate only over the domains that already pass conservation threhsold
      # filteredMatchingCandidates <- foreach(i_query=1:length(matchingCandidates)) %dopar% { # do this after filtering a 1st for ratio of conservation 
      # iterate over the possibly conserved TADs, to check how much the matching domains from other tissues cover it
      filteredMatchingCandidates <- foreach(i_query=idxConservedTADs) %dopar% {
        # retrieve info about the current domain
        domain_end <- chromo_allDT$end[chromo_allDT$idx == i_query]
        domain_start <- chromo_allDT$start[chromo_allDT$idx == i_query]
        domain_size <- domain_end - domain_start + 1
        # no check each candidate
        ratioBpMatch <- unlist(sapply(matchingCandidates[[i_query]], function(x)  {
          # retrieve the matching domain(s) in the current tissue that overlap with the queryIdx TAD
          matchIdx <- which( (chromo_allDT$start >= domain_start-tolRad) & 
                               (chromo_allDT$end <= domain_end+tolRad) & 
                                  chromo_allDT$tissue == x)
          sum(chromo_allDT$end[matchIdx] - chromo_allDT$start[matchIdx] + 1)/domain_size
        }))
        # additional check: I should have one value for each tissue matching, and at this stage, it should be >= than the conservation threshold
        stopifnot(length(ratioBpMatch)/nTissue >= conservThresh )
        trueMatching <- matchingCandidates[[i_query]][ratioBpMatch >= coverageThresh]
        if(length(trueMatching) == 0) {
          return(NA)
        }else{
          return(trueMatching)
        }
      }
      names(filteredMatchingCandidates) <- idxConservedTADs
      # now, filteredMatchingCandidates is similar to matchingCandidates, 
      # but has replaced with NA the tissue that did not match the coverage threshold
      # so check again the ratio of conservation after removing these NA
      conservCoverageTissueMatching <- setNames(foreach(i = filteredMatchingCandidates, .combine='c') %dopar% 
                                                  {length(na.omit(i))/nTissue}, names(filteredMatchingCandidates))
      idxConservedTADs <- names(conservCoverageTissueMatching)[conservCoverageTissueMatching >= conservThresh ]
      finalTissueMatching <- conservCoverageTissueMatching
      
#      if(!is.null(logFile)) {
#        tmpDT <- data.frame(variable="nbrTADs_matching_after_coverThresh", value = length(idxConservedTADs), stringsAsFactors = FALSE)
#        outDT <- rbind(outDT, tmpDT)
#      }
      
    # end-if (coverageThresh > 0) 
    # else: do not check the coverage, no need to update the ratio of tissues that match 
    } else {
      finalTissueMatching <- conservTissueMatching
    }
  # end-if (is.null(weightValue) | conservThresh > 0) 
  # else: no threshold at all -> take all the possible queries
  } else{
    idxConservedTADs <- as.character(unique(startOverlapDT$queryIdx))

  }



  if(!is.null(logFile)) {
    tmpDT <- data.frame(variable="nbrTADs_matching_after_coverThresh", value = length(idxConservedTADs), stringsAsFactors = FALSE)
    outDT <- rbind(outDT, tmpDT)
  }



  if(length(idxConservedTADs) == 0){
    output_TAD_DT <- data.frame()
    cat(paste0("... number of consensus TADs found: ", nrow(output_TAD_DT), " (total input: ", nrow(chromo_allDT), " from ", nTissue, " datasets)\n"))
    return(output_TAD_DT)
  }
  conservedTADs_DT <- chromo_allDT
  conservedTADs_DT <- conservedTADs_DT[conservedTADs_DT$idx %in% idxConservedTADs,]

# ADDED 09.10.2018 for check
#cat(paste0("CHECK B: nrow(conservedTADs_DT) = ", nrow(conservedTADs_DT), "\n"))

  
  # if weightValue is NULL => the weight is the ratio of conservation across tissues
  if(is.null(weightValue)) {
    conservedTADs_DT$score <- finalTissueMatching[as.character(conservedTADs_DT$idx)]
  }else{
    conservedTADs_DT$score <- weightValue  
  }
  
  # order by end position, then start position for consistency over argument order passing
  conservedTADs_DT <- conservedTADs_DT[order(conservedTADs_DT$end, conservedTADs_DT$start),]


# ADDED 09.10.2018 for check
#cat(paste0("CHECK C: nrow(conservedTADs_DT) = ", nrow(conservedTADs_DT), "\n"))

  
  if(!is.null(logFile)) {
    tmpDT <- data.frame(variable="nbrTADs_beforeWSI", value = nrow(conservedTADs_DT), stringsAsFactors = FALSE)
    outDT <- rbind(outDT, tmpDT)
  }
  
  if(nrow(conservedTADs_DT) == 0) {
      output_TAD_DT <- data.frame()
      cat(paste0("... number of consensus TADs found: ", nrow(output_TAD_DT), " (total input: ", nrow(chromo_allDT), " from ", nTissue, " datasets)\n"))
      return(output_TAD_DT)
  } 
  
  # maybe at some point would like to update/do something with the tissue column ?? # not interested for the moment
  # [-> not taking into account tissue column allows to remove duplicates based on chromo/start/end at this point]
  dupTADs <- which(duplicated(conservedTADs_DT[,c("chromo", "start", "end")]))
  # tissue column in final data frame does not make sense after having removed the duplicates
  conservedTADs_DT$tissue <- NULL
  
  cat(paste0("... number of duplicated conserved TADs: ", length(dupTADs), "\n"))
  if(length(dupTADs) > 0)
    conservedTADs_DT <- conservedTADs_DT[-dupTADs,]

  
  if(!is.null(logFile)) {
    tmpDT <- data.frame(variable="nbrTADs_beforeWSI_nodup", value = nrow(conservedTADs_DT), stringsAsFactors = FALSE)
    outDT <- rbind(outDT, tmpDT)
  }

write.table(chromo_allDT[1:3,],file="", sep="\t", quote=F, row.names=F, col.names=T)
  
  cat(paste0("... number of TADs before WSI algorithm (i.e. possibly overlapping): ", nrow(conservedTADs_DT), "\n"))
  
  ##############################################################################################################################
  # find the biggest end that does not overlapp with start
  # "follows" returns the index of the range in subject that a query range in x directly follows.
  cat("... find last non overlapping\n")
  lastNonOverlapping <- follow(IRanges(start = conservedTADs_DT$start, end=conservedTADs_DT$end),
                               subject=IRanges(start = conservedTADs_DT$start, end=conservedTADs_DT$end))
  lastNonOverlapping[is.na(lastNonOverlapping)] <- 0

#cat("lastNonOverlapping = ", lastNonOverlapping, "\n")
  
  # prepare empty vectors for DP
  cat("... start DP to retrieve optimal partition\n")  
  v_vect <- c(0, conservedTADs_DT$score)
  # p_vect should be 1 based as it is used as index of M_vect
  # and the domains should start at the 2nd position
  # -> shift all by +1
  p_vect <- c(1,lastNonOverlapping+1)
  stopifnot(all(p_vect < 1+1:length(p_vect)))
  M_vect <- c(0, rep(NA, nrow(conservedTADs_DT)))
  stopifnot(length(p_vect) == length(M_vect))
  stopifnot(length(M_vect) == length(v_vect))
  
  cat("... start WSI - forward\n")
  for(i in 2:length(M_vect)) {
    if( M_vect[p_vect[i]] + v_vect[i] > M_vect[i-1]) {
      M_vect[i] <- M_vect[p_vect[i]] + v_vect[i]
    } else{
      M_vect[i] <- M_vect[i-1]
    }
  }
  cat("... start WSI - backward\n")
  # backward to get optimal solution
  solutionVect <- c()
  i <- length(M_vect)
  while(i > 1){
    if (v_vect[i]+M_vect[p_vect[i]] >= M_vect[i-1]){
      solutionVect <- c(solutionVect, i)
      i <-  p_vect[i]
    } else{
      i <- i-1
    }
  }

# ADDED 09.10.2018 for check
#save(lastNonOverlapping, file="v2_lastNonOverlapping.Rdata")
#save(conservedTADs_DT, file="v2_conservedTADs_DT.Rdata")
#save(v_vect, file="v2_v_vect.Rdata")
#save(p_vect, file="v2_p_vect.Rdata")
#save(M_vect, file="v2_M_vect.Rdata")
#cat("Rdata files written\n")
#cat(paste0("CHECK D: nrow(conservedTADs_DT) = ", nrow(conservedTADs_DT), "\n"))
#cat(paste0("CHECK D: length(solutionVect) = ", length(solutionVect), "\n"))
#cat(paste0("CHECK D: length(v_vect) = ", length(v_vect), "\n"))
#cat(paste0("CHECK D: length(p_vect) = ", length(p_vect), "\n"))
#cat(paste0("CHECK D: length(M_vect) = ", length(M_vect), "\n"))

  cat("... build consensus set of TADs\n")
  output_TAD_DT <- conservedTADs_DT[rev(solutionVect-1),]
  
  if(!is.null(logFile)) {
    tmpDT <- data.frame(variable="nbrTADs_afterWSI", value = nrow(output_TAD_DT), stringsAsFactors = FALSE)
    outDT <- rbind(outDT, tmpDT)
  }

  if(nrow(output_TAD_DT) == 0) {
    output_TAD_DT <- data.frame()
    cat(paste0("... number of consensus TADs found: ", nrow(output_TAD_DT), " (total input: ", nrow(chromo_allDT), " from ", nTissue, " datasets)\n"))
    return(output_TAD_DT)
  } 

  ### SOME CHECKS BEFORE EXITING  
  # check that the following start is always bigger than previous end
  stopifnot(all(output_TAD_DT$start[-1] > output_TAD_DT$end[-length(output_TAD_DT$end)]))
  # check that all ends are bigger than the starts
  stopifnot(all(output_TAD_DT$end > output_TAD_DT$start))
  ## check the conservation (check that the TADs find conserved indeed fullfill the conservation threshold criterion !)
  foo <- foreach(i = 1:nrow(output_TAD_DT)) %dopar% { 
    start_matchDT <- chromo_allDT[abs(chromo_allDT$start-output_TAD_DT$start[i] )<=tolRad,]
    end_matchDT <- chromo_allDT[abs(chromo_allDT$end-output_TAD_DT$end[i])<=tolRad,]
    tissueMatching <- intersect(start_matchDT$tissue, end_matchDT$tissue)
    stopifnot( (length(tissueMatching) / nTissue) >= conservThresh)
  }
  cat(paste0("... number of consensus TADs found: ", nrow(output_TAD_DT), " (total input: ", nrow(chromo_allDT), " from ", nTissue, " datasets)\n"))
  rownames(output_TAD_DT) <- NULL



  return(output_TAD_DT[,c("chromo", "start", "end")])
}

test_function <- function() {cat("hello\n")}
