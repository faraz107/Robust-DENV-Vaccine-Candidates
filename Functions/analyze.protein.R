analyze.protein <- function(PROTEIN, NUM_MAX_MISMATCH, GAPratio, ...) {

  # T epitopes of all DENV in DENV1 #######
  
  data_file <- here("Data/MSA_data", paste0(PROTEIN, " DENV1 MSA.fasta"))
  align_MSA <- read.alignment(file = data_file, format = "fasta", forceToLower = FALSE)
  
  cnt <- vector()
  for (i in seq.int(1, align_MSA$nb)) {
    if (sum(s2c(align_MSA$seq[i])=="-")/length(s2c(align_MSA$seq[i])) <= GAPratio) {
      cnt <- c(cnt, i)
    } 
  }
  v <- as.alignment(nb = length(cnt), nam = align_MSA$nam[cnt], 
                    seq = align_MSA$seq[cnt], com = align_MSA$com[cnt])
  align_MSA <- v
  viprMat <- tolower(as.matrix.alignment(x = align_MSA))
  
  viprCons <- consensus(matali = viprMat, method = "majority")
  
  viprMat <- viprMat[, which(viprCons != "-")]
  viprCons <- viprCons[which(viprCons != "-")]
  
  viprMat <- viprMat[, which(viprCons != "x")]
  viprCons <- viprCons[which(viprCons != "x")]
  
  numSeq <- dim(viprMat)[1]
  numPos <- dim(viprMat)[2]
  
  numSeq1 <- numSeq
  numPos1 <- numPos
  
  y <- list()
  for (i in seq.int(1,dim(viprMat)[1])) {
    y[[i]] <- viprMat[i,]
  }
  
  write.fasta(sequences = y, names = rownames(viprMat), file.out = here("Data", paste0("outMSA_", PROTEIN, "_DENV1.fa")))
  
  align_MSA <- read.alignment(file = here("Data", paste0("outMSA_", PROTEIN, "_DENV1.fa")), format = "fasta", forceToLower = FALSE)
  
  x <- BALCONY::calculate_AA_variation(alignment = align_MSA)
  
  # T epitopes of DENV
  {
    data_file <- here("Data", "TepitopesDENV.fa")
    Tepitopes_DENV <- read.fasta(file = data_file, seqtype = "AA", as.string = "TRUE")
    temp <- tolower(as.character(Tepitopes_DENV))
    teiptope_dict <- AAStringSet(temp)
    
    y <- data.frame(start=rep(0,length(teiptope_dict)), end=rep(0,length(teiptope_dict)))
    for (i in 1:length(teiptope_dict)) {
      temp <- (vmatchPattern(pattern = as.character(teiptope_dict[i]), 
                             subject = align_MSA$seq[numSeq], 
                             max.mismatch = NUM_MAX_MISMATCH))
      y$start[i] <- as.double(IRanges::start(temp))
      y$end[i] <- as.double(IRanges::end(temp))
    }
    
    y$length <- y$end - y$start + 1
    
    tep_ind = vector()
    for (i in 1:length(y$length)) {
      if (!is.na(y[i,1])){
        tep_ind = union(tep_ind, seq.int(from = y$start[i], to = y$end[i]))}
    }
    
    y1 <- y
    tep_ind1 <- tep_ind
    
    tep_avg_cons = rep(NA, length(y$start))
    tep_min_cons = rep(NA, length(y$start))
    tep_med_cons = rep(NA, length(y$start))
    
    for (i in 1:length(y$start)) {
      if (!is.na(y[i,1])) {
        tep_avg_cons[i] <- mean(as.double(x$percentage[1, y[i,1]:y[i,2]]), na.rm = TRUE)
        tep_min_cons[i] <- min(as.double(x$percentage[1, y[i,1]:y[i,2]]), na.rm = TRUE)
        tep_med_cons[i] <- median(as.double(x$percentage[1, y[i,1]:y[i,2]]), na.rm = TRUE)
      }
    }
    
    temp1 <- (vcountPDict(pdict = (teiptope_dict), 
                          subject = AAStringSet(align_MSA$seq), 
                          max.mismatch = NUM_MAX_MISMATCH))
    
    temp1 <- rowSums(temp1)
    
    ExEpiTcell_DENV1 <- data.frame(IdD1 = temp1, 
                                   AverageEpConsD1 = tep_avg_cons,
                                   MinEpConsD1 = tep_min_cons,
                                   D1start = y1$start,
                                   D1end = y1$end,
                                   D1length = y1$length)
    
  }
  
  # T epitopes of all DENV in DENV2 #######
  
  data_file <- here("Data/MSA_data", paste0(PROTEIN, " DENV2 MSA.fasta"))
  align_MSA <- read.alignment(file = data_file, format = "fasta", forceToLower = FALSE)
  
  cnt <- vector()
  for (i in seq.int(1, align_MSA$nb)) {
    if (sum(s2c(align_MSA$seq[i])=="-")/length(s2c(align_MSA$seq[i])) <= GAPratio) {
      cnt <- c(cnt, i)
    } 
  }
  v <- as.alignment(nb = length(cnt), nam = align_MSA$nam[cnt], 
                    seq = align_MSA$seq[cnt], com = align_MSA$com[cnt])
  align_MSA <- v
  viprMat <- tolower(as.matrix.alignment(x = align_MSA))
  
  viprCons <- consensus(matali = viprMat, method = "majority")
  
  viprMat <- viprMat[, which(viprCons != "-")]
  viprCons <- viprCons[which(viprCons != "-")]
  
  viprMat <- viprMat[, which(viprCons != "x")]
  viprCons <- viprCons[which(viprCons != "x")]
  
  numSeq <- dim(viprMat)[1]
  numPos <- dim(viprMat)[2]
  
  numSeq2 <- numSeq
  numPos2 <- numPos
  
  y <- list()
  for (i in seq.int(1,dim(viprMat)[1])) {
    y[[i]] <- viprMat[i,]
  }
  
  write.fasta(sequences = y, names = rownames(viprMat), file.out = here("Data", paste0("outMSA_", PROTEIN, "_DENV2.fa")))
  
  align_MSA <- read.alignment(file = here("Data", paste0("outMSA_", PROTEIN, "_DENV2.fa")), format = "fasta", forceToLower = FALSE)
  
  x <- BALCONY::calculate_AA_variation(alignment = align_MSA)
  
  # T epitopes for DENV
  {
    data_file <- here("Data", "TepitopesDENV.fa")
    Tepitopes_DENV <- read.fasta(file = data_file, seqtype = "AA", as.string = "TRUE")
    temp <- tolower(as.character(Tepitopes_DENV))
    teiptope_dict <- AAStringSet(temp)
    
    y <- data.frame(start=rep(0,length(teiptope_dict)), end=rep(0,length(teiptope_dict)))
    for (i in 1:length(teiptope_dict)) {
      temp <- (vmatchPattern(pattern = as.character(teiptope_dict[i]), 
                             subject = align_MSA$seq[numSeq], 
                             max.mismatch = NUM_MAX_MISMATCH))
      y$start[i] <- as.double(IRanges::start(temp))
      y$end[i] <- as.double(IRanges::end(temp))
    }
    
    y$length <- y$end - y$start + 1
    
    tep_ind = vector()
    for (i in 1:length(y$length)) {
      if (!is.na(y[i,1])){
        tep_ind = union(tep_ind, seq.int(from = y$start[i], to = y$end[i]))}
    }
    
    y1 <- y
    tep_ind1 <- tep_ind
    
    tep_avg_cons = rep(NA, length(y$start))
    tep_min_cons = rep(NA, length(y$start))
    tep_med_cons = rep(NA, length(y$start))
    
    for (i in 1:length(y$start)) {
      if (!is.na(y[i,1])) {
        tep_avg_cons[i] <- mean(as.double(x$percentage[1, y[i,1]:y[i,2]]), na.rm = TRUE)
        tep_min_cons[i] <- min(as.double(x$percentage[1, y[i,1]:y[i,2]]), na.rm = TRUE)
        tep_med_cons[i] <- median(as.double(x$percentage[1, y[i,1]:y[i,2]]), na.rm = TRUE)
      }
    }
    
    temp1 <- (vcountPDict(pdict = (teiptope_dict), 
                          subject = AAStringSet(align_MSA$seq), 
                          max.mismatch = NUM_MAX_MISMATCH))
    
    temp1 <- rowSums(temp1)
    
    ExEpiTcell_DENV2 <- data.frame(IdD2 = temp1, 
                                   AverageEpConsD2 = tep_avg_cons,
                                   MinEpConsD2 = tep_min_cons,
                                   D2start = y1$start,
                                   D2end = y1$end,
                                   D2length = y1$length)
    
  }
  
  
  
  # T epitopes of all DENV in DENV3 #######
  
  data_file <- here("Data/MSA_data", paste0(PROTEIN, " DENV3 MSA.fasta"))
  align_MSA <- read.alignment(file = data_file, format = "fasta", forceToLower = FALSE)
  
  cnt <- vector()
  for (i in seq.int(1, align_MSA$nb)) {
    if (sum(s2c(align_MSA$seq[i])=="-")/length(s2c(align_MSA$seq[i])) <= GAPratio) {
      cnt <- c(cnt, i)
    } 
  }
  v <- as.alignment(nb = length(cnt), nam = align_MSA$nam[cnt], 
                    seq = align_MSA$seq[cnt], com = align_MSA$com[cnt])
  align_MSA <- v
  viprMat <- tolower(as.matrix.alignment(x = align_MSA))
  
  viprCons <- consensus(matali = viprMat, method = "majority")
  
  viprMat <- viprMat[, which(viprCons != "-")]
  viprCons <- viprCons[which(viprCons != "-")]
  
  viprMat <- viprMat[, which(viprCons != "x")]
  viprCons <- viprCons[which(viprCons != "x")]
  
  numSeq <- dim(viprMat)[1]
  numPos <- dim(viprMat)[2]
  
  numSeq3 <- numSeq
  numPos3 <- numPos
  
  y <- list()
  for (i in seq.int(1,dim(viprMat)[1])) {
    y[[i]] <- viprMat[i,]
  }
  
  write.fasta(sequences = y, names = rownames(viprMat), file.out = here("Data", paste0("outMSA_", PROTEIN, "_DENV3.fa")))
  
  align_MSA <- read.alignment(file = here("Data", paste0("outMSA_", PROTEIN, "_DENV3.fa")), format = "fasta", forceToLower = FALSE)
  
  x <- BALCONY::calculate_AA_variation(alignment = align_MSA)
  
  # T epitopes for DENV
  {
    data_file <- here("Data", "TepitopesDENV.fa")
    Tepitopes_DENV <- read.fasta(file = data_file, seqtype = "AA", as.string = "TRUE")
    temp <- tolower(as.character(Tepitopes_DENV))
    teiptope_dict <- AAStringSet(temp)
    
    y <- data.frame(start=rep(0,length(teiptope_dict)), end=rep(0,length(teiptope_dict)))
    for (i in 1:length(teiptope_dict)) {
      temp <- (vmatchPattern(pattern = as.character(teiptope_dict[i]), 
                             subject = align_MSA$seq[numSeq], 
                             max.mismatch = NUM_MAX_MISMATCH))
      y$start[i] <- as.double(IRanges::start(temp))
      y$end[i] <- as.double(IRanges::end(temp))
    }
    
    y$length <- y$end - y$start + 1
    
    tep_ind = vector()
    for (i in 1:length(y$length)) {
      if (!is.na(y[i,1])){
        tep_ind = union(tep_ind, seq.int(from = y$start[i], to = y$end[i]))}
    }
    
    y1 <- y
    tep_ind1 <- tep_ind
    
    tep_avg_cons = rep(NA, length(y$start))
    tep_min_cons = rep(NA, length(y$start))
    tep_med_cons = rep(NA, length(y$start))
    
    for (i in 1:length(y$start)) {
      if (!is.na(y[i,1])) {
        tep_avg_cons[i] <- mean(as.double(x$percentage[1, y[i,1]:y[i,2]]), na.rm = TRUE)
        tep_min_cons[i] <- min(as.double(x$percentage[1, y[i,1]:y[i,2]]), na.rm = TRUE)
        tep_med_cons[i] <- median(as.double(x$percentage[1, y[i,1]:y[i,2]]), na.rm = TRUE)
      }
    }
    
    temp1 <- (vcountPDict(pdict = (teiptope_dict), 
                          subject = AAStringSet(align_MSA$seq), 
                          max.mismatch = NUM_MAX_MISMATCH))
    
    temp1 <- rowSums(temp1)
    
    ExEpiTcell_DENV3 <- data.frame(IdD3 = temp1, 
                                   AverageEpConsD3 = tep_avg_cons,
                                   MinEpConsD3 = tep_min_cons,
                                   D3start = y1$start,
                                   D3end = y1$end,
                                   D3length = y1$length)
    
  }
  
  
  
  # T epitopes of all DENV in DENV4 #######
  
  data_file <- here("Data/MSA_data", paste0(PROTEIN, " DENV4 MSA.fasta"))
  align_MSA <- read.alignment(file = data_file, format = "fasta", forceToLower = FALSE)
  
  cnt <- vector()
  for (i in seq.int(1, align_MSA$nb)) {
    if (sum(s2c(align_MSA$seq[i])=="-")/length(s2c(align_MSA$seq[i])) <= GAPratio) {
      cnt <- c(cnt, i)
    } 
  }
  v <- as.alignment(nb = length(cnt), nam = align_MSA$nam[cnt], 
                    seq = align_MSA$seq[cnt], com = align_MSA$com[cnt])
  align_MSA <- v
  viprMat <- tolower(as.matrix.alignment(x = align_MSA))
  
  viprCons <- consensus(matali = viprMat, method = "majority")
  
  viprMat <- viprMat[, which(viprCons != "-")]
  viprCons <- viprCons[which(viprCons != "-")]
  
  viprMat <- viprMat[, which(viprCons != "x")]
  viprCons <- viprCons[which(viprCons != "x")]
  
  numSeq <- dim(viprMat)[1]
  numPos <- dim(viprMat)[2]
  
  numSeq4 <- numSeq
  numPos4 <- numPos
  
  y <- list()
  for (i in seq.int(1,dim(viprMat)[1])) {
    y[[i]] <- viprMat[i,]
  }
  
  write.fasta(sequences = y, names = rownames(viprMat), file.out = here("Data", paste0("outMSA_", PROTEIN, "_DENV4.fa")))
  
  align_MSA <- read.alignment(file = here("Data", paste0("outMSA_", PROTEIN, "_DENV4.fa")), format = "fasta", forceToLower = FALSE)
  
  x <- BALCONY::calculate_AA_variation(alignment = align_MSA)
  
  # T epitopes for DENV
  {
    data_file <- here("Data", "TepitopesDENV.fa")
    Tepitopes_DENV <- read.fasta(file = data_file, seqtype = "AA", as.string = "TRUE")
    temp <- tolower(as.character(Tepitopes_DENV))
    teiptope_dict <- AAStringSet(temp)
    
    y <- data.frame(start=rep(0,length(teiptope_dict)), end=rep(0,length(teiptope_dict)))
    for (i in 1:length(teiptope_dict)) {
      temp <- (vmatchPattern(pattern = as.character(teiptope_dict[i]), 
                             subject = align_MSA$seq[numSeq], 
                             max.mismatch = NUM_MAX_MISMATCH))
      y$start[i] <- as.double(IRanges::start(temp))
      y$end[i] <- as.double(IRanges::end(temp))
    }
    
    y$length <- y$end - y$start + 1
    
    tep_ind = vector()
    for (i in 1:length(y$length)) {
      if (!is.na(y[i,1])){
        tep_ind = union(tep_ind, seq.int(from = y$start[i], to = y$end[i]))}
    }
    
    y1 <- y
    tep_ind1 <- tep_ind
    
    tep_avg_cons = rep(NA, length(y$start))
    tep_min_cons = rep(NA, length(y$start))
    tep_med_cons = rep(NA, length(y$start))
    
    for (i in 1:length(y$start)) {
      if (!is.na(y[i,1])) {
        tep_avg_cons[i] <- mean(as.double(x$percentage[1, y[i,1]:y[i,2]]), na.rm = TRUE)
        tep_min_cons[i] <- min(as.double(x$percentage[1, y[i,1]:y[i,2]]), na.rm = TRUE)
        tep_med_cons[i] <- median(as.double(x$percentage[1, y[i,1]:y[i,2]]), na.rm = TRUE)
      }
    }
    
    temp1 <- (vcountPDict(pdict = (teiptope_dict), 
                          subject = AAStringSet(align_MSA$seq), 
                          max.mismatch = NUM_MAX_MISMATCH))
    
    temp1 <- rowSums(temp1)
    
    ExEpiTcell_DENV4 <- data.frame(IdD4 = temp1, 
                                   AverageEpConsD4 = tep_avg_cons,
                                   MinEpConsD4 = tep_min_cons,
                                   D4start = y1$start,
                                   D4end = y1$end,
                                   D4length = y1$length)
    
  }
  
  
  
  
  # All DENV epitopes combining######
  
  ExEpiTcell_DENVall <- bind_cols(ExEpiTcell_DENV1,
                                  ExEpiTcell_DENV2,
                                  ExEpiTcell_DENV3,
                                  ExEpiTcell_DENV4)
  
  ExEpiTcell_DENVall$`IEDB ID` <- ExEpiTcell_DENVall$`IEDB ID`
  ExEpiTcell_DENVall$`Epitope Sequence` <- ExEpiTcell_DENVall$`Epitope Sequence`
  
  ExEpiTcell_DENVall %>% mutate(AvgScore = 
                                  rowSums(cbind(AverageEpConsD1, AverageEpConsD2, AverageEpConsD3, AverageEpConsD4), na.rm = TRUE)/4) -> ExEpiTcell_DENVall
  
  ExEpiTcell_DENVall %>% mutate(MaxScore = 
                                  pmax(AverageEpConsD1, AverageEpConsD2, AverageEpConsD3, AverageEpConsD4, na.rm = TRUE)) -> ExEpiTcell_DENVall
  
  ExEpiTcell_DENVall %>% mutate(MinScore = 
                                  pmin(AverageEpConsD1, AverageEpConsD2, AverageEpConsD3, AverageEpConsD4, na.rm = TRUE)) -> ExEpiTcell_DENVall
  
  ExEpiTcell_DENVall$IdD1 <- ExEpiTcell_DENVall$IdD1 / numSeq1
  ExEpiTcell_DENVall$IdD2 <- ExEpiTcell_DENVall$IdD2 / numSeq2
  ExEpiTcell_DENVall$IdD3 <- ExEpiTcell_DENVall$IdD3 / numSeq3
  ExEpiTcell_DENVall$IdD4 <- ExEpiTcell_DENVall$IdD4 / numSeq4
  
  numSeq <- c(numSeq1, numSeq2, numSeq3, numSeq4)
  
  list(ExEpiTcell_DENVall, numSeq) 
  
}