prepare.datatable <- function(ExEpiTcell, AnalyzedProteins, NameProteins, ID_THRESH, ...) {
  
  for (i in seq.int(1, length(NameProteins))) {
    x <- AnalyzedProteins[[i]]
    x <- bind_cols(ExEpiTcell, x)
    x <- x %>% filter(!(IdD1==0 & IdD2==0 & IdD3==0 & IdD4==0))
    x$Protein <- NameProteins[[i]]
    if (i==1) {
      y <- x
    } 
    if (i>1) {
     y <- bind_rows(y, x)  
    }
  }

  y %>% 
    distinct(`IEDB ID`, .keep_all = TRUE) %>% 
    group_by(`IEDB ID`) %>% 
    mutate(meanID = mean(c(IdD1, IdD2, IdD3, IdD4))) %>% 
    mutate(epiLength = median(na.omit(c(D1length, D2length, D3length, D4length)))) %>% 
    ungroup() %>% 
    mutate(meanIDrank = dense_rank(dplyr::desc(meanID))) %>% 
    mutate(AvgSrank = dense_rank(dplyr::desc(AvgScore))) -> y
  
  y$order4 <- if_else(
    ((y$IdD1>ID_THRESH)&(y$IdD2>ID_THRESH)&(y$IdD3>ID_THRESH)&(y$IdD4>ID_THRESH)), "T", "F")
  
  y$order3 <- if_else(
    ((y$IdD1>ID_THRESH)&(y$IdD2>ID_THRESH)&(y$IdD3>ID_THRESH)) | 
      ((y$IdD1>ID_THRESH)&(y$IdD2>ID_THRESH)&(y$IdD4>ID_THRESH)) |
      ((y$IdD1>ID_THRESH)&(y$IdD3>ID_THRESH)&(y$IdD4>ID_THRESH)) |  
      ((y$IdD2>ID_THRESH)&(y$IdD3>ID_THRESH)&(y$IdD4>ID_THRESH)), "T", "F")
  
  y$order2 <- if_else(
    ((y$IdD1>=ID_THRESH & y$IdD2>=ID_THRESH) | 
       (y$IdD1>=ID_THRESH & y$IdD3>=ID_THRESH) |
       (y$IdD1>=ID_THRESH & y$IdD4>=ID_THRESH) |
       (y$IdD2>=ID_THRESH & y$IdD3>=ID_THRESH) |
       (y$IdD2>=ID_THRESH & y$IdD4>=ID_THRESH) |
       (y$IdD3>=ID_THRESH & y$IdD4>=ID_THRESH)), "T", "F")
  
  y$order1 <- if_else(
    y$IdD1>ID_THRESH | 
      y$IdD2>ID_THRESH |
      y$IdD3>ID_THRESH |  
      y$IdD4>ID_THRESH, "T", "F")
  
 as.data.frame(y)


}