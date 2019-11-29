generate.figure2 <- function(ProteinDatatableOut, ...) {

df <- ProteinDatatableOut  
  
grobdt_Ep_mapping <- rasterGrob(image_read(here("Figure/Utils", "Ep_mapping.tif")), just = "center")  
  
df <- df %>% 
  # filter(`MHC Allele Classes` == "I") %>% 
  group_by(`MHC Allele Classes`, Protein) %>% summarise(n = n()) %>% 
  spread(`MHC Allele Classes`, n) %>% 
  mutate("Total" = rowSums(.[2:3])) %>%  
  ungroup() %>% 
  as.data.frame() 

df <- t(rbind(df, c("Total", colSums(df[2:4]))))
toname <- df[1,]
dt <- df[2:4,]
dt <- mapply(dt, FUN=as.numeric)
dt <- as.data.frame(matrix(dt, ncol=11, nrow=3))
rownames(dt) <- c("HLA Class I restricted", "HLA Class II restricted", "Total")
colnames(dt) <- toname

dt_all <- formattable(dt, align = "c", list(
  area(row = 1, col = 1:10) ~ custom_color_tile(brewer.pal(8, "Blues")[1], brewer.pal(8, "Set1")[2]),
  area(row = 2, col = 1:10) ~ custom_color_tile(brewer.pal(8, "Greens")[1], brewer.pal(8, "Set1")[3]),
  area(row = c(3)) ~ custom_title(), area(col = 11) ~ custom_title()))

heading <- formatter("span", style = ~style(font.weight = "bold", font.size = "1.25em"))
colnames(dt_all) <- heading(colnames(dt))
rownames(dt_all) <- heading(rownames(dt))

export_formattable(dt_all, here("Figure", "temp_fig.png"), 10, width = "100%", height = "100%")

grobdt_numTeps <- rasterGrob(image_read(here("Figure", "temp_fig.png")), just = "center")


#DENV1 mapped
{
  
  #class I
  {
    df <- ProteinDatatableOut %>% 
      filter(`MHC Allele Classes` == "I") -> df0
    
    
    # tep_ind
    {  df <- df0 %>% filter(Protein == "NS3")
      
      tep_ind_NS3 = vector()
      for (i in 1:length(df$`IEDB ID`)) {
        if (!is.na(df$D1start[i])){
          tep_ind_NS3 = c(tep_ind_NS3, seq.int(
            from = df$D1start[i], 
            to = df$D1end[i]))}
      }
      
      df <- df0 %>% filter(Protein == "NS5")
      
      tep_ind_NS5 = vector()
      for (i in 1:length(df$`IEDB ID`)) {
        if (!is.na(df$D1start[i])){
          tep_ind_NS5 = c(tep_ind_NS5, seq.int(
            from = df$D1start[i], 
            to = df$D1end[i]))}
      }
      
      df <- df0 %>% filter(Protein == "NS4b")
      
      tep_ind_NS4b = vector()
      for (i in 1:length(df$`IEDB ID`)) {
        if (!is.na(df$D1start[i])){
          tep_ind_NS4b = c(tep_ind_NS4b, seq.int(
            from = df$D1start[i], 
            to = df$D1end[i]))}
      }
      
      df <- df0 %>% filter(Protein == "E")
      
      tep_ind_E = vector()
      for (i in 1:length(df$`IEDB ID`)) {
        if (!is.na(df$D1start[i])){
          tep_ind_E = c(tep_ind_E, seq.int(
            from = df$D1start[i], 
            to = df$D1end[i]))}
      }
      
      df <- df0 %>% filter(Protein == "NS4a")
      
      tep_ind_NS4a = vector()
      for (i in 1:length(df$`IEDB ID`)) {
        if (!is.na(df$D1start[i])){
          tep_ind_NS4a = c(tep_ind_NS4a, seq.int(
            from = df$D1start[i], 
            to = df$D1end[i]))}
      }
      
      df <- df0 %>% filter(Protein == "NS2a")
      
      tep_ind_NS2a = vector()
      for (i in 1:length(df$`IEDB ID`)) {
        if (!is.na(df$D1start[i])){
          tep_ind_NS2a = c(tep_ind_NS2a, seq.int(
            from = df$D1start[i], 
            to = df$D1end[i]))}
      }
      
      df <- df0 %>% filter(Protein == "NS2b")
      
      tep_ind_NS2b = vector()
      for (i in 1:length(df$`IEDB ID`)) {
        if (!is.na(df$D1start[i])){
          tep_ind_NS2b = c(tep_ind_NS2b, seq.int(
            from = df$D1start[i], 
            to = df$D1end[i]))}
      }
      
      df <- df0 %>% filter(Protein == "NS1")
      
      tep_ind_NS1 = vector()
      for (i in 1:length(df$`IEDB ID`)) {
        if (!is.na(df$D1start[i])){
          tep_ind_NS1 = c(tep_ind_NS1, seq.int(
            from = df$D1start[i], 
            to = df$D1end[i]))}
      }
      
      df <- df0 %>% filter(Protein == "prM")
      
      tep_ind_prM = vector()
      for (i in 1:length(df$`IEDB ID`)) {
        if (!is.na(df$D1start[i])){
          tep_ind_prM = c(tep_ind_prM, seq.int(
            from = df$D1start[i], 
            to = df$D1end[i]))}
      }
      
      df <- df0 %>% filter(Protein == "C")
      
      tep_ind_C = vector()
      for (i in 1:length(df$`IEDB ID`)) {
        if (!is.na(df$D1start[i])){
          tep_ind_C = c(tep_ind_C, seq.int(
            from = df$D1start[i], 
            to = df$D1end[i]))}
      }
      
    }

    tep_df_I <- data.frame(
      values = c(tep_ind_E, tep_ind_prM, tep_ind_C, tep_ind_NS1, tep_ind_NS2a, tep_ind_NS2b, tep_ind_NS3, tep_ind_NS4a, tep_ind_NS4b, tep_ind_NS5),
      prot = c(rep("E", length(tep_ind_E)), rep("prM", length(tep_ind_prM)), rep("C", length(tep_ind_C)), 
               rep("NS1", length(tep_ind_NS1)), rep("NS2a", length(tep_ind_NS2a)), rep("NS2b", length(tep_ind_NS2b)),
               rep("NS3", length(tep_ind_NS3)), rep("NS4a", length(tep_ind_NS4a)), rep("NS4b", length(tep_ind_NS4b)),
               rep("NS5", length(tep_ind_NS5))))
    tep_df_I$class <- "I"
    
    
  }
  #classII
  {
    df <- ProteinDatatableOut %>% 
      filter(`MHC Allele Classes` == "II") -> df0
    
    # tep_ind
    {  df <- df0 %>% filter(Protein == "NS3")
      
      tep_ind_NS3 = vector()
      for (i in 1:length(df$`IEDB ID`)) {
        if (!is.na(df$D1start[i])){
          tep_ind_NS3 = c(tep_ind_NS3, seq.int(
            from = df$D1start[i], 
            to = df$D1end[i]))}
      }
      
      df <- df0 %>% filter(Protein == "NS5")
      
      tep_ind_NS5 = vector()
      for (i in 1:length(df$`IEDB ID`)) {
        if (!is.na(df$D1start[i])){
          tep_ind_NS5 = c(tep_ind_NS5, seq.int(
            from = df$D1start[i], 
            to = df$D1end[i]))}
      }
      
      df <- df0 %>% filter(Protein == "NS4b")
      
      tep_ind_NS4b = vector()
      for (i in 1:length(df$`IEDB ID`)) {
        if (!is.na(df$D1start[i])){
          tep_ind_NS4b = c(tep_ind_NS4b, seq.int(
            from = df$D1start[i], 
            to = df$D1end[i]))}
      }
      
      df <- df0 %>% filter(Protein == "E")
      
      tep_ind_E = vector()
      for (i in 1:length(df$`IEDB ID`)) {
        if (!is.na(df$D1start[i])){
          tep_ind_E = c(tep_ind_E, seq.int(
            from = df$D1start[i], 
            to = df$D1end[i]))}
      }
      
      df <- df0 %>% filter(Protein == "NS4a")
      
      tep_ind_NS4a = vector()
      for (i in 1:length(df$`IEDB ID`)) {
        if (!is.na(df$D1start[i])){
          tep_ind_NS4a = c(tep_ind_NS4a, seq.int(
            from = df$D1start[i], 
            to = df$D1end[i]))}
      }
      
      df <- df0 %>% filter(Protein == "NS2a")
      
      tep_ind_NS2a = vector()
      for (i in 1:length(df$`IEDB ID`)) {
        if (!is.na(df$D1start[i])){
          tep_ind_NS2a = c(tep_ind_NS2a, seq.int(
            from = df$D1start[i], 
            to = df$D1end[i]))}
      }
      
      df <- df0 %>% filter(Protein == "NS2b")
      
      tep_ind_NS2b = vector()
      for (i in 1:length(df$`IEDB ID`)) {
        if (!is.na(df$D1start[i])){
          tep_ind_NS2b = c(tep_ind_NS2b, seq.int(
            from = df$D1start[i], 
            to = df$D1end[i]))}
      }
      
      df <- df0 %>% filter(Protein == "NS1")
      
      tep_ind_NS1 = vector()
      for (i in 1:length(df$`IEDB ID`)) {
        if (!is.na(df$D1start[i])){
          tep_ind_NS1 = c(tep_ind_NS1, seq.int(
            from = df$D1start[i], 
            to = df$D1end[i]))}
      }
      
      df <- df0 %>% filter(Protein == "prM")
      
      tep_ind_prM = vector()
      for (i in 1:length(df$`IEDB ID`)) {
        if (!is.na(df$D1start[i])){
          tep_ind_prM = c(tep_ind_prM, seq.int(
            from = df$D1start[i], 
            to = df$D1end[i]))}
      }
      
      df <- df0 %>% filter(Protein == "C")
      
      tep_ind_C = vector()
      for (i in 1:length(df$`IEDB ID`)) {
        if (!is.na(df$D1start[i])){
          tep_ind_C = c(tep_ind_C, seq.int(
            from = df$D1start[i], 
            to = df$D1end[i]))}
      }
      
    }

    tep_df_II <- data.frame(
      values = c(tep_ind_E, tep_ind_prM, tep_ind_C, tep_ind_NS1, tep_ind_NS2a, tep_ind_NS2b, tep_ind_NS3, tep_ind_NS4a, tep_ind_NS4b, tep_ind_NS5),
      prot = c(rep("E", length(tep_ind_E)), rep("prM", length(tep_ind_prM)), rep("C", length(tep_ind_C)), 
               rep("NS1", length(tep_ind_NS1)), rep("NS2a", length(tep_ind_NS2a)), rep("NS2b", length(tep_ind_NS2b)),
               rep("NS3", length(tep_ind_NS3)), rep("NS4a", length(tep_ind_NS4a)), rep("NS4b", length(tep_ind_NS4b)),
               rep("NS5", length(tep_ind_NS5))))
    tep_df_II$class <- "II"
    
  }
  
  # plotting 
  {   
    tep_df <- rbind(tep_df_I, tep_df_II)

    tep_df$prot <- factor(tep_df$prot, levels = c("E", "prM", "C", "NS1", "NS2a", "NS2b", "NS3", "NS4a", "NS4b", "NS5"), 
                          labels = c("E", "prM", "C", "NS1", "NS2a", "NS2b", "NS3", "NS4a", "NS4b", "NS5"))
    
    { g <- ggplot(data = tep_df %>% filter(prot %in% c("NS5")), mapping = aes(x = values, color = class))
      # g <- g + geom_histogram(binwidth = 1)
      # g <- g + geom_freqpoly(binwidth = 1, size = 0.5, position = "identity")
      # g <- g + geom_area(stat = "bin", binwidth = 1, alpha = 0.25, position = "identity", aes(fill = class, y = ..count..))
      g <- g + geom_step(aes(y = ..count..), stat = "bin", binwidth = 1, position = "identity", size = 0.5, show.legend = FALSE)
      g <- g + scale_color_manual(values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]))
      
      g <- g + scale_fill_manual(labels = c("I", "II"), values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]), name = "MHC Allele Class")
      
      g <- g + scale_x_continuous(breaks = seq.int(50, 900, 50), limits = c(1, 900), expand = c(0,5))
      g <- g + scale_y_continuous(breaks = seq.int(0, 12, 4), limits = c(0,12))
      g <- g + ylab("") + xlab("Position") 
      # + ggtitle("T cell epitopes coverage of DENV proteins")
      # g <- g + facet_grid(rows = vars(prot), margins = FALSE, shrink = TRUE)
      g <- g + facet_wrap(facets = vars(prot), shrink = TRUE, scales = "free_x", nrow = 1)
      g <- g + theme_minimal() + theme(strip.background = element_rect(color = "white", fill = "grey95"), 
                                       legend.position = "bottom", 
                                       axis.title.y = element_text(margin = margin(t=50, r = 1, b = 1, l = 1), hjust = 0.5),
                                       axis.ticks = element_line(), panel.grid.minor = element_blank(),
                                       strip.text = element_text(face = "bold", size = 10), 
                                       plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
      gh_NS5 <- g}
    
    {    g <- ggplot(data = tep_df %>% filter(prot %in% c("NS4a")), mapping = aes(x = values, y = ..count.., color = class))
      # g <- g + geom_histogram(binwidth = 1, color = brewer.pal(9, "Set1")[2])
      # g <- g + geom_freqpoly(binwidth = 1, size = 0.5, position = "identity")
      # g <- g + scale_color_manual(labels = c("I", "II"), values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]), name = "MHC Allele Class")
      g <- g + geom_step(aes(y = ..count..), stat = "bin", binwidth = 1, position = "identity", size = 0.5, show.legend = FALSE)
      g <- g + scale_color_manual(values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]))
      
      g <- g + scale_fill_manual(labels = c("I", "II"), values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]), name = "MHC Allele Class")
      g <- g + scale_x_continuous(breaks = seq.int(50, 900, 50), limits = c(1, 127), expand = c(0,5))
      g <- g + scale_y_continuous(breaks = seq.int(0, 12, 4), limits = c(0,12))
      g <- g + ylab("") + xlab("") 
      # + ggtitle("T cell epitopes coverage of DENV proteins")
      # g <- g + facet_grid(rows = vars(prot), margins = FALSE, shrink = TRUE)
      g <- g + facet_wrap(facets = vars(prot), shrink = TRUE, scales = "free_x", nrow = 1)
      g <- g + theme_minimal() + theme(strip.background = element_rect(color = "white", fill = "grey95"), 
                                       legend.position = "bottom", 
                                       axis.title.y = element_text(margin = margin(t=50, r = 1, b = 1, l = 1), hjust = 0.5),
                                       axis.ticks = element_line(), panel.grid.minor = element_blank(),
                                       strip.text = element_text(face = "bold", size = 10), 
                                       plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                                       axis.text.y = element_blank(), axis.ticks.y = element_blank())
      
      gh_NS4a <- g}    
    
    {    g <- ggplot(data = tep_df %>% filter(prot %in% c("NS4b")), mapping = aes(x = values, y = ..count.., color = class))
      # g <- g + geom_histogram(binwidth = 1, color = brewer.pal(9, "Set1")[2])
      # g <- g + geom_freqpoly(binwidth = 1, size = 0.5, position = "identity")
      # g <- g + scale_color_manual(labels = c("I", "II"), values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]), name = "MHC Allele Class")
      g <- g + geom_step(aes(y = ..count..), stat = "bin", binwidth = 1, position = "identity", size = 0.5, show.legend = FALSE)
      g <- g + scale_color_manual(values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]))
      
      g <- g + scale_fill_manual(labels = c("I", "II"), values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]), name = "MHC Allele Class")
      g <- g + scale_x_continuous(breaks = seq.int(50, 900, 50), limits = c(1, 250), expand = c(0,5))
      g <- g + scale_y_continuous(breaks = seq.int(0, 12, 4), limits = c(0,12))
      g <- g + ylab("") + xlab("") 
      # + ggtitle("T cell epitopes coverage of DENV proteins")
      # g <- g + facet_grid(rows = vars(prot), margins = FALSE, shrink = TRUE)
      g <- g + facet_wrap(facets = vars(prot), shrink = TRUE, scales = "free_x", nrow = 1)
      g <- g + theme_minimal() + theme(strip.background = element_rect(color = "white", fill = "grey95"), 
                                       legend.position = "bottom", 
                                       axis.title.y = element_text(margin = margin(t=50, r = 1, b = 1, l = 1), hjust = 0.5),
                                       axis.ticks = element_line(), panel.grid.minor = element_blank(),
                                       strip.text = element_text(face = "bold", size = 10), 
                                       plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
      
      gh_NS4b <- g}
    
    {    g <- ggplot(data = tep_df %>% filter(prot %in% c("NS3")), mapping = aes(x = values, y = ..count.., color = class))
      # g <- g + geom_histogram(binwidth = 1, color = brewer.pal(9, "Set1")[2])
      # g <- g + geom_freqpoly(binwidth = 1, size = 0.5, position = "identity")
      # g <- g + scale_color_manual(labels = c("I", "II"), values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]), name = "MHC Allele Class")
      g <- g + geom_step(aes(y = ..count..), stat = "bin", binwidth = 1, position = "identity", size = 0.5, show.legend = FALSE)
      g <- g + scale_color_manual(values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]))
      
      g <- g + scale_fill_manual(labels = c("I", "II"), values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]), name = "MHC Allele Class")
      g <- g + scale_x_continuous(breaks = seq.int(50, 900, 50), limits = c(1, 620), expand = c(0,5))
      g <- g + scale_y_continuous(breaks = seq.int(0, 12, 4), limits = c(0,12))
      g <- g + ylab("") + xlab("") 
      # + ggtitle("T cell epitopes coverage of DENV proteins")
      # g <- g + facet_grid(rows = vars(prot), margins = FALSE, shrink = TRUE)
      g <- g + facet_wrap(facets = vars(prot), shrink = TRUE, scales = "free_x", nrow = 1)
      g <- g + theme_minimal() + theme(strip.background = element_rect(color = "white", fill = "grey95"), 
                                       legend.position = "bottom", 
                                       axis.title.y = element_text(margin = margin(t=50, r = 1, b = 1, l = 1), hjust = 0.5),
                                       axis.ticks = element_line(), panel.grid.minor = element_blank(),
                                       strip.text = element_text(face = "bold", size = 10), 
                                       plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
      
      gh_NS3 <- g}
    
    {    g <- ggplot(data = tep_df %>% filter(prot %in% c("NS2b")), mapping = aes(x = values, y = ..count.., color = class))
      # g <- g + geom_histogram(binwidth = 1, color = brewer.pal(9, "Set1")[2])
      # g <- g + geom_freqpoly(binwidth = 1, size = 0.5, position = "identity")
      # g <- g + scale_color_manual(labels = c("I", "II"), values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]), name = "MHC Allele Class")
      g <- g + geom_step(aes(y = ..count..), stat = "bin", binwidth = 1, position = "identity", size = 0.5, show.legend = FALSE)
      g <- g + scale_color_manual(values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]))
      
      g <- g + scale_fill_manual(labels = c("I", "II"), values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]), name = "MHC Allele Class")
      g <- g + scale_x_continuous(breaks = seq.int(50, 900, 50), limits = c(1, 130), expand = c(0,5))
      g <- g + scale_y_continuous(breaks = seq.int(0, 12, 4), limits = c(0,12))
      g <- g + ylab("") + xlab("") 
      # + ggtitle("T cell epitopes coverage of DENV proteins")
      # g <- g + facet_grid(rows = vars(prot), margins = FALSE, shrink = TRUE)
      g <- g + facet_wrap(facets = vars(prot), shrink = TRUE, scales = "free_x", nrow = 1)
      g <- g + theme_minimal() + theme(strip.background = element_rect(color = "white", fill = "grey95"), 
                                       legend.position = "bottom", 
                                       axis.title.y = element_text(margin = margin(t=50, r = 1, b = 1, l = 1), hjust = 0.5),
                                       axis.ticks = element_line(), panel.grid.minor = element_blank(),
                                       strip.text = element_text(face = "bold", size = 10), 
                                       plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
      
      gh_NS2b <- g}
    
    {    g <- ggplot(data = tep_df %>% filter(prot %in% c("NS2a")), mapping = aes(x = values, y = ..count.., color = class))
      # g <- g + geom_histogram(binwidth = 1, color = brewer.pal(9, "Set1")[2])
      # g <- g + geom_freqpoly(binwidth = 1, size = 0.5, position = "identity")
      # g <- g + scale_color_manual(labels = c("I", "II"), values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]), name = "MHC Allele Class")
      g <- g + geom_step(aes(y = ..count..), stat = "bin", binwidth = 1, position = "identity", size = 0.5, show.legend = FALSE)
      g <- g + scale_color_manual(values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]))
      
      g <- g + scale_fill_manual(labels = c("I", "II"), values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]), name = "MHC Allele Class")
      g <- g + scale_x_continuous(breaks = seq.int(50, 900, 50), limits = c(1, 218), expand = c(0,5))
      g <- g + scale_y_continuous(breaks = seq.int(0, 12, 4), limits = c(0,12))
      g <- g + ylab("") + xlab("Position") 
      # + ggtitle("T cell epitopes coverage of DENV proteins")
      # g <- g + facet_grid(rows = vars(prot), margins = FALSE, shrink = TRUE)
      g <- g + facet_wrap(facets = vars(prot), shrink = TRUE, scales = "free_x", nrow = 1)
      g <- g + theme_minimal() + theme(strip.background = element_rect(color = "white", fill = "grey95"), 
                                       legend.position = "bottom", 
                                       axis.title.y = element_text(margin = margin(t=50, r = 1, b = 1, l = 1), hjust = 0.5),
                                       axis.ticks = element_line(), panel.grid.minor = element_blank(),
                                       strip.text = element_text(face = "bold", size = 10), 
                                       plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
      
      gh_NS2a <- g}
    
    {    g <- ggplot(data = tep_df %>% filter(prot %in% c("NS1")), mapping = aes(x = values, y = ..count.., color = class))
      # g <- g + geom_histogram(binwidth = 1, color = brewer.pal(9, "Set1")[2])
      # g <- g + geom_freqpoly(binwidth = 1, size = 0.5, position = "identity")
      # g <- g + scale_color_manual(labels = c("I", "II"), values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]), name = "MHC Allele Class")
      g <- g + geom_step(aes(y = ..count..), stat = "bin", binwidth = 1, position = "identity", size = 0.5, show.legend = FALSE)
      g <- g + scale_color_manual(values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]))
      
      g <- g + scale_fill_manual(labels = c("I", "II"), values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]), name = "MHC Allele Class")
      g <- g + scale_x_continuous(breaks = seq.int(50, 900, 50), limits = c(1, 352), expand = c(0,5))
      g <- g + scale_y_continuous(breaks = seq.int(0, 12, 4), limits = c(0,12))
      g <- g + ylab("") + xlab("") 
      # + ggtitle("T cell epitopes coverage of DENV proteins")
      # g <- g + facet_grid(rows = vars(prot), margins = FALSE, shrink = TRUE)
      g <- g + facet_wrap(facets = vars(prot), shrink = TRUE, scales = "free_x", nrow = 1)
      g <- g + theme_minimal() + theme(strip.background = element_rect(color = "white", fill = "grey95"), 
                                       legend.position = "bottom", 
                                       axis.title.y = element_text(margin = margin(t=50, r = 1, b = 1, l = 1), hjust = 0.5),
                                       axis.ticks = element_line(), panel.grid.minor = element_blank(),
                                       strip.text = element_text(face = "bold", size = 10), 
                                       plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
      
      gh_NS1 <- g}
    
    {    g <- ggplot(data = tep_df %>% filter(prot %in% c("C")), mapping = aes(x = values, y = ..count.., color = class))
      # g <- g + geom_histogram(binwidth = 1, color = brewer.pal(9, "Set1")[2])
      # g <- g + geom_freqpoly(binwidth = 1, size = 0.5, position = "identity")
      # g <- g + scale_color_manual(labels = c("I", "II"), values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]), name = "MHC Allele Class")
      g <- g + geom_step(aes(y = ..count..), stat = "bin", binwidth = 1, position = "identity", size = 0.5, show.legend = FALSE)
      g <- g + scale_color_manual(values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]))
      
      g <- g + scale_fill_manual(labels = c("I", "II"), values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]), name = "MHC Allele Class")
      g <- g + scale_x_continuous(breaks = seq.int(50, 900, 50), limits = c(1, 100), expand = c(0,5))
      g <- g + scale_y_continuous(breaks = seq.int(0, 12, 4), limits = c(0,12))
      g <- g + ylab("") + xlab("") 
      # + ggtitle("T cell epitopes coverage of DENV proteins")
      # g <- g + facet_grid(rows = vars(prot), margins = FALSE, shrink = TRUE)
      g <- g + facet_wrap(facets = vars(prot), shrink = TRUE, scales = "free_x", nrow = 1)
      g <- g + theme_minimal() + theme(strip.background = element_rect(color = "white", fill = "grey95"), 
                                       legend.position = "bottom", 
                                       axis.title.y = element_text(margin = margin(t=50, r = 1, b = 1, l = 1), hjust = 0.5),
                                       axis.ticks = element_line(), panel.grid.minor = element_blank(),
                                       strip.text = element_text(face = "bold", size = 10), 
                                       plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
      
      gh_C <- g}
    
    {    g <- ggplot(data = tep_df %>% filter(prot %in% c("prM")), mapping = aes(x = values, y = ..count.., color = class))
      # g <- g + geom_histogram(binwidth = 1, color = brewer.pal(9, "Set1")[2])
      # g <- g + geom_freqpoly(binwidth = 1, size = 0.5, position = "identity")
      # g <- g + scale_color_manual(labels = c("I", "II"), values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]), name = "MHC Allele Class")
      g <- g + geom_step(aes(y = ..count..), stat = "bin", binwidth = 1, position = "identity", size = 0.5, show.legend = FALSE)
      g <- g + scale_color_manual(values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]))
      
      g <- g + scale_fill_manual(labels = c("I", "II"), values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]), name = "MHC Allele Class")
      g <- g + scale_x_continuous(breaks = seq.int(50, 900, 50), limits = c(1, 166), expand = c(0,5))
      g <- g + scale_y_continuous(breaks = seq.int(0, 12, 4), limits = c(0,12))
      g <- g + ylab("") + xlab("") 
      # + ggtitle("T cell epitopes coverage of DENV proteins")
      # g <- g + facet_grid(rows = vars(prot), margins = FALSE, shrink = TRUE)
      g <- g + facet_wrap(facets = vars(prot), shrink = TRUE, scales = "free_x", nrow = 1)
      g <- g + theme_minimal() + theme(strip.background = element_rect(color = "white", fill = "grey95"), 
                                       legend.position = "bottom", 
                                       axis.title.y = element_text(margin = margin(t=50, r = 1, b = 1, l = 1), hjust = 0.5),
                                       axis.ticks = element_line(), panel.grid.minor = element_blank(),
                                       strip.text = element_text(face = "bold", size = 10), 
                                       plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                                       axis.text.y = element_blank(), axis.ticks.y = element_blank())
      
      gh_prM <- g}
    
    {    g <- ggplot(data = tep_df %>% filter(prot %in% c("E")), mapping = aes(x = values, y = ..count.., color = class))
      # g <- g + geom_histogram(binwidth = 1, color = brewer.pal(9, "Set1")[2])
      # g <- g + geom_freqpoly(binwidth = 1, size = 0.5, position = "identity")
      # g <- g + scale_color_manual(labels = c("I", "II"), values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]), name = "MHC Allele Class")
      g <- g + geom_step(aes(y = ..count..), stat = "bin", binwidth = 1, position = "identity", size = 0.5, show.legend = FALSE)
      g <- g + scale_color_manual(values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]))
      
      g <- g + scale_fill_manual(labels = c("I", "II"), values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]), name = "MHC Allele Class")
      g <- g + scale_x_continuous(breaks = seq.int(50, 900, 50), limits = c(1, 495), expand = c(0,5))
      g <- g + scale_y_continuous(breaks = seq.int(0, 12, 4), limits = c(0,12))
      g <- g + ylab("") + xlab("") 
      # + ggtitle("T cell epitopes coverage of DENV proteins")
      # g <- g + facet_grid(rows = vars(prot), margins = FALSE, shrink = TRUE)
      g <- g + facet_wrap(facets = vars(prot), shrink = TRUE, scales = "free_x", nrow = 1)
      g <- g + theme_minimal() + theme(strip.background = element_rect(color = "white", fill = "grey95"), 
                                       legend.position = "bottom", 
                                       axis.title.y = element_text(margin = margin(t=50, r = 1, b = 1, l = 1), hjust = 0.5),
                                       axis.ticks = element_line(), panel.grid.minor = element_blank(),
                                       strip.text = element_text(face = "bold", size = 10), 
                                       plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
      
      gh_E <- g}
    
    {g2 <- ggplot(data = data.frame(x = c(" Class I restricted", " Class II restricted"), y = 0), 
                  mapping = aes(x = x, y = y, color = x))
      g2 <- g2 + geom_text(aes(label = x, size = 2),  fontface = "bold") 
      g2 <- g2 + scale_color_manual(values = c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[3]))
      g2 <- g2 + theme_void() + theme(legend.position = "none")
      
      g2 <- g2 + theme(plot.margin = unit(c(-10,0,0,0), "pt")) +  theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
      )}
    
    gh_C <- gh_C + geom_hline(yintercept = 0, show.legend = FALSE, color = "grey60", size = 0.5)
    gh_E <- gh_E + geom_hline(yintercept = 0, show.legend = FALSE, color = "grey60", size = 0.5)
    gh_prM <- gh_prM + geom_hline(yintercept = 0, show.legend = FALSE, color = "grey60", size = 0.5)
    gh_NS1 <- gh_NS1 + geom_hline(yintercept = 0, show.legend = FALSE, color = "grey60", size = 0.5)
    gh_NS2a <- gh_NS2a + geom_hline(yintercept = 0, show.legend = FALSE, color = "grey60", size = 0.5)
    gh_NS2b <- gh_NS2b + geom_hline(yintercept = 0, show.legend = FALSE, color = "grey60", size = 0.5)
    gh_NS3 <- gh_NS3 + geom_hline(yintercept = 0, show.legend = FALSE, color = "grey60", size = 0.5)
    gh_NS4a <- gh_NS4a + geom_hline(yintercept = 0, show.legend = FALSE, color = "grey60", size = 0.5)
    gh_NS4b <- gh_NS4b + geom_hline(yintercept = 0, show.legend = FALSE, color = "grey60", size = 0.5)
    gh_NS5 <- gh_NS5 + geom_hline(yintercept = 0, show.legend = FALSE, color = "grey60", size = 0.5)
    
    fig1 <- ggarrange(gh_C + theme(panel.border = element_blank(), plot.margin = margin(0, 10, 0, 0)), 
                      gh_prM + theme(panel.border = element_blank(), plot.margin = margin(0, 10, 0, 0)), 
                      ncol = 2, nrow = 1, common.legend = FALSE, legend = "none")
    fig2 <- ggarrange(gh_NS2b + theme(panel.border = element_blank(), plot.margin = margin(0, 10, 0, 0)), 
                      gh_NS4a + theme(panel.border = element_blank(), plot.margin = margin(0, 10, 0, 0)), 
                      ncol = 2, nrow = 1, common.legend = FALSE, legend = "none")
    
    fig3 <- ggarrange(fig1, 
                      fig2, 
                      gh_NS4b + theme(panel.border = element_blank(), plot.margin = margin(0, 10, 0, 0)), 
                      gh_NS2a + theme(panel.border = element_blank(), plot.margin = margin(0, 10, 0, 0)), 
                      nrow = 4, ncol = 1, legend = "none", common.legend = FALSE)    
    
    fig11 <- ggarrange(gh_C + theme(panel.border = element_blank(), plot.margin = margin(0, 10, 0, 0)), 
                       gh_prM + theme(panel.border = element_blank(), plot.margin = margin(0, 10, 0, 0)),
                       gh_NS2a + theme(panel.border = element_blank(), plot.margin = margin(0, 10, 0, 0)) + 
                         theme(axis.title.x = element_text(colour = "white")) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()),
                       gh_NS2b + theme(panel.border = element_blank(), plot.margin = margin(0, 10, 0, 0)) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()),
                       gh_NS4a + theme(panel.border = element_blank(), plot.margin = margin(0, 10, 0, 0)),
                       ncol = 5, nrow = 1, widths = c(1, 1.5, 2.25, 1.25, 1.25), common.legend = FALSE, legend = "none")
    
    fig12 <- ggarrange(gh_NS1 + theme(panel.border = element_blank(), plot.margin = margin(0, 10, 0, 0)), 
                       gh_E + theme(panel.border = element_blank(), plot.margin = margin(0, 10, 0, 0)) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()),
                       ncol = 2, nrow = 1, widths = c(1, 1.4), common.legend = FALSE, legend = "none")
    
    fig13 <- ggarrange(gh_NS4b + theme(panel.border = element_blank(), plot.margin = margin(0, 10, 0, 0)), 
                       gh_NS3 + theme(panel.border = element_blank(), plot.margin = margin(0, 10, 0, 0)) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()),
                       ncol = 2, nrow = 1, widths = c(1, 2.4), common.legend = FALSE, legend = "none")
    gempty <- ggplot(df) + theme_void()
    fig_tmp <- ggarrange(gempty, 
                         fig11,
                         fig12,
                         fig13,
                         gh_NS5,
                         ggarrange(gempty, g2, gempty, ncol = 3, nrow = 1, widths = c(0.25,1,0.25), legend = FALSE),
                         ncol = 1, nrow = 6, heights = c(0.5, 2.5, 2.5, 2.5, 2.5, 1), legend = "none")
    
    fig_tmp <- annotate_figure(fig_tmp,
                               left = text_grob("Number of epitopes", face = "bold", size = 11, rot = 90))

    g_polygons_DENV1 <- fig_tmp         
    
  }  
}

fig_tmp <- ggarrange(
  grobdt_numTeps, 
  grobdt_Ep_mapping, 
  g_polygons_DENV1 + theme(plot.margin = margin(0), plot.title = element_blank()), 
  nrow = 3, ncol = 1, labels = c("A", "B", "C"), 
  font.label = list(size=14), hjust = 0, vjust = 1.5, 
  heights = c(0.75, 0.75, 2.25))

file.remove(here("Figure", "temp_fig.png"))

fig_tmp

}