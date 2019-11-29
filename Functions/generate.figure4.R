generate.figure4 <- function(ProteinDatatableOut, ...) {

    df <- ProteinDatatableOut %>% 
      filter(Protein == "NS3") %>% filter(order3 == "T")
        
    {
         t_ind_NS3_1 = vector()
        for (i in 1:length(df$`IEDB ID`)) {
          if (!is.na(df$D1start[i])){
            t_ind_NS3_1 = c(t_ind_NS3_1, seq.int(
              from = df$D1start[i], 
              to = df$D1end[i]))}
        }
        t_ind_NS3_2 = vector()
        for (i in 1:length(df$`IEDB ID`)) {
          if (!is.na(df$D2start[i])){
            t_ind_NS3_2 = c(t_ind_NS3_2, seq.int(
              from = df$D2start[i], 
              to = df$D2end[i]))}
        }
        t_ind_NS3_3 = vector()
        for (i in 1:length(df$`IEDB ID`)) {
          if (!is.na(df$D3start[i])){
            t_ind_NS3_3 = c(t_ind_NS3_3, seq.int(
              from = df$D3start[i], 
              to = df$D3end[i]))}
        }
        t_ind_NS3_4 = vector()
        for (i in 1:length(df$`IEDB ID`)) {
          if (!is.na(df$D4start[i])){
            t_ind_NS3_4 = c(t_ind_NS3_4, seq.int(
              from = df$D4start[i], 
              to = df$D4end[i]))}
        }
      
    tep_df <- data.frame(
        values = c(t_ind_NS3_1, t_ind_NS3_2, t_ind_NS3_3, t_ind_NS3_4),
        serotype = c(rep("DENV1", length(t_ind_NS3_1)), rep("DENV2", length(t_ind_NS3_2)),
                     rep("DENV3", length(t_ind_NS3_3)), rep("DENV4", length(t_ind_NS3_4))))
      
      g <- ggplot(data = tep_df %>% filter(serotype == "DENV4"), mapping = aes(x = values, y = ..count.., color = serotype))
      
      g <- g + geom_step(aes(y = ..count..), stat = "bin", binwidth = 1, position = "identity", size = 0.5, show.legend = FALSE)
      g <- g + scale_color_manual(values = c(brewer.pal(9, "Set1")[3]))
    
      g <- g + scale_x_continuous(limit = c(1,620), breaks = seq.int(100,600,100), expand=c(0,0))
      g <- g + scale_y_continuous(breaks = seq.int(1,8), limits = c(0,8.5), expand=c(0,0))
      g <- g + ylab("No. of epitopes") + xlab("Amino acid position") + ggtitle("")
      
      g <- g + theme_bw() + theme(strip.background = element_blank(), 
                                  legend.position = "none", 
                                  axis.ticks = element_line(), 
                                  panel.grid.minor = element_blank(),
                                  strip.text = element_blank(), 
                                  plot.title = element_blank(),
                                  axis.title = element_text(size = 10),
                                  axis.text = element_text(size = 8),
                                  axis.line = element_blank()) + theme(panel.border = element_blank())
      g_NS3 <- g + geom_hline(yintercept = 0, show.legend = FALSE, color = "grey60", size = 0.5)
      
    }

    df <- ProteinDatatableOut %>% 
      filter(Protein == "NS5") %>% filter(order3 == "T")
    
    {
     t_ind_NS5_1 = vector()
    for (i in 1:length(df$`IEDB ID`)) {
      if (!is.na(df$D1start[i])){
        t_ind_NS5_1 = c(t_ind_NS5_1, seq.int(
          from = df$D1start[i], 
          to = df$D1end[i]))}
    }
    t_ind_NS5_2 = vector()
    for (i in 1:length(df$`IEDB ID`)) {
      if (!is.na(df$D2start[i])){
        t_ind_NS5_2 = c(t_ind_NS5_2, seq.int(
          from = df$D2start[i], 
          to = df$D2end[i]))}
    }
    t_ind_NS5_3 = vector()
    for (i in 1:length(df$`IEDB ID`)) {
      if (!is.na(df$D3start[i])){
        t_ind_NS5_3 = c(t_ind_NS5_3, seq.int(
          from = df$D3start[i], 
          to = df$D3end[i]))}
    }
    t_ind_NS5_4 = vector()
    for (i in 1:length(df$`IEDB ID`)) {
      if (!is.na(df$D4start[i])){
        t_ind_NS5_4 = c(t_ind_NS5_4, seq.int(
          from = df$D4start[i], 
          to = df$D4end[i]))}
    }
    
    tep_df <- data.frame(
        values = c(t_ind_NS5_1, t_ind_NS5_2, t_ind_NS5_3, t_ind_NS5_4),
        serotype = c(rep("DENV1", length(t_ind_NS5_1)), rep("DENV2", length(t_ind_NS5_2)),
                     rep("DENV3", length(t_ind_NS5_3)), rep("DENV4", length(t_ind_NS5_4))))
      
      g <- ggplot(data = tep_df %>% filter(serotype == "DENV3"), mapping = aes(x = values, y = ..count.., color = serotype))

      g <- g + geom_step(aes(y = ..count..), stat = "bin", binwidth = 1, position = "identity", size = 0.5, show.legend = FALSE)
      g <- g + scale_color_manual(values = c(brewer.pal(9, "Set1")[3]))
      
      g <- g + scale_x_continuous(limit = c(1,900), breaks = seq.int(150,900,150), expand=c(0,0))
      g <- g + scale_y_continuous(breaks = seq.int(1,8), limits = c(0,8.5), expand=c(0,0))
      g <- g + ylab("No. of epitopes") + xlab("Amino acid position") + ggtitle("")
    
      g <- g + theme_bw() + theme(strip.background = element_blank(), 
                                  legend.position = "none", 
                                  axis.ticks = element_line(), 
                                  panel.grid.minor = element_blank(),
                                  strip.text = element_blank(), 
                                  plot.title = element_blank(),
                                  axis.title = element_text(size = 10),
                                  axis.text = element_text(size = 8),
                                  axis.line = element_blank()) + theme(panel.border = element_blank())
      g_NS5 <- g + geom_hline(yintercept = 0, show.legend = FALSE, color = "grey60", size = 0.5)
      
    } 
    

    {
      grobHelicase_0 <- image_crop(image_read(
        here("Figure/Utils", "NS3_Helicase_RNA_0deg.png")), 
        "2400x2350+100+255")
      
      grobHelicase_180Y <- image_crop(image_read(
        here("Figure/Utils", "NS3_Helicase_RNA_180deg.png")), 
        "2400x2350+200+255")
      
      grobHelicase_0 <- rasterGrob(grobHelicase_0, just = "center")
      
      grobHelicase_180Y <- rasterGrob(grobHelicase_180Y, just = "center")
      
      grob180 <- rasterGrob(image_read(
        here("Figure/Utils", "180_rot.tif")), just = "center")
      
      grobPolymerase_0 <- image_crop(image_read(
        here("Figure/Utils", "NS5_0deg.png")), 
        "2400x2350+100+200")
      
      grobPolymerase_180Y <- image_crop(image_read(
        here("Figure/Utils", "NS5_180deg.png")), 
        "2400x2350+200+200")
      
      grobPolymerase_0 <- rasterGrob(grobPolymerase_0, just = "center")
      
      grobPolymerase_180Y <- rasterGrob(grobPolymerase_180Y, just = "center")
      
      grob180 <- rasterGrob(image_read(here("Figure/Utils", "180_rot.tif")), just = "center")
      
      g_prot_NS3 <- ggarrange(grobHelicase_0, 
                              grob180,
                              grobHelicase_180Y,
                              nrow = 1, ncol = 3,
                              widths = c(1,0.2,1), align = "hv") + 
        theme(plot.margin = unit(c(0,0,0,0), "pt")) + 
        theme(panel.grid = element_blank()) + 
        theme(axis.title = element_blank()) + 
        theme(axis.text = element_blank()) + 
        theme(axis.ticks = element_blank())
      
      
      g_prot_NS5 <- ggarrange(grobPolymerase_0, 
                              grob180,
                              grobPolymerase_180Y,
                              nrow = 1, ncol = 3,
                              widths = c(1,0.2,1), align = "hv") + 
        theme(plot.margin = unit(c(0,0,0,0), "pt")) + 
        theme(panel.grid = element_blank()) + 
        theme(axis.title = element_blank()) + 
        theme(axis.text = element_blank()) + 
        theme(axis.ticks = element_blank())
      
      gempty <- ggplot(data = tep_df) + theme_void()
      
      fig_NS5 <- ggarrange(g_NS5,
                           g_prot_NS5,
                           ncol = 2, nrow = 1, widths = c(0.5, 1), align = "hv", legend = "none",
                           labels = c("", ""), 
                           font.label = list(font = "plain", size = 10))

      fig_NS5 <- annotate_figure(fig_NS5, 
                                 top = text_grob("DENV NS5 protien", color = "black", face = "bold", size = 12, just = "center"))
      
      fig_NS3 <- ggarrange(g_NS3,
                           g_prot_NS3,
                           ncol = 2, nrow = 1, widths = c(0.5, 1), align = "hv", legend = "none",
                           labels = c("", ""), 
                           font.label = list(font = "plain", size = 10))
      
      fig_NS3 <- annotate_figure(fig_NS3, 
                                  top = text_grob("DENV NS3 protien", color = "black", face = "bold", size = 12, just = "center"))
      
      fig_tmp <- ggarrange(fig_NS5,
                           fig_NS3,
                           ncol = 1, nrow = 2, heights = c(1, 1), align = "hv", legend = "none",
                           labels = c("A", "B"),
                           font.label = list(font = "plain", size = 10))
    }
      
fig_tmp

}