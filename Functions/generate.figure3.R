generate.figure3 <- function(ProteinDatatableOut, ...) {

  df <- ProteinDatatableOut %>% filter(order3 == "T") 
  
  ceil <- ceiling(df$meanID %>% length() / 1)
  ceil <- rep(seq.int(1,ceil), each=1)
  
  df %>% 
    arrange(meanID) %>% 
    mutate(facetID = ceil) %>% 
    select(`IEDB ID`, `Epitope Sequence`, Protein, 
           IdD1, IdD2, IdD3, IdD4, 
           meanID, meanIDrank, facetID) -> df
  
  df %>% pull(`Epitope Sequence`) %>% 
    factor(levels = df$`Epitope Sequence`) -> df$`Epitope Sequence`
  
  df %>% gather(key = IdDkey, value = IdDvalue, c(IdD1, IdD2, IdD3, IdD4)) -> df2
  
  df2$IdDkey <- df2 %>% pull(IdDkey) %>% 
    factor(levels = c("IdD1", "IdD2", "IdD3", "IdD4"), 
           labels = c("DENV1", "DENV2", "DENV3", "DENV4"))

  colorings <- c(brewer.pal(12, "Paired")[1], brewer.pal(9, "Set1")[c(5,4)], brewer.pal(12, "Paired")[c(5, 12)], 
                 brewer.pal(9, "Set1")[2], brewer.pal(12, "Paired")[9], brewer.pal(9, "Set1")[9], 
                 brewer.pal(9, "Set1")[1], brewer.pal(12, "Paired")[7])
  
  names(colorings) <- c("C", "E", "NS1", "NS2a", "NS2b", "NS3", "NS4a", "NS4b", "NS5", "prM")
  
  # Split-tiles
  
  splitings <- max(df2$facetID)/3 * c(1,2,3)
  
  df2$splitings <- if_else(condition = (df2$facetID >= splitings[2]), true = 1, 
                           false = if_else(condition = (df2$facetID >= splitings[1]) & (df2$facetID < splitings[2]), 
                                           true = 2, false = 3))
  
  suppressWarnings({
    
    gf <- ggplot(data = df2, mapping = aes(x = IdDkey, y = `Epitope Sequence`, fill = IdDvalue))
    
    gf <- gf + geom_tile(color = "grey60", width = 0.40, height = 0.60)
    
    gf <- gf + scale_fill_gradient(low = "white",# brewer.pal(9, "PuBu")[1],
                                   high = brewer.pal(9, "PuBu")[9],
                                   guide = "colorbar", 
                                   limits = c(0,1),
                                   breaks = c(0, 0.5, 1),
                                   labels = c("0", "0.5", "1"),
                                   name = "Fraction of sequences") 
    
    
    gf <- gf + geom_rect(xmin = -2, xmax = 0, ymin = 0, ymax = nrow(df2) + 1, fill = "white") + expand_limits(x = -2)
    
    gf <- gf + geom_text(aes(color = Protein, label = `Epitope Sequence`, angle = 0, hjust = 0), size = 2, fontface = "plain", x = -2, fill = NA, show.legend = TRUE) 
    
    gf <- gf + scale_colour_manual(values = colorings, name = "Protein")
    
    gf <- gf + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), panel.background = element_rect(fill = "white"))
    
    gf <- gf + ylab("") + xlab("") + scale_x_discrete(position = "top")
    
    gf <- gf + theme(axis.text.x = element_text(angle = 90, size = 8, hjust = 0.5,  vjust = 0.5, color = "grey5", face = "bold"), axis.ticks.x = element_blank())
    
    gf <- gf + theme(panel.grid = element_blank(),
                     strip.background = element_blank(),
                     strip.text = element_blank(), axis.ticks.y = element_blank(),
                     axis.text.y = element_blank(),
                     plot.title = element_text(size = 15, face = "bold", hjust = 1, vjust = 0.5),
                     plot.subtitle = element_blank(),
                     axis.title.y = element_text(face = "bold", hjust = 0.5),
                     axis.title.x = element_blank(),
                     legend.title = element_text(size = 9),
                     legend.text = element_text(size = 8), legend.spacing.x = unit(5, "mm"),
                     legend.key = element_rect(fill = "white", color = "white"),
                     legend.key.size = unit(x = 10, units = "pt"), legend.box = "vertical",
                     legend.position = "bottom", legend.spacing.y = unit(-1, "pt"))
    
    gf <- gf + theme(panel.spacing.x = unit(0, "lines")) +
      guides(fill=guide_colorbar(title.vjust = 0.85, 
                                 label.theme = element_text(size = 8, margin = margin(t=2)))) +
      guides(color=guide_legend(override.aes = 
                                  list(labels = c(paste0("E (", {df %>% filter(Protein == "E") %>% n_distinct()}, ")"), 
                                                  paste0("NS1 (", {df %>% filter(Protein == "NS1") %>% n_distinct()}, ")"),
                                                  paste0("NS3 (", {df %>% filter(Protein == "NS3") %>% n_distinct()}, ")"),
                                                  paste0("NS4b (", {df %>% filter(Protein == "NS4b") %>% n_distinct()}, ")"),
                                                  paste0("NS5 (", {df %>% filter(Protein == "NS5") %>% n_distinct()}, ")")), 
                                       size = 2.25), 
                                label = FALSE, title = NULL)) 
    
    gempty <- ggplot(data = df) + theme_void()
    
    gf <- gf + facet_wrap(facets = vars(splitings), scales = "free")
      
  })
  
  fig_tmp <- gf
  
  
fig_tmp

}