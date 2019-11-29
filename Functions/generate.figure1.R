generate.figure1 <- function(ProteinNumSeqs, NameProteins, ProteinTypes, ...) {

x <- matrix(nrow = length(ProteinTypes), 
            ncol = length(NameProteins))

for (i in seq.int(1, length(NameProteins))) {
  x[,i] <- ProteinNumSeqs[[i]]
}

x <- rbind(x, colSums(x))
x <- cbind(x, rowSums(x))

numSeqs <- as.data.frame(x = x)
colnames(numSeqs) <- c(NameProteins, "Total")
numSeqs <- rownames_to_column(.data = as.data.frame(t(numSeqs)))
colnames(numSeqs) <- c("Protein", unlist(ProteinTypes), "All serotypes")

list.files(here("Data"), pattern = "outMSA") -> z
year <- vector()
country_names <- vector()
for (k in seq.int(1, length(z))) {
  aln <- read.alignment(file = here("Data", z[k]), format = "fasta")
  country_names <- c(country_names, 
                     str_split(string = aln$nam, pattern = "\\|", simplify = TRUE)[,4])
  year <- c(year, 
            str_split(string = aln$nam, pattern = "\\|", simplify = TRUE)[,3])
}

year[year == "NA"] <- NA
year <- na.omit(year)

country_names[country_names == "NA"] <- NA
country_names  %>% na.omit(country_names) %>% 
  str_replace_all(pattern = "_", replacement = " ") -> country_names

country_names <- data.frame(country_names, stringsAsFactors = FALSE)

df_seq_countries <- country_names %>% group_by(country_names) %>% summarise(n=n())


df_map <- map_data("world")

diff <- setdiff(df_seq_countries$country_names, df_map$region)

# diff
# "Antigua and Barbuda"              "Borneo"                           "British Virgin Islands"          
# "Cote dIvoire"                     "East Timor"                       "Pacific Ocean"                   
# "Saint Kitts and Nevis"            "Saint Vincent and the Grenadines" "Timor Leste"                     
# "Trinidad and Tobago"              "Tuvalu"                           "Viet Nam"              

df_seq_countries <- rbind(df_seq_countries, c("Antigua", df_seq_countries %>% 
                                                  filter(country_names == diff[1]) %>% pull(n)))

df_seq_countries <- rbind(df_seq_countries, c("Barbuda", df_seq_countries %>% 
                                                filter(country_names == diff[1]) %>% pull(n)))

df_seq_countries <- rbind(df_seq_countries, c("Virgin Islands", df_seq_countries %>% 
                                                filter(country_names == diff[3]) %>% pull(n)))

df_seq_countries <- rbind(df_seq_countries, c("Ivory Coast", df_seq_countries %>% 
                                                filter(country_names == diff[4]) %>% pull(n)))

df_seq_countries <- rbind(df_seq_countries, c("Timor-Leste", df_seq_countries %>% 
                                                filter(country_names == diff[5]) %>% pull(n)))

df_seq_countries <- rbind(df_seq_countries, c("Saint Kitts", df_seq_countries %>% 
                                                filter(country_names == diff[7]) %>% pull(n)))

df_seq_countries <- rbind(df_seq_countries, c("Nevis", df_seq_countries %>% 
                                                filter(country_names == diff[7]) %>% pull(n)))

df_seq_countries <- rbind(df_seq_countries, c("Saint Vincent", df_seq_countries %>% 
                                                filter(country_names == diff[8]) %>% pull(n)))

df_seq_countries <- rbind(df_seq_countries, c("Grenadines", df_seq_countries %>% 
                                                filter(country_names == diff[8]) %>% pull(n)))

df_seq_countries <- rbind(df_seq_countries, c("Timor-Leste", df_seq_countries %>% 
                                                filter(country_names == diff[9]) %>% pull(n)))

df_seq_countries <- rbind(df_seq_countries, c("Trinidad", df_seq_countries %>% 
                                                filter(country_names == diff[10]) %>% pull(n)))

df_seq_countries <- rbind(df_seq_countries, c("Tobago", df_seq_countries %>% 
                                                filter(country_names == diff[10]) %>% pull(n)))

df_seq_countries <- rbind(df_seq_countries, c("Vietnam", df_seq_countries %>% 
                                                filter(country_names == diff[12]) %>% pull(n)))

df_seq_countries <- df_seq_countries %>% filter(!(country_names %in% diff))

df_seq_countries$n <- as.double(df_seq_countries$n)

df_seq_countries <- df_seq_countries %>% group_by(country_names) %>% summarise_all(.funs=sum)

df_seq_countries$region <- df_seq_countries$country_names
df_seq_countries$samples <- df_seq_countries$n

df <- left_join(df_map, df_seq_countries, "region")


g <- ggplot(df, aes(x = long, y = lat, label = region)) + 
  geom_polygon(aes( group = group, fill = samples))

g <- g + scale_fill_gradient(low = brewer.pal(9, "YlOrRd")[2], 
                             high = brewer.pal(9, "YlOrRd")[9], 
                             na.value = "grey85", 
                             guide = "colorbar")

g <- g + theme(legend.position = "bottom", 
               panel.background = element_blank(), 
               panel.grid = element_blank(),
               axis.title = element_blank(), 
               axis.text = element_blank(), 
               axis.ticks = element_blank())

g <- g + guides(fill = guide_colorbar(title = "Number of DENV protein sequences", 
                                      title.theme = element_text(size = 10, hjust = 0.5), title.position = "top", 
                                      label.theme = element_text(size = 8), nbin = 5000, 
                                      barwidth = 12, barheight = 0.85, raster = TRUE, draw.llim = FALSE, draw.ulim = FALSE))

dt_numSeqs <- formattable(numSeqs, 
                          align = "c", 
                          list(area(col = 2:5, row = 1:10) ~ custom_color_tile(brewer.pal(9, "YlOrRd")[2], 
                                                                               brewer.pal(8, "YlOrRd")[4]), 
                               area(col = c(1,6)) ~ custom_title(), area(row = 11) ~ custom_title()))

heading <- formatter("span", style = ~style(font.weight = "bold", font.size = "2.0em"))
colnames(dt_numSeqs) <- heading(names(numSeqs))

export_formattable(dt_numSeqs, here("Figure", "temp_fig.png"), 10, width = "100%", height = "100%")

grobdt_numSeqs <- rasterGrob(image_read(here("Figure", "temp_fig.png")), just = "center")

fig_tmp <- ggarrange(grobdt_numSeqs, g,
                     nrow = 2, ncol = 1, 
                     labels=c("A", "B"), heights = c(0.8, 1),
                     font.label = list(size=14), hjust = 0, vjust = 2)

file.remove(here("Figure", "temp_fig.png"))

fig_tmp

}