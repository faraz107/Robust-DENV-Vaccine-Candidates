

list_packages <-  c("here", "rmarkdown", "RColorBrewer", "tidyverse", "ggpubr", "seqinr",
                    "readxl",  "readr", "DT", "htmltools", "webshot", "magick", "grid",
                    "gridExtra", "formattable", "glue", "maps", "reticulate", "BALCONY", "BiocManager")


install.packages(list_packages)

BiocManager::install()
BiocManager::install(c("msa", "Biostrings"))

webshot::install_phantomjs()