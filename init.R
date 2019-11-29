## R Packages
# Checks and installs the required R packages.

list_packages <-  c("here", "rmarkdown", "RColorBrewer", "tidyverse", "ggpubr", "seqinr",
                    "readxl",  "readr", "DT", "htmltools", "webshot", "magick", "grid",
                    "gridExtra", "formattable", "glue", "maps", "reticulate", "BALCONY", "BiocManager")

new_packages <- list_packages[!(list_packages %in% installed.packages()[,"Package"])]

if(length(new_packages)) install.packages(new_packages)

list_packages <-  c("msa", "Biostrings")

new_packages <- list_packages[!(list_packages %in% installed.packages()[,"Package"])]

if(length(new_packages)){
  BiocManager::install()
  BiocManager::install(c("msa", "Biostrings"))
}

if(is.null(webshot:::find_phantom())){
  webshot::install_phantomjs()
}


