if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.14")

BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")

# Install the optparse package, which allows better commandline arg parsing
install.packages('optparse',dependencies=TRUE, repos='http://cran.rstudio.com/')
