#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
library(magrittr)
source("matrix_operations.R")

# ydir <- "15b_ms1_sparsebinned_singlescans_tannotated/"
# outdir <- "15b_peptide_abundances/"
ydir <- args[1L]
xdir <- ydir
outdir <- args[2L]

xfiles <- list.files(xdir) %>%
	stringr::str_subset("xdecoy\\.txt")

idx_ordered <- stringr::str_extract(xfiles, "^\\d+") %>%
	as.integer() %>%
	order()

xfiles <- xfiles[idx_ordered]
yfiles <- stringr::str_replace(xfiles, "xdecoy", "y")
t <- as.integer(stringr::str_extract(xfiles, "^\\d+"))

if (!dir.exists(outdir)) {
	dir.create(outdir)
}

numjobs <- 0

input_files <- purrr::map2(xfiles, yfiles, c) %>%
	purrr::map(~ paste0(xdir, .x))

for (i in seq_along(input_files)) {
	outputdir <- paste0(outdir, t[i], '/')
	if (!dir.exists(outputdir)) dir.create(outputdir)
	cat("Computing OLS on", xfiles[i], yfiles[i], "\n")
	bB <- glmnet_ols(input_files[[i]][1], input_files[[i]][2])
	write(bB, paste0(outputdir, "bB.txt"))
	numjobs <- numjobs + 1
}

cat(numjobs, "\n")