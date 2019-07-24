source("regression_functions.R")

ydir <- "15b_ms1_sparsebinned_singlescans_tannotated/"
xdir <- ydir
outdir <- "15b_peptide_abundances/"

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

subX_list <- list()
subY_list <- list()
sub_binbounds_list <- list()

for (i in seq_along(input_files)) {
	x <- readLines(input_files[[i]][1])
	y <- readLines(input_files[[i]][2])
	scans <- get_scannum(x)
	unique_scans <- unique(scans)
	binbounds <- get_binbounds(x)
	X <- sparse_to_regular(x)
	Y <- sparse_to_regular(y)
	subX <- vector("list", length(unique_scans)) %>% set_names(unique_scans)
	subY <- vector("list", length(unique_scans)) %>% set_names(unique_scans)
	sub_binbounds <- vector("list", length(unique_scans)) %>% set_names(unique_scans)

	for (j in unique_scans) {
		idx_rows <- scans[scans == j] %>% names() %>% as.integer()
		subX[[as.character(j)]] <- X[idx_rows, ]
		subY[[as.character(j)]] <- Y[idx_rows, , drop = FALSE]
		sub_binbounds[[as.character(j)]] <- binbounds[rownames(binbounds) %in% idx_rows, ]
	}

	subX_list <- c(subX_list, subX)
	subY_list <- c(subY_list, subY)
	sub_binbounds_list <- c(sub_binbounds_list, sub_binbounds)
}

mdl_list <- map2(subX_list, subY_list,
				 ~glmnet::cv.glmnet(x = .x, y = .y, intercept = FALSE, lower.limits = 0))

coef_list <- map2(mdl_list, subX_list, ~predict(.x, .y, s ="lambda.min", type = "coef")) %>%
	map(as.vector) %>%
	map(~.x[-1])

output_mats <- map2(subX_list, coef_list, ~ t(apply(.x, 1, function(x) x * .y)))

save(subX_list, subY_list, sub_binbounds_list, coef_list, output_mats,
	 file = "regression_outputs.RData")


# generate user friendly outputs ------------------------------------------

sampleX <- subX_list[[1]]
idx_nonzero_cols <- apply(sampleX, 2, function(x) any(x > 0)) %>% which()
sampleX <- sampleX[, idx_nonzero_cols]
sampleY <- subY_list[[1]]
sample_coef <- coef_list[[1]][idx_nonzero_cols]
sample_output <- output_mats[[1]][, idx_nonzero_cols]
