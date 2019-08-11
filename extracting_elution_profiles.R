source("matrix_operations.R")

# set parameters ----------------------------------------------------------

idx_file <- 1
scan_num <- "0"
ms1_file <- "15c.ms1"
dir_data <- str_replace(ms1_file, ".ms1$", "_ms1_sparsebinned_singlescans_tannotated/")
mass_accuracy <- 6e-6

# read matrices from file -------------------------------------------------

xfiles <- list.files(dir_data) %>%
	stringr::str_subset("xdecoy\\.txt")

idx_ordered <- stringr::str_extract(xfiles, "^\\d+") %>%
	as.integer() %>%
	order()

xfiles <- xfiles[idx_ordered]
yfiles <- stringr::str_replace(xfiles, "xdecoy", "y")
annfiles <- stringr::str_replace(xfiles, "xdecoy", "annfile")

input_files <- purrr::pmap(list(xfiles, yfiles, annfiles), c) %>%
	purrr::map(~ paste0(dir_data, .x))

subX_list <- list()
subY_list <- list()
sub_binbounds_list <- list()
sub_features_list <- list()

for (i in seq_along(input_files)) {
	x <- readLines(input_files[[i]][1])
	y <- readLines(input_files[[i]][2])
	sub_features_df <- read_tsv(input_files[[i]][3], col_names = c("scan", "minmz", "maxmz", "charge"),
								skip = 1, col_types = "iddi")

	scans <- get_scannum(x)
	unique_scans <- unique(scans)
	binbounds <- get_binbounds(x)
	X <- sparse_to_regular(x)
	X <- X[, 1:(ncol(X) / 2)]
	Y <- sparse_to_regular(y)
	subX <- vector("list", length(unique_scans)) %>% set_names(unique_scans)
	subY <- vector("list", length(unique_scans)) %>% set_names(unique_scans)
	sub_binbounds <- vector("list", length(unique_scans)) %>% set_names(unique_scans)
	sub_features <- vector("list", length(unique_scans)) %>% set_names(unique_scans)

	for (j in unique_scans) {
		idx_rows <- scans[scans == j] %>% names() %>% as.integer()
		subX[[as.character(j)]] <- X[idx_rows, ]
		idx_nonzero_cols <- colSums(subX[[as.character(j)]]) > 0
		subX[[as.character(j)]] <- subX[[as.character(j)]][, idx_nonzero_cols]
		subY[[as.character(j)]] <- Y[idx_rows, , drop = FALSE]
		sub_binbounds[[as.character(j)]] <- binbounds[rownames(binbounds) %in% idx_rows, ]
		sub_features[[as.character(j)]] <- filter(sub_features_df, idx_nonzero_cols)
	}

	subX_list <- c(subX_list, subX)
	subY_list <- c(subY_list, subY)
	sub_binbounds_list <- c(sub_binbounds_list, sub_binbounds)
	sub_features_list <- c(sub_features_list, sub_features)
}

sub_features_list <- map2(
	sub_features_list, subX_list,
	~add_column(.x, isop_pattern = apply(.y, 2,
										 function(x) {idx <- which(x > 0); x[idx[1]:idx[length(idx)]]})
	)
)

coefs_list <- vector("list", length(subX_list))
for (i in seq_along(subX_list)) {
	X <- subX_list[[i]]
	Y <- subY_list[[i]]
	binbounds <- sub_binbounds_list[[i]]

	mzbins <- rowMeans(binbounds)
	idx_nonzero_y <- Y[, 1] > 0
	Y <- Y[idx_nonzero_y, , drop = FALSE]
	X <- X[idx_nonzero_y, ]
	mzbins <- mzbins[idx_nonzero_y]
	rownames(Y) <-  mzbins

	X_df <- as_tibble(X) %>%
		mutate(mzbin = mzbins) %>%
		tidyr::gather("feature", "value", -mzbin) %>%
		filter(value > 0) %>%
		mutate(feature = str_remove(feature, "^V") %>% as.integer())

	graph_df <- X_df %>%
		select(mzbin, feature) %>%
		group_by(mzbin) %>%
		nest(.key = "features") %>%
		mutate(features = map(features, ~.x[[1]])) %>%
		mutate(len = map_int(features, length)) %>%
		filter(len > 1) %>%
		mutate(features = map(features, combn, m = 2, simplify = FALSE)) %>%
		unnest(features) %>%
		mutate(
			from = map_int(features, 1),
			to = map_int(features, 2)
		) %>%
		select(mzbin, from, to)

	g <- igraph::graph_from_data_frame(graph_df[, c(2, 3)],
									   vertices = unique(X_df$feature), directed = FALSE)

	ceb <- igraph::cluster_edge_betweenness(g)

	feature_cluster <- stack(ceb) %>%
		rename(feature = values, cluster = ind) %>%
		mutate(feature = as.integer(feature))

	multi_clusters <- feature_cluster %>%
		group_by(cluster) %>%
		summarise(n_features= n()) %>%
		filter(n_features > 1) %>%
		pull(cluster) %>% as.integer()

	mzbin_cluster <- X_df %>%
		left_join(feature_cluster, "feature") %>%
		select(mzbin, feature, cluster)

	coefs <- vector("numeric", length = ncol(X))
	for (clust in multi_clusters) {
		idx_col <- filter(mzbin_cluster, cluster == clust) %>% pull(feature) %>% unique()
		idx_row <- rowSums(X[, idx_col, drop = FALSE]) > 0
		Xc <- X[idx_row, idx_col, drop = FALSE]
		Yc <- Y[idx_row, ,drop = FALSE]

		if (nrow(Xc) < 4) next()

		mdl_cv_wt <- glmnet::cv.glmnet(Xc, Yc, intercept = TRUE, lower.limits = 0,
									   nfolds = nrow(Xc), grouped = FALSE,
									   weights = 1 / (10 + drop(Yc)))

		idx_col_sub <- predict(mdl_cv_wt, type = "coef", s = "lambda.min") %>%
			as.vector() %>% (function(x) which(x[-1] > 0))

		if (length(idx_col_sub) == 0) next()
		mdl <- glmnet::glmnet(Xc[, idx_col_sub, drop = FALSE], Yc, lambda = 0,
							  intercept = FALSE, lower.limits = 0)

		coefs_sub <- predict(mdl, type = "coef") %>% as.vector() %>% .[-1]
		coefs[idx_col[idx_col_sub]] <- coefs_sub
	}

	coefs_list[[i]] <- coefs
}

sub_features_list <- map2(sub_features_list, coefs_list,
	 ~add_column(.x, peak = .y))

sub_features_list <- map(sub_features_list, filter, peak > 0) %>%
	map(~ mutate(.x, scan_peak = map2(scan, peak, ~list(c(.x, .y))))) %>%
	map(select, minmz, maxmz, charge, isop_pattern, scan_peak)

features <- sub_features_list %>%
	map(~ mutate(.x, meanmz = ((minmz + maxmz) / 2))) %>%
	map(select, meanmz, charge, isop_pattern, scan_peak) %>%
	bind_rows() %>%
	pmap(list)

features_aligned <- features[1]
for (feature_curr in features[-1]) {
	add_new <- TRUE
	for (i in seq_along(features_aligned)) {
		feature <- features_aligned[[i]]
		if (
			feature_curr$charge == feature$charge &&
			between(feature_curr$meanmz, feature$meanmz * (1 - mass_accuracy), feature$meanmz * (1 + mass_accuracy)) &&
			identical(feature_curr$isop_pattern, feature$isop_pattern)
		) {
			feature$scan_peak <- c(feature$scan_peak, feature_curr$scan_peak)
			add_new <- FALSE
		}
		features_aligned[[i]] <- feature
	}

	if (add_new) features_aligned <- c(features_aligned, list(feature_curr))
}

idx <-  map_int(features_aligned, ~length(.x$scan_peak)) %>%
	(function(x) x > 5)

features_aligned2 <- features_aligned[idx]
idx_ord <- map_dbl(features_aligned2, "meanmz") %>%
	order()

features_aligned2 <- features_aligned2[idx_ord]

profiles_mat <- matrix(NA_real_, ncol = length(features_aligned2), nrow = length(subX_list))
for (i in seq_along(features_aligned2)) {
	peaks <- features_aligned2[[i]][["scan_peak"]]
	for (peak in peaks) {
		profiles_mat[as.integer(peak[1]) + 1, i] <- peak[2]
	}
}