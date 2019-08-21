regressor <- function(X, Y, clustering_info) {
	valid_clusters <- distinct(clustering_info[c("feature", "cluster")]) %>%
		count(cluster) %>%
		filter(n > 1) %>% pull(cluster) %>% as.integer()

	coefs <- vector("numeric", length = ncol(X))
	for (clust in valid_clusters) {
		idx_col <- filter(clustering_info, cluster == clust) %>% pull(feature) %>% unique()
		idx_row <- rowSums(X[, idx_col, drop = FALSE]) > 0
		Xc <- X[idx_row, idx_col, drop = FALSE]
		Yc <- Y[idx_row, , drop = FALSE]

		if (nrow(Xc) < 4) next()
		idx_col_sub <- feature_selector(Xc, Yc)
		if (length(idx_col_sub) == 0) next()

		coefs_sub <- regressor_per_clust(Xc[, idx_col_sub, drop = FALSE], Yc)
		coefs[idx_col[idx_col_sub]] <- coefs_sub
	}
	coefs
}

feature_selector <- function(X, Y) {
	mdl <- glmnet::cv.glmnet(X, Y, intercept = TRUE, lower.limits = 0,
							 nfolds = nrow(X), grouped = FALSE, type.measure = "mae",
							 weights = 1 / (10 + drop(Y)))

	idx_col<- predict(mdl, type = "coef", s = "lambda.min") %>%
		as.vector() %>% (function(x) which(x[-1] > 0))
	idx_col
}

regressor_per_clust <- function(X, Y) {
	mdl <- glmnet::glmnet(X, Y, intercept = FALSE, lower.limits = 0, lambda = 0)

	coefs <- predict(mdl, type = "coef") %>% as.vector() %>% .[-1]
	coefs
}