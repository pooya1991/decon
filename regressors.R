regressor <- function(X, Y, clustering_info, drop_zero_rows = TRUE) {
	idx_nonzero_y <- Y[, 1] > 0
	Y <- Y[idx_nonzero_y, , drop = FALSE]
	X <- X[idx_nonzero_y, ]

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
		browser()
		if (length(idx_col_sub) == 0) next()
		Xc <- Xc[, idx_col_sub, drop = FALSE]

		if (drop_zero_rows) {
			idx_row_sub <- rowSums(Xc) > 0
			Xc <- Xc[idx_row_sub, , drop = FALSE]
			Yc <- Yc[idx_row_sub, , drop = FALSE]
		}

		coefs_sub <- regressor_per_clust(Xc, Yc)
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
	idx_col <- 1:ncol(X)
	idx_col
}

regressor_per_clust <- function(X, Y) {
	mdl <- glmnet::glmnet(X, Y, intercept = FALSE, lower.limits = 0,
						  weights = 1 / (0.05 * max(Y) + drop(Y)),
						  lambda = 0)

	coefs <- predict(mdl, type = "coef") %>% as.vector() %>% .[-1]
	coefs
}

# regressor_per_clust <- function(X, Y) {
# 	beta <- CVXR::Variable(ncol(X))
# 	objective <- CVXR::Minimize(sum(CVXR::huber(Y - X %*% beta, M = 0.01)))
# 	problem <- CVXR::Problem(objective)
# 	result <- solve(problem, solver = "SCS")
# 	result$getValue(beta)
# }
#
# regressor_per_clust <- function(X, Y) {
# 	mdl <- lsfit(X, Y, intercept = FALSE)
# 	mdl$coefficients
# }
#
# regressor_per_clust <- function(X, Y) {
# 	mdl <- MASS::rlm(X, Y, intercept = FALSE)
# 	mdl$coefficients
# }
