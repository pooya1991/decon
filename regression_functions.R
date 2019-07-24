source("matrix_operations.R")

glmnet_ols <- function(xfile, yfile) {
	X <- sparse_to_regular(readLines(xfile))
	Y <- sparse_to_regular(readLines(yfile))
	mdl <- glmnet::glmnet(X, Y, lambda = 0, standardize = TRUE, lower.limits = 0, thresh = 1e-3, intercept = FALSE)
	as.matrix(mdl$beta) %>%
		regular_to_sparse()
}

glmnet_batch <- function(x, y, ...) {
	scans <- unique(get_scannum(x)) %>% as.character()
	lookup_row_ranges <- get_scan_row_ranges(x)
	lookup_col_ranges <- get_scan_nonzero_col_ranges(x)
	X <- sparse_to_regular(x)
	Y <- sparse_to_regular(y)
	subX <- vector("list", length(scans))
	subY <- vector("list", length(scans))
	names(subX) <- scans
	names(subY) <- scans

	for (i in scans) {
		row_range <- lookup_row_ranges[i, ]
		col_range <- lookup_col_ranges[i, ]
		# subX[[i]] <- X[row_range[1]:row_range[2], col_range[1]:col_range[2]]
		subX[[i]] <- X[row_range[1]:row_range[2], ]
		subY[[i]] <- Y[row_range[1]:row_range[2], , drop = FALSE]
	}

	cv_glmnet_partial <- purrr::partial(glmnet::cv.glmnet, ...)
	mdl_list <- map2(subX, subY, ~cv_glmnet_partial(x = .x, y = .y))
	map2(mdl_list, subX, ~predict(.x, .y, s ="lambda.min", type = "coef"))
}
