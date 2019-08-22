peaklist_to_clusters <- function(peaks_df, mass_accuracy, .sort = FALSE) {
	if (.sort == TRUE) {
		peaks_df <- arrange(peaks_df, mz)
	}

	peaks_list <- pmap(peaks_df, list)
	clusters <- vector("list", length(peaks_list))
	clusters[[1]] <- list(
		minmz = peaks_list[[1]]$mz,
		maxmz = peaks_list[[1]]$mz,
		intensities = peaks_list[[1]]$intensity
	)

	i <- 1
	for (peak in peaks_list[-1]) {
		mz <- peak$mz; intensity <- peak$intensity
		clust_curr <- clusters[[i]]
		if (clust_curr$maxmz > mz * (1 - mass_accuracy)) {
			clust_curr$maxmz <- mz
			clust_curr$intensities <- c(clust_curr$intensities, intensity)
			clusters[[i]] <- clust_curr
		} else {
			i <- i + 1
			clusters[[i]] <- list(minmz = mz, maxmz = mz, intensities = intensity)
		}
	}

	compact(clusters)
}

clusters_to_features <- function(clusters, maxcharge, mass_accuracy) {
	features <- vector("list", maxcharge * length(clusters))
	for (i in seq_along(clusters)) {
		clust_curr <- clusters[[i]]
		isops <- map(maxcharge:1L, ~ list(
			minmz = clust_curr$minmz + 1.003355 / .x,
			maxmz = clust_curr$maxmz + 1.003355 / .x,
			charge = .x
		))

		for (isop in isops) {

			for (clust_next in clusters[-1:-i]) {
				if (clust_next$maxmz * (1 + mass_accuracy) < isop$minmz) next()
				if (isop$maxmz * (1 + mass_accuracy) < clust_next$minmz) break()
				# There is a match. Make this feature and add it to a list
				features[[(i - 1) * maxcharge + isop$charge]] <- list(
					charge = isop$charge,
					minmz  = clust_curr$minmz * (1 - mass_accuracy),
					maxmz  = clust_curr$maxmz * (1 + mass_accuracy)
				)

				break()
			}
		}
	}

	compact(features)
}

noisetol <- 5e-3
mass_accuracy <- 6e-6
maxcharge <- 6L


ms1_file <- readLines("15c.ms1")
ms1_file <- str_subset(ms1_file, "^(D|H|I|Z)", negate = TRUE)
split_scan <- str_detect(ms1_file, "^S") %>% cumsum()
peaks_dfs <- split(ms1_file, split_scan) %>%
	map(read_delim, " ", skip = 1, col_names = c("mz", "intensity"), col_types = "dd")

idx_unduplicated <- map(peaks_dfs, ~rep(FALSE, nrow(.x)))
x <- peaks_dfs[[1]][["mz"]]
for (t in 2:length(peaks_dfs)) {
	y <- peaks_dfs[[t]][["mz"]]
	i <- j <- 1
	I <- length(x)
	J <- length(y)

	while (i <= I && j <= J) {
		diff <- x[i] - y[j]

		if (abs(diff) <= noisetol) {
			idx_unduplicated[[t - 1]][i] <- TRUE
			idx_unduplicated[[t]][j] <- TRUE
		}

		if (diff > 0) {
			j <- j + 1
		} else {
			i <- i + 1
		}
	}
	x <- y
}

peaks_dfs <- map2(peaks_dfs, idx_unduplicated, ~filter(.x, .y))

clusters_list <- map(peaks_dfs, peaklist_to_clusters, mass_accuracy, .sort = TRUE)
features_list <- map(clusters_list, clusters_to_features, maxcharge, mass_accuracy)


df %>% mutate(cluster = mdl$cluster) %>%
	select(mz, cluster) %>%
	group_by(cluster) %>%
	nest(.key = "mz") %>%
	mutate(mz = map(mz, 1)) %>%
	mutate(minmz = map_dbl(mz, min),
		   maxmz = map_dbl(mz, max)) %>%
	select(minmz, maxmz) %>%
	write_csv("dbscan_mz_range.csv")

