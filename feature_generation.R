tic <- proc.time()
ms1_file <- readLines("15c.ms1")
ms1_file <- str_subset(ms1_file, "^(D|H|I|Z)", negate = TRUE)
split_scan <- str_detect(ms1_file, "^S") %>% cumsum()
peaks_df <- split(ms1_file, split_scan) %>%
	map(read_delim, " ", skip = 1, col_names = c("mz", "intensity"), col_types = "dd")

# peaklist <- map(b15_splitted,
# 				~transmute(.x, mz_intensity = map2(mz, intensity, c)) %>% pull(mz_intensity))

tol <- 5e-3
idx_unduplicated <- map(peaks_df, ~rep(FALSE, nrow(.x)))
x <- peaks_df[[1]][["mz"]]
for (t in 2:length(idx_unduplicated)) {
	y <- peaks_df[[t]][["mz"]]
	i <- j <- 1
	I <- length(x)
	J <- length(y)

	while (i <= I && j <= J) {
		diff <- x[i] - y[j]
		if (abs(diff) <= tol) {
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

peaks_df <- map2(peaks_df, idx_unduplicated, ~filter(.x, .y))
toc <- proc.time()

df <- bind_rows(peaks_df, .id = "scan")

peakslist <- map(peaks_df,
				~transmute(.x, mz_intensity = map2(mz, intensity, c)) %>% pull(mz_intensity))

peaklist <- peakslist[[1]]
tol = 0.01
mz <- peaklist[[1]][1]
intensity <- peaklist[[1]][2]
clusters <- list(list(minmz = mz, maxmz = mz, peaks = intensity))

for (vec in peaklist[-1]) {
	mz <- vec[1]; intensity <- vec[2]
	clust_curr <- clusters[[length(clusters)]]
	if (clust_curr$maxmz > (mz * (1 - 6e-6))) {
		clust_curr$maxmz <- mz
		clust_curr$peaks <- c(clust_curr$peaks, intensity)
		clusters[[length(clusters)]] <- clust_curr
	} else {
		clusters <- c(clusters, list(list(minmz = mz, maxmz = mz, peaks = intensity)))
	}
}

len_clusters <- length(clusters)
features <- list()
maxcharge <- 4L

for (n in seq_len(len_clusters)) {
	clust_curr <- clusters[[n]]
	n2 <- n + 1L

	for (charge in maxcharge:1L) {
		isop <- c(minmz = clust_curr$minmz + 1 / charge, maxmz = clust_curr$maxmz + 1 / charge)
		while (n2 <= len_clusters) {
			if (clusters[[n2]]$maxmz + tol < isop["minmz"]) {
				n2 <- n2 + 1L
				next()
			}

			if (isop["maxmz"] + tol < clusters[[n2]]$minmz) break()
			# There is a match. Make this feature and add it to a list
			features <- c(features, list(c(n, charge, clust_curr$minmz - tol, clust_curr$maxmz + tol)))
			break()
		}
	}
}


df %>% mutate(cluster = mdl$cluster) %>%
	select(mz, cluster) %>%
	group_by(cluster) %>%
	nest(.key = "mz") %>%
	mutate(mz = map(mz, 1)) %>%
	mutate(minmz = map_dbl(mz, min),
		   maxmz = map_dbl(mz, max)) %>%
	select(minmz, maxmz) %>%
	write_csv("dbscan_mz_range.csv")

