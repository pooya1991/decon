library(dplyr)
library(stringr)
source("peakclusterforbinning.R")
source("read_rowmajor_matrix.r")

# alignment_status <- function(prec, ann, ixes, i, minlen, t) {
# 	charge <- ann["charge"]
# 	minmz <- ann["min.m.z"]
# 	maxmz <- ann["max.m.z"]
# 	if ((minmz > prec$maxmz & charge == pci$charge) | charge > prec$charge) {
# 		if (!i %in% ixes_to_continue) {
# 			flen <- prec$endscan - prec$startstan
# 			if (flan >= minlen - 1) {
# 				if (prec$to_output) {
# 					return(1L)
# 				} else {
# 					return(2L)
# 				}
# 			} else {
# 				return(3L)
# 			}
# 		} else {
# 		return(4L)
# 		}
# 	}
#
# 	if ((maxmz < prec$minmz & charge == pci$charge) | charge < prec$charge) return(5L)
#
# 	if (prec$endscan == t) {
# 		return(6L)
# 	}
#
# 	return(7L)
# }

alignment_status <- function(prec, ann, ixes, i, minlen, t) {
	charge <- ann["charge"]
	minmz <- ann["min.m.z"]
	maxmz <- ann["max.m.z"]
	if ((minmz > prec$maxmz & charge == pci$charge) | charge > prec$charge) return(1L)
	if ((maxmz < prec$minmz & charge == pci$charge) | charge < prec$charge) return(2L)
	return(3L)
}

sort_dict_charge_minmz <- function(d) {
	result <- order(d[, "t"], d[, "charge"], d[, "min.m.z"])
	# do.call(order, purrr::flatten(apply(d[, c("t", "charge", "min.m.z")]), 1, list))
	return(result)
}

output_precursor <- function(fout, pci) {
	res <- paste(pci$minmz, '\t', pci$maxmz, '\t',  pci$charge, '\t', pci$startscan, '\t', pci$endscan, '\t', sep = "")
	indices <- pci$indices
	indices <- paste(indices, collapse = ",")
	res <- paste0(res, indices)
	write(res, fout, append = TRUE)
}

d <- "15b_ms1_sparsebinned_singlescans_tannotated/"
d <- paste(getwd(), d, sep = "/")

massHMono <- 1.0078246

allfiles_annfile <- list.files(d, pattern = "annfile")
allfiles_swp <- list.files(d, pattern = ".swp")
allfiles_joint <- list.files(d, pattern = "joint")

annotationfiles <- allfiles_annfile[!allfiles_annfile %in% allfiles_swp & !allfiles_annfile %in% allfiles_joint]
annotationfiles <- annotationfiles %>%
	as.data.frame() %>%
	mutate(number = as.numeric(str_extract(., "[0-9]+")),
		   text = str_extract(., "[aA-zZ]+")) %>%
	arrange(number) %>%
	pull(1) %>% as.character()

cat(annotationfiles, "\n")

af <- annotationfiles[1]
af <- paste(d, af, sep = "/")
anns <- read.table(af, header = TRUE, sep = "\t")

tanns <- segregate_annotations_by_time(anns)
tanns <- purrr::map(tanns, as.matrix)
ts = names(tanns) %>% as.numeric() %>% sort()
t = ts[1]
anns = tanns[[t]]
N = nrow(anns)
order = sort_dict_charge_minmz(anns)
precursors <- vector("list", N)
for (i in 1:N) {
	precursors[[i]] = precursor$new(anns[i,"min.m.z"], anns[i,"max.m.z"], anns[i, "charge"], t, order[i])
}

numfinished = 0
lens <- vector()
minlen = 4

tindex = 2
afindex = 1

outfile = paste0(d, 'joint_annfile_minlen', minlen, '.txt')

fout <- file(outfile, "w")

write("start m/z\tend m/z\tcharge\tstart file\tend file\tindices", fout, append = TRUE)
total = length(precursors)
unaligned = 0
ttol = 1

while (tindex <= length(ts)) {
	t <- ts[tindex]
	prevt <- as.numeric(t) - ttol
	anns <- tanns[[t]]
	total <- total + length(anns[,"charge"])
	tindex <- tindex + 1

	if (tindex %% 50 == 0) {
		cat('Num features evaluated:', numfinished, "\n")
		cat('Features that do not persist over ', minlen, ' scans\n')
		cat(unaligned, '/', total, ',', unaligned/total)
	}

	if (tindex == length(ts) & afindex < length(annotationfiles)) {
		af <- annotationfiles[afindex]
		af <- paste0(d, af)
		afindex <- afindex + 1
		tanns <- segregate_annotations_by_time(read_tab_delimited(af))
		ts <- names(tanns) %>% as.numeric() %>% sort() %>% as.character()
	}

	precursors <- precursors[order(purrr::map_dbl(precursors, "charge"),
								   purrr::map_dbl(precursors, "minmz"))]
	order = sort_dict_charge_minmz(anns)
	anns <- anns[order, ]
	J <- nrow(anns)
	I <- length(precursors)
	ixes_to_continue <- integer()

	i <- j <- 1

	while (i <= I & j <= J) {
		curr_ann <- anns[j, ]
		pci <- precursors[[i]]

		alignment_status_code <- alignment_status(pci, curr_ann)

		# They don't align
		if (alignment_status_code == 1L) {
			if (!i %in% ixes_to_continue) {
				flen <- pci$endscan - pci$startscan
				if (flen >= minlen - 1) {
					if (pci$to_output) {
						output_precursor(fout, pci)
						k <- length(lens) + 1
						lens[k] <- flen
					}
					numfinished = numfinished + 1
				} else {
					unaligned = unaligned + length(pci$indices)
				}
			}
			i = i + 1
			next()
		}

		# They don't align
		if (alignment_status_code == 2L) {
			ixes_to_continue <- c(ixes_to_continue, I + 1)
			precursors[[I + 1]] <- precursor$new(anns[j,"min.m.z"], anns[j,"max.m.z"], anns[j, "charge"], t, order[j])
			j = j + 1
			next()
		}

		# They do align! add it to an old one
		if (pci$endscan == t){
			oldix <- pci$indices[-1]
			pci$indices <- c(pci$indices, -order[j])
		} else {
			pci$indices <- c(pci$indices, order[j])
		}
		pci$endscan <- t
		pci$minmz <- min(pci$minmz, curr_ann['min.m.z'])
		pci$maxmz <- min(pci$maxmz, curr_ann['max.m.z'])
		ixes_to_continue <- c(ixes_to_continue, i)

		# See if anns[*][j] connects to subsequent Is. If the subsequent ones do, combine them with this
		# original i.
		k <- i + 1
		while (k <= I) {
			pck <- precursors[[k]]
			if ((pck$minmz <= curr_ann['max.m.z']) && (pck$charge == curr_ann['charge'])) {
				# Combine pci and precursors[k]!
				#print i, k
				# browser()
				pci$incorporate( pck )
			}

			k = k + 1
		}
		j = j + 1
	}

	while (j < J) {
		ixes_to_continue <- c(ixes_to_continue, length(precursors))
		dummy_counter <- length(precursors) + 1
		precursors[[dummy_counter]] <- precursor$new(anns[j,"min.m.z"], anns[j,"max.m.z"], anns[j, "charge"], t, order[j])
		j = j + 1
	}
	precursors = precursors[ixes_to_continue]
}

cat("Unaligned")

if (TRUE) {
	numscansperblock = 5
	prec <- data.frame(charge = unlist(lapply(precursors, "[[", "charge")),
					   minmz = unlist(lapply(precursors, "[[", "minmz")))
	prec_order <- order(prec$charge, prec$minmz)
	precursors <- precursors[prec_order]
	prec_len <- length(precursors)
	for (i in 1:prec_len) {
		pci <- precursors[[i]]
		flen = as.numeric(pci$endscan) - as.numeric(pci$startscan)
		if (flen >= minlen - 1){
			output_precursor(fout, pci)
			lens <- c(lens, flen)
			numfinished = numfinished + 1
		}
	}
}

close(fout)