get_elemental_composition <- function(avgmz, charge, masses, amass, averagine) {
    residuemass <- (avgmz - masses["H"]) * charge
    aunits <- residuemass / amass
    acomp <- as.integer(round(averagine * aunits))
    names(acomp) <- names(masses)
    averagine_mass <- sum(acomp * masses)
    numH <- as.integer((residuemass - averagine_mass) / masses["H"])
    acomp["H"] <- acomp["H"] + numH
    averagine_mass <- averagine_mass + numH * masses["H"]
    return(list(acomp, averagine_mass))
}

normalize_peaks <- function(peaks) {
    scalars <- purrr::map(peaks, ~purrr::map_dbl(.x, 2)) %>%
        purrr::map_dbl(function(x) 1 / sqrt(sum(x ^ 2)))
    purrr::map2(peaks, scalars, ~purrr::map(.x, function(x) x * c(1, .y)))
}

run_computems1_on_features <- function(featureinfo, d, ffile, sfile, ...) {
    to_print <- TRUE
    ffile <- paste0(d, ffile)
    sfile <- paste0(d, sfile)
    submitted <- vector("list", length(featureinfo))
    formulanum <- -1
    duplicates <- 0
    formulas <- character()
    compstring <- character(length(featureinfo))
    # fout <- file(ffile, open = "w", raw = TRUE)

    i <- 1
    for (vec in featureinfo) {
        n <- vec[1]; charge <- vec[2]; mz1 <- vec[3]; mz2 <- vec[4]
        smz <- (mz1 + mz2) / 2
        comp_mass <- get_elemental_composition(smz, charge, ...)
        comp <- comp_mass[[1]]
        mass <- comp_mass[[2]]
        compstring[i] <- paste0(names(comp), comp, collapse = "") %>%
            paste0("\t", charge)

        if (!compstring[i] %in% formulas) {
            formulanum <- formulanum + 1
        } else {
            duplicates <- duplicates + 1
        }

        submitted[[i]] <- c(smz, charge, formulanum)
        i <- i + 1
    }

    write(compstring, file = ffile, sep = "\n", append = FALSE)
    formulanum <- formulanum + 1

    if (toprint) {
        cat("Submitted", formulanum, "\n")
        cat("Duplicates", duplicates, "\n")
        cat("Running computems1 command to", sfile, "\n")
    }

    sfile_entry <- system(paste("./ComputeMS1.exe", ffile), intern = TRUE)
    write(sfile_entry, sfile, append = FALSE)
    # system(paste0("cmd /c echo|set /p=", sfile_entry, ">>", sfile))

    fin <- readLines(sfile)
    charges <- integer()
    peaks <- vector("list", formulanum)
    i <- 0
    isvals <- FALSE

    for (line in fin) {
        if (str_detect(line, "successful")) {
            isvals <- FALSE
            next()
        }

        if (str_detect(line, "Sequence")) {
            i <- i + 1
        }

        if (str_detect(line, "Average Integer")) {
            mass <- as.numeric(str_extract(line, "(?<=MW:\\s).*(?=,)"))
            mz <- as.numeric(str_extract(line, "(?<=m/z:\\s).*"))
            charges <- c(charges, as.integer(round(mass / mz)))
        }

        if (str_detect(line, "Calculation")) {
            isvals <- TRUE
            next()
        }

        if (isvals) {
            tokens <- str_split(line, "\t")[[1]]
            mz <- as.numeric(tokens[1])
            intensity <- as.numeric(tokens[2])
            if (intensity < 0.1) {
                next()
            }
            peaks[[i]] <- c(peaks[[i]], list(c(mz, intensity)))
        }
    }
    peaks <- normalize_peaks(peaks)
    if (toprint) {
        cat("Specfile:", sfile, "lenpeaks:", length(peaks), "\n")
    }
    return(peaks)
}
