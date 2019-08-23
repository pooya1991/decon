get_elemental_composition <- function(features) {
    masses <- c(C = 12, H = 1.0078246, N = 14.0030732, O = 15.9949141, S = 31.972072)
    # Elemental composition of the average amino acid resideu
    averagine <- c(C = 4.9384, H = 7.7583, N = 1.3577, O = 1.4773, S = 0.0417)
    # Mass of monoisotopic averagine residue
    amass <- sum(masses * averagine)
    bind_rows(features) %>%
        mutate(meanmz = (minmz + maxmz) / 2,
               residumass = (meanmz - masses["H"]) * charge,
               aunits = residumass / amass,
               acomp = map(aunits, ~ averagine * .x)
        ) %>% pull(acomp) %>%
        map(~as.integer(round(.x))) %>%
        map(set_names, names(masses))
}

run_computems1_on_features <- function(features) {
    comp <- get_elemental_composition(features)
    comp_str <- map_chr(comp, ~paste0(names(.x), .x, collapse = "")) %>%
        paste(map_chr(features, "charge"), sep = "\t")

    ffile <- tempfile()
    writeLines(comp_str, ffile)
    sfile <- system(paste("./ComputeMS1", ffile), intern = TRUE)

    mass_mz <- str_subset(sfile, "^Average Integer") %>%
        stringr::str_split(":|,") %>%
        map(`[`, c(2, 4)) %>%
        do.call(what = rbind) %>%
        apply(2, as.numeric)

    charges <- (mass_mz[, 1] / mass_mz[, 2]) %>%
        round() %>% as.integer()

    sfile <- str_subset(sfile, "^(Sequence|\\d+)")
    split_sfile <- str_detect(sfile, "^Sequence") %>% cumsum()
    peaks_dfs <- split(sfile, split_sfile) %>%
        map(read_delim, "\t", skip = 1, col_names = c("mz", "intensity"), col_types = "dd")

    map(peaks_dfs, ~mutate(.x, intensity = intensity / sqrt(sum(intensity ^ 2))))
}
