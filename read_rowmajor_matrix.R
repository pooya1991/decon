read_tab_delimited <- function(f, sep = "\t", header = TRUE) {
  data <- readLines(f)
  if(header) {
    header = data[1]
    data = data[-1]
    if (is.na(header)) {
      return()
    } else {
      header <- strsplit(header, sep)[[1]]
      if(is.na(data[1])) {
        return(header)
      } else {
        result <- do.call(rbind, strsplit(data, sep))
        colnames(result) <- header
        return(result)
      }
    }
  } else {
    if(is.na(data[1])) {
      return()
    } else {
      result <- do.call(rbind, strsplit(data, sep))
      return(result)
    }
  }
}

