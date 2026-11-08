# apply rank-based inverse normal transformation to each phenotype across samples
# input and output are TSV matrices: trait ID followed by sample columns
# missing values remain NA; rows with fewer than two observed values become NA

normalize_row <- function(values) {
    observed <- which(!is.na(values))
    result <- rep(NA_real_, length(values))
    n <- length(observed)
    if (n > 1L) {
        ranks <- rank(values[observed], ties.method = "average")
        result[observed] <- qnorm((ranks - 0.5) / n)
    }
    result
}

main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) != 4L || !setequal(args[c(1L, 3L)], c("--input", "--output"))) {
        stop("usage: Rscript normalize_phenotypes.R --input phenotypes.tsv --output normalized.tsv",
             call. = FALSE)
    }
    input <- normalizePath(args[match("--input", args) + 1L], mustWork = TRUE)
    output <- args[match("--output", args) + 1L]
    data <- read.delim(input, colClasses = "character", check.names = FALSE,
                       quote = "", comment.char = "", na.strings = NULL,
                       row.names = NULL, fill = FALSE)
    header <- strsplit(readLines(input, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1L]]
    if (nrow(data) < 1L || ncol(data) < 2L || !identical(names(data), header)
            || any(!nzchar(header)) || anyDuplicated(header)
            || any(grepl("\\s", header[-1L], perl = TRUE))) {
        stop("expected a nonempty TSV matrix with distinct sample columns", call. = FALSE)
    }
    missing_tokens <- c("", "NA", "NaN", "nan", "N/A", ".")
    ids <- data[[1L]]
    if (any(ids %in% missing_tokens) || anyDuplicated(ids) || any(grepl("\\s", ids, perl = TRUE))) {
        stop("phenotype IDs must be nonempty, distinct and contain no whitespace", call. = FALSE)
    }
    raw <- as.matrix(data[-1L])
    raw[] <- trimws(raw)
    missing <- raw %in% missing_tokens
    values <- suppressWarnings(matrix(as.numeric(raw), nrow = nrow(raw), dimnames = dimnames(raw)))
    if (any(is.na(values) & !missing) || any(is.infinite(values))) {
        stop("phenotypes must contain finite numbers or missing values", call. = FALSE)
    }
    values[missing] <- NA_real_
    normalized <- t(vapply(seq_len(nrow(values)), function(i) normalize_row(values[i, ]),
                           numeric(ncol(values))))
    data[-1L] <- as.data.frame(normalized)
    dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
    output <- file.path(normalizePath(dirname(output), mustWork = TRUE), basename(output))
    if (input == normalizePath(output, mustWork = FALSE)) {
        stop("input and output must be different files", call. = FALSE)
    }
    temporary <- tempfile(".normalized-", tmpdir = dirname(output))
    on.exit(unlink(temporary), add = TRUE)
    write.table(data, temporary, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    if (!file.rename(temporary, output)) {
        stop("could not write normalized phenotype matrix", call. = FALSE)
    }
    cat(sprintf("normalized %d phenotypes across %d samples: %s\n",
                nrow(data), ncol(data) - 1L, output))
}

main()
