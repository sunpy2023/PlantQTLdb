# estimate PEER factors from a filtered FPKM TSV matrix: gene ID followed by samples
# apply log2(FPKM + 1), filter genes by missingness and check sample missingness
# fill remaining missing values with gene means only in the PEER input matrix
# all sample IDs and their order are retained; phenotype input files are unchanged

main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    usage <- paste("usage: Rscript estimate_peer_factors.R --input expression.filtered.tsv",
                   "--output peer_factors.tsv [--factor-count 20]",
                   "[--max-gene-missing 0.05] [--max-sample-missing 0.05]")
    if (identical(args, "--help")) {
        cat(usage, "\n")
        return(invisible(NULL))
    }
    allowed <- c("--input", "--output", "--factor-count", "--max-gene-missing", "--max-sample-missing")
    if (length(args) < 4L || length(args) %% 2L) stop(usage, call. = FALSE)
    flags <- args[seq.int(1L, length(args), 2L)]
    if (anyDuplicated(flags) || any(!flags %in% allowed)
            || !all(c("--input", "--output") %in% flags)) stop(usage, call. = FALSE)
    options <- c("--factor-count" = "20", "--max-gene-missing" = "0.05", "--max-sample-missing" = "0.05")
    options[flags] <- args[seq.int(2L, length(args), 2L)]
    k <- suppressWarnings(as.numeric(options[["--factor-count"]]))
    limits <- suppressWarnings(as.numeric(options[c("--max-gene-missing", "--max-sample-missing")]))
    if (!is.finite(k) || k < 1 || k != floor(k) || any(!is.finite(limits) | limits < 0 | limits > 1)) {
        stop("factor count must be a positive integer; missingness limits must be within [0, 1]",
             call. = FALSE)
    }
    input <- normalizePath(options[["--input"]], mustWork = TRUE)
    output <- options[["--output"]]
    dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
    output <- file.path(normalizePath(dirname(output), mustWork = TRUE), basename(output))
    if (input == normalizePath(output, mustWork = FALSE)) {
        stop("input and output must be different files", call. = FALSE)
    }
    data <- read.delim(input, colClasses = "character", check.names = FALSE,
                       quote = "", comment.char = "", na.strings = NULL,
                       row.names = NULL, fill = FALSE)
    header <- strsplit(readLines(input, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1L]]
    if (nrow(data) < 1L || ncol(data) < 3L || !identical(names(data), header)
            || any(!nzchar(header)) || anyDuplicated(header)
            || any(grepl("\\s", header[-1L], perl = TRUE))) {
        stop("expected a nonempty FPKM TSV with distinct sample columns", call. = FALSE)
    }
    missing_tokens <- c("", "NA", "NaN", "nan", "N/A", ".")
    ids <- data[[1L]]
    if (any(ids %in% missing_tokens) || anyDuplicated(ids) || any(grepl("\\s", ids, perl = TRUE))) {
        stop("gene IDs must be nonempty, distinct and contain no whitespace", call. = FALSE)
    }
    raw <- as.matrix(data[-1L])
    raw[] <- trimws(raw)
    missing <- raw %in% missing_tokens
    expression <- suppressWarnings(matrix(as.numeric(raw), nrow = nrow(raw), dimnames = dimnames(raw)))
    if (any(is.na(expression) & !missing) || any(is.infinite(expression))
            || any(expression < 0, na.rm = TRUE)) {
        stop("FPKM must contain finite nonnegative numbers or missing values", call. = FALSE)
    }
    expression[missing] <- NA_real_
    expression <- expression[rowMeans(is.na(expression)) <= limits[1L]
                             & rowSums(!is.na(expression)) >= 2L, , drop = FALSE]
    if (k >= min(dim(expression))) {
        stop("factor count must be smaller than the sample count and retained gene count", call. = FALSE)
    }
    failed <- colMeans(is.na(expression)) > limits[2L]
    if (any(failed)) {
        stop("samples exceed the PEER missingness limit: ",
             paste(colnames(expression)[failed], collapse = ", "),
             "; no samples were removed", call. = FALSE)
    }
    expression <- log2(expression + 1)
    missing <- which(is.na(expression), arr.ind = TRUE)
    if (nrow(missing)) expression[missing] <- rowMeans(expression, na.rm = TRUE)[missing[, 1L]]
    suppressPackageStartupMessages(library(peer))
    model <- PEER()
    PEER_setPhenoMean(model, t(expression))
    PEER_setNk(model, as.integer(k))
    PEER_update(model)
    factors <- PEER_getX(model)
    if (!identical(dim(factors), c(ncol(expression), as.integer(k))) || any(!is.finite(factors))) {
        stop("PEER returned invalid factor dimensions or non-finite values", call. = FALSE)
    }
    result <- data.frame(id = paste0("PEER", seq_len(k)), t(factors), check.names = FALSE)
    names(result) <- c("id", colnames(expression))
    temporary <- tempfile(".peer-", tmpdir = dirname(output))
    on.exit(unlink(temporary), add = TRUE)
    write.table(result, temporary, sep = "\t", quote = FALSE, row.names = FALSE)
    if (!file.rename(temporary, output)) stop("could not write PEER factors", call. = FALSE)
    cat(sprintf("estimated %d PEER factors from %d/%d genes and all %d samples: %s\n",
                k, nrow(expression), nrow(data), ncol(expression), output))
    cat(sprintf("filled %d missing log2(FPKM + 1) values for PEER only\n", nrow(missing)))
}

main()
