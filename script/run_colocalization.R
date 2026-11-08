# run coloc.abf on molecular features selected by PLINK GWAS lead variants
# inputs are written by prepare_coloc_inputs.py; both phenotypes are quantitative
# select features with a cis-QTL variant in a lead window, then use all shared cis variants
# MAF and N estimate sdY, as in the original workflow; do not assume sdY = 1
# outputs: <prefix>.all.tsv, .suggestive.tsv, .strong.tsv and .status.tsv

suppressPackageStartupMessages({
    library(coloc)
    library(data.table)
})

read_statistics <- function(path, qtl = FALSE) {
    columns <- c("SNP", "chr", "position", "beta", "varbeta", "MAF", "p", "n_samples")
    if (qtl) columns <- c("feature", columns)
    types <- c(SNP = "character", chr = "character", if (qtl) c(feature = "character"))
    table <- fread(path, colClasses = types, na.strings = "NA")
    if (anyDuplicated(names(table)) || !all(columns %in% names(table))) {
        stop("missing or duplicate statistics columns: ", path)
    }
    ids <- c("SNP", "chr", if (qtl) "feature")
    for (column in ids) {
        values <- as.character(table[[column]])
        if (anyNA(values) || any(!nzchar(values) | grepl("\\s", values) |
                                values %in% c(".", "NaN", "nan", "N/A", "None"))) {
            stop("missing or invalid ", column, " IDs: ", path)
        }
        set(table, j = column, value = values)
    }
    for (column in setdiff(columns, ids)) {
        values <- suppressWarnings(as.numeric(table[[column]]))
        if (any(!is.finite(values))) stop("nonfinite or missing ", column, ": ", path)
        set(table, j = column, value = values)
    }
    if (any(!grepl("^[1-9][0-9]*$", table$chr)) ||
            any(table$position <= 0 | table$position != floor(table$position)) ||
            any(table$varbeta <= 0 | table$MAF <= 0 | table$MAF > 0.5) ||
            any(table$p < 0 | table$p > 1) ||
            any(table$n_samples <= 2 | table$n_samples != floor(table$n_samples))) {
        stop("invalid coordinates, variances, MAF, P-values or sample counts: ", path)
    }
    variants <- unique(table[, .(SNP, chr, position)])
    if (anyDuplicated(variants$SNP)) {
        stop("each variant ID must identify one chromosome and position: ", path)
    }
    keys <- c(if (qtl) "feature", "SNP")
    if (anyDuplicated(table[, ..keys])) stop("duplicate association records: ", path)
    if (qtl) {
        metadata <- table[, .(n = uniqueN(n_samples), chromosomes = uniqueN(chr)), by = feature]
        if (any(metadata$n != 1L | metadata$chromosomes != 1L)) {
            stop("each molecular feature requires one sample count and chromosome")
        }
    } else if (nrow(table) && uniqueN(table$n_samples) != 1L) {
        stop("GWAS sample count must match the observed trait count")
    }
    table
}

lead_variants <- function(prefix, gwas) {
    log <- readLines(paste0(prefix, ".log"), warn = FALSE)
    if (!any(grepl("^End time:", log)) || any(grepl("^Error:", log))) {
        stop("PLINK clumping did not complete successfully; inspect ", prefix, ".log")
    }
    if (any(grepl("No significant --clump results", log, fixed = TRUE))) return(character())
    lead <- fread(paste0(prefix, ".clumped"), colClasses = c(SNP = "character"))
    if (!"SNP" %in% names(lead) || !nrow(lead) || anyNA(lead$SNP) || anyDuplicated(lead$SNP) ||
            !all(lead$SNP %in% gwas$SNP)) stop("invalid or unmatched PLINK lead variants")
    lead$SNP
}

candidate_features <- function(qtl, leads, window) {
    windows <- unique(leads[, .(chr, start = pmax(1, position - window), end = position + window)])
    setorder(windows, chr, start, end)
    windows[, previous_end := shift(cummax(end), fill = -Inf), by = chr]
    windows[, group := cumsum(start > previous_end), by = chr]
    windows <- windows[, .(start = min(start), end = max(end)), by = .(chr, group)]
    unique(qtl[windows, on = .(chr, position >= start, position <= end),
               .(feature = x.feature), nomatch = 0L]$feature)
}

empty_summary <- function() {
    data.table(feature = character(), n_shared = integer(), PP.H0 = numeric(), PP.H1 = numeric(),
               PP.H2 = numeric(), PP.H3 = numeric(), PP.H4 = numeric(), top_variant = character(),
               top_variant_PP.H4 = numeric(), top_variant_p_qtl = numeric(),
               top_variant_p_gwas = numeric(), n_qtl = integer(), n_gwas = integer())
}

run_analysis <- function(qtl, gwas, features, options) {
    setkey(qtl, feature)
    rows <- vector("list", length(features))
    status <- data.table(feature = features, status = rep("ok", length(features)),
                         n_shared = rep(NA_integer_, length(features)), detail = rep("", length(features)))
    for (index in seq_along(features)) {
        feature_id <- features[index]
        pairs <- merge(qtl[J(feature_id)], gwas, by = c("SNP", "chr", "position"),
                       suffixes = c("_qtl", "_gwas"))
        status[index, n_shared := nrow(pairs)]
        if (nrow(pairs) < options$min_shared) {
            status[index, status := "insufficient_shared_variants"]
            next
        }
        dataset <- function(suffix) {
            list(beta = pairs[[paste0("beta_", suffix)]],
                 varbeta = pairs[[paste0("varbeta_", suffix)]],
                 MAF = pairs[[paste0("MAF_", suffix)]],
                 N = unique(pairs[[paste0("n_samples_", suffix)]]), type = "quant",
                 snp = pairs$SNP, position = pairs$position)
        }
        result <- tryCatch(coloc.abf(dataset("gwas"), dataset("qtl"),
                                    p1 = options$p1, p2 = options$p2, p12 = options$p12),
                           error = function(error) error)
        if (inherits(result, "error")) {
            status[index, `:=`(status = "coloc_error", detail = conditionMessage(result))]
            next
        }
        top <- result$results[which.max(result$results$SNP.PP.H4), ]
        association <- pairs[SNP == top$snp]
        rows[[index]] <- data.table(
            feature = feature_id, n_shared = nrow(pairs),
            PP.H0 = unname(result$summary["PP.H0.abf"]),
            PP.H1 = unname(result$summary["PP.H1.abf"]),
            PP.H2 = unname(result$summary["PP.H2.abf"]),
            PP.H3 = unname(result$summary["PP.H3.abf"]),
            PP.H4 = unname(result$summary["PP.H4.abf"]),
            top_variant = top$snp, top_variant_PP.H4 = top$SNP.PP.H4,
            top_variant_p_qtl = association$p_qtl, top_variant_p_gwas = association$p_gwas,
            n_qtl = pairs$n_samples_qtl[1L], n_gwas = pairs$n_samples_gwas[1L])
    }
    summary <- rbindlist(rows)
    if (!nrow(summary)) summary <- empty_summary()
    setorder(summary, -PP.H4, -n_shared, feature)
    list(summary = summary, status = status)
}

main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    required <- c("--qtl", "--gwas", "--clump-prefix", "--out-prefix")
    defaults <- c("--window-bp" = "200000", "--min-shared" = "50", "--p1" = "1e-4",
                  "--p2" = "1e-4", "--p12" = "1e-5", "--pph4-min" = "0.5", "--pph4-strong" = "0.8")
    usage <- paste("usage: Rscript script/run_colocalization.R --qtl qtl.coloc.tsv.gz",
                   "--gwas gwas.coloc.tsv.gz --clump-prefix gwas --out-prefix results/trait")
    if (identical(args, "--help")) { cat(usage, "\n"); return(invisible(NULL)) }
    if (!length(args) || length(args) %% 2L) stop(usage, call. = FALSE)
    flags <- args[seq.int(1L, length(args), 2L)]
    if (anyDuplicated(flags) || !all(required %in% flags) ||
            !all(flags %in% c(required, names(defaults)))) stop(usage, call. = FALSE)
    options <- defaults
    options[flags] <- args[seq.int(2L, length(args), 2L)]
    values <- suppressWarnings(as.numeric(options[names(defaults)]))
    names(values) <- sub("^--", "", names(defaults))
    if (any(!is.finite(values)) || values[1L] < 0 || values[1L] != floor(values[1L]) ||
            values[2L] < 1 || values[2L] != floor(values[2L]) ||
            any(values[3:7] <= 0 | values[3:7] > 1) || values[6L] > values[7L] ||
            values[5L] > min(values[3:4]) || values[3L] + values[4L] - values[5L] >= 1) {
        stop("invalid window, shared-variant count, prior probabilities or PP.H4 thresholds")
    }
    inputs <- c(options[["--qtl"]], options[["--gwas"]], paste0(options[["--clump-prefix"]], ".log"))
    clumped <- paste0(options[["--clump-prefix"]], ".clumped")
    if (file.exists(clumped)) inputs <- c(inputs, clumped)
    inputs <- normalizePath(inputs, mustWork = TRUE)
    prefix <- options[["--out-prefix"]]
    outputs <- paste0(prefix, c(".all.tsv", ".suggestive.tsv", ".strong.tsv", ".status.tsv"))
    if (any(normalizePath(outputs, mustWork = FALSE) %in% inputs)) {
        stop("output files must not overwrite inputs")
    }
    qtl <- read_statistics(options[["--qtl"]], TRUE)
    gwas <- read_statistics(options[["--gwas"]])
    if (!nrow(gwas)) stop("GWAS input has no valid records")
    shared <- merge(unique(qtl[, .(SNP, chr, position)]), gwas[, .(SNP, chr, position)], by = "SNP")
    if (any(shared$chr.x != shared$chr.y | shared$position.x != shared$position.y)) {
        stop("QTL and GWAS coordinates disagree for shared variant IDs")
    }
    rm(shared)
    lead_ids <- lead_variants(options[["--clump-prefix"]], gwas)
    features <- if (length(lead_ids) && nrow(qtl)) {
        candidate_features(qtl, gwas[SNP %in% lead_ids], values[["window-bp"]])
    } else character()
    parameters <- list(min_shared = values[["min-shared"]], p1 = values[["p1"]],
                       p2 = values[["p2"]], p12 = values[["p12"]])
    result <- run_analysis(qtl, gwas, features, parameters)
    if (!length(features)) {
        reason <- if (!length(lead_ids)) "no_gwas_leads" else "no_features_in_lead_windows"
        result$status <- data.table(feature = NA_character_, status = reason, n_shared = 0L, detail = "")
    }
    tables <- list(result$summary, result$summary[PP.H4 >= values[["pph4-min"]]],
                   result$summary[PP.H4 >= values[["pph4-strong"]]], result$status)
    dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
    for (index in seq_along(outputs)) {
        temporary <- tempfile(".coloc-", tmpdir = dirname(prefix))
        fwrite(tables[[index]], temporary, sep = "\t", na = "NA")
        if (!file.rename(temporary, outputs[index])) stop("cannot save ", outputs[index])
    }
    message("coloc ", packageVersion("coloc"), ": ", length(lead_ids), " GWAS leads; ",
            length(features), " candidate features; ", nrow(result$summary), " successful analyses; ",
            nrow(tables[[2L]]), " suggestive; ", nrow(tables[[3L]]), " strong")
    if (any(result$status$status == "coloc_error")) stop("coloc errors recorded in ", outputs[4L])
}

if (sys.nframe() == 0L) main()
