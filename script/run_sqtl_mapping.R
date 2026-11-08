# map cis/trans-sQTLs from a QC VCF, normalized junction PSI and complete covariates
# shared observed-sample regression also supports run_apaqtl_mapping.R
# run from the pipeline root; reuse table/genotype helpers without running eQTL mapping
# raw PSI QC and rank normalization are performed in phenotype preparation
# PSI NA are preserved; each junction is fitted using its observed samples only
# events with the same observed samples are fitted together; covariates must remain full rank
# genotype QC precedes grouping; remaining GT NA use each group's observed variant mean
# beta defaults to per REF allele; --counted-allele ALT retains the VCF GT dosage direction
# cis/trans significance threshold: 0.05 / retained variants, including SNPs and indels
# outputs: cis.significant.tsv, trans.significant.tsv, cis_splicing_traits.normalized.tsv,
# cis_splicing_traits.all_cis.tsv; association tables include n_samples and df

source("script/run_eqtl_mapping.R")

observed_groups <- function(values, covariates) {
    keys <- apply(!is.na(values), 1L, function(values) paste(which(values), collapse = ","))
    groups <- unname(split(seq_len(nrow(values)), keys))
    for (rows in groups) {
        observed <- !is.na(values[rows[1L], ])
        tryCatch(check_design(covariates[, observed, drop = FALSE]),
                 error = function(e) stop("cannot fit phenotype ", rownames(values)[rows[1L]], ": ",
                                          conditionMessage(e), call. = FALSE))
    }
    groups
}

empty_associations <- function(id_column = "junction_id") {
    table <- data.frame(SNP = character(), junction_id = character(), beta = numeric(),
                        "t-stat" = numeric(), "p-value" = numeric(), FDR = numeric(),
                        n_samples = integer(), df = integer(), check.names = FALSE)
    names(table)[2L] <- id_column
    table
}

append_cis <- function(path, output, n_samples, df, id_column = "junction_id") {
    connection <- file(path, "r")
    on.exit(close(connection))
    header <- strsplit(readLines(connection, n = 1L), "\t", fixed = TRUE)[[1L]]
    if (!identical(header, c("SNP", "gene", "beta", "t-stat", "p-value"))) {
        stop("unexpected MatrixEQTL all-cis columns")
    }
    repeat {
        lines <- readLines(connection, n = 10000L, warn = FALSE)
        if (!length(lines)) break
        table <- fread(text = paste0(paste(lines, collapse = "\n"), "\n"), sep = "\t",
                       header = FALSE, col.names = header, colClasses = list(character = 1:2),
                       data.table = FALSE, showProgress = FALSE)
        names(table)[2L] <- id_column
        table$n_samples <- n_samples
        table$df <- df
        fwrite(table, output, sep = "\t", quote = FALSE, append = TRUE)
    }
}

run_observed_models <- function(snps, positions, values, covariates, locations, distance, threshold,
                                all_cis_file = NULL, id_column = "junction_id") {
    groups <- observed_groups(values, covariates)
    parts <- list(cis = vector("list", length(groups)), trans = vector("list", length(groups)))
    tests <- c(cis = 0, trans = 0)
    if (!is.null(all_cis_file)) {
        table <- empty_associations(id_column)
        fwrite(table[, names(table) != "FDR"], all_cis_file, sep = "\t", quote = FALSE)
    }
    for (i in seq_along(groups)) {
        rows <- groups[[i]]
        observed <- !is.na(values[rows[1L], ])
        n_samples <- sum(observed)
        df <- n_samples - nrow(covariates) - 2L
        group_snps <- snps$Clone()
        group_snps$ColumnSubsample(which(observed))
        export <- if (is.null(all_cis_file)) "" else paste0(all_cis_file, ".group")
        result <- Matrix_eQTL_main(
            snps = group_snps, gene = sliced_matrix(values[rows, observed, drop = FALSE]),
            cvrt = sliced_matrix(covariates[, observed, drop = FALSE]), useModel = modelLINEAR,
            output_file_name = "", pvOutputThreshold = if (is.null(all_cis_file)) threshold else 0,
            output_file_name.cis = export,
            pvOutputThreshold.cis = if (is.null(all_cis_file)) threshold else 1,
            snpspos = positions, genepos = locations[rows, , drop = FALSE], cisDist = distance,
            pvalue.hist = FALSE, min.pv.by.genesnp = FALSE,
            noFDRsaveMemory = !is.null(all_cis_file), verbose = FALSE)
        if (is.null(all_cis_file)) {
            for (kind in names(parts)) {
                table <- as.data.frame(result[[kind]]$eqtls)
                tests[kind] <- tests[kind] + result[[kind]]$ntests
                if (!nrow(table)) {
                    table <- empty_associations(id_column)
                } else {
                    names(table)[match(c("snps", "gene", "statistic", "pvalue"), names(table))] <-
                        c("SNP", id_column, "t-stat", "p-value")
                    table$n_samples <- n_samples
                    table$df <- df
                    table <- table[, names(empty_associations(id_column))]
                }
                parts[[kind]][[i]] <- table
            }
        } else {
            append_cis(export, all_cis_file, n_samples, df, id_column)
            unlink(export)
        }
        cat(sprintf("observed-sample group %d/%d: %d phenotypes; %d samples; df = %d\n",
                    i, length(groups), length(rows), n_samples, df))
    }
    if (!is.null(all_cis_file)) return(invisible(NULL))
    for (kind in names(parts)) {
        table <- as.data.frame(rbindlist(parts[[kind]]))
        if (nrow(table)) {
            table <- table[order(table[["p-value"]]), , drop = FALSE]
            # use the total tests across groups, rather than each group's FDR
            table$FDR <- p.adjust(table[["p-value"]], method = "BH", n = tests[kind])
        }
        parts[[kind]] <- table
    }
    parts
}

main <- function(phenotype = "splicing") {
    settings <- switch(phenotype,
        splicing = c(value = "--psi", location = "--junction-locations", label = "PSI",
                     trait = "junctions", id = "junction_id", task = "sqtl"),
        apa = c(value = "--pdui", location = "--apa-locations", label = "PDUI",
                trait = "APA traits", id = "apa_id", task = "apaqtl"))
    if (is.null(settings)) stop("unsupported phenotype type")
    value_flag <- settings[["value"]]
    location_flag <- settings[["location"]]
    id_column <- settings[["id"]]
    label <- settings[["label"]]
    args <- commandArgs(trailingOnly = TRUE)
    required <- c("--vcf", value_flag, location_flag, "--covariates", "--out-dir")
    defaults <- c("--cis-distance" = "200000", "--min-maf" = "0.05", "--max-variant-missing" = "0.10",
                  "--counted-allele" = "REF")
    usage <- paste("usage: Rscript script/run_", settings[["task"]], "_mapping.R --vcf genotype.qc.vcf.gz ",
                   value_flag, " ", phenotype, ".normalized.tsv ", location_flag, " ", phenotype,
                   ".locations.tsv --covariates covariates.tsv --out-dir ", settings[["task"]],
                   " [--cis-distance 200000] [--min-maf 0.05] [--max-variant-missing 0.10]",
                   " [--counted-allele REF|ALT]", sep = "")
    if (identical(args, "--help")) { cat(usage, "\n"); return(invisible(NULL)) }
    if (length(args) < 10L || length(args) %% 2L) stop(usage, call. = FALSE)
    flags <- args[seq.int(1L, length(args), 2L)]
    if (anyDuplicated(flags) || !all(required %in% flags)
            || any(!flags %in% c(required, names(defaults)))) stop(usage, call. = FALSE)
    options <- defaults
    options[flags] <- args[seq.int(2L, length(args), 2L)]
    parameters <- suppressWarnings(as.numeric(options[names(defaults)[1:3]]))
    if (any(!is.finite(parameters)) || parameters[1L] < 0 || parameters[1L] != floor(parameters[1L])
            || parameters[2L] < 0 || parameters[2L] > 0.5 || parameters[3L] < 0 || parameters[3L] > 1) {
        stop("invalid cis distance, MAF or missingness threshold")
    }
    counted_allele <- options[["--counted-allele"]]
    if (!counted_allele %in% c("REF", "ALT")) stop("counted allele must be REF or ALT")
    if (!nzchar(Sys.which("bcftools"))) stop("bcftools is required on PATH")
    inputs <- vapply(options[required[1:4]], normalizePath, character(1), mustWork = TRUE)
    out_dir <- options[["--out-dir"]]
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_dir <- normalizePath(out_dir, mustWork = TRUE)
    output_names <- c("cis.significant.tsv", "trans.significant.tsv",
                      paste0("cis_", phenotype, "_traits.normalized.tsv"),
                      paste0("cis_", phenotype, "_traits.all_cis.tsv"))
    outputs <- file.path(out_dir, output_names)
    if (any(vapply(outputs, normalizePath, character(1), mustWork = FALSE) %in% inputs)) {
        stop("output files must not overwrite inputs")
    }
    values <- read_matrix(options[[value_flag]], allow_missing = TRUE)
    covariates <- read_matrix(options[["--covariates"]])
    samples <- colnames(values)
    if (!setequal(samples, colnames(covariates))) stop(label, " and covariate sample sets differ")
    covariates <- covariates[, match(samples, colnames(covariates)), drop = FALSE]
    check_design(covariates)
    locations <- read_locations(options[[location_flag]], rownames(values))
    keep <- apply(values, 1L, function(values) {
        observed <- values[!is.na(values)]
        length(observed) >= 2L && max(observed) > min(observed)
    })
    cat(sprintf("%s %s: %d; removed all-missing/invariant rows: %d\n",
                label, settings[["trait"]], nrow(values), sum(!keep)))
    values <- values[keep, , drop = FALSE]
    locations <- locations[keep, , drop = FALSE]
    if (!nrow(values)) stop("no variable ", label, " phenotypes remain")
    groups <- observed_groups(values, covariates)
    cat(sprintf("remaining %s NA: %d; observed-sample groups: %d; phenotype NA are not imputed\n",
                label, sum(is.na(values)), length(groups)))
    loaded <- load_genotypes(options[["--vcf"]], samples, parameters[2L], parameters[3L],
                             counted_allele = counted_allele)
    threshold <- 0.05 / loaded$snps$nRows()
    cat(sprintf("MatrixEQTL %s; samples: %d; phenotypes: %d; covariates: %d; cis distance: %g bp\n",
                packageVersion("MatrixEQTL"), length(samples), nrow(values),
                nrow(covariates), parameters[1L]))
    cat(sprintf("cis/trans P threshold: 0.05 / %d = %.17g\n", loaded$snps$nRows(), threshold))
    cat("genotype dosage:", counted_allele, "allele count; beta is per counted allele\n")
    cat("remaining GT NA use per-variant observed means within each observed-sample group\n")
    cat(sprintf("phenotypes on chromosomes without retained variants: %d (trans tests only)\n",
                sum(!locations$chr %in% loaded$positions$chr)))
    working <- tempfile(paste0(".", settings[["task"]], "-"), tmpdir = out_dir)
    dir.create(working)
    on.exit(unlink(working, recursive = TRUE), add = TRUE)
    temporary <- file.path(working, output_names)
    result <- run_observed_models(loaded$snps, loaded$positions, values, covariates, locations,
                                  parameters[1L], threshold, id_column = id_column)
    fwrite(result$cis, temporary[1L], sep = "\t", quote = FALSE)
    fwrite(result$trans, temporary[2L], sep = "\t", quote = FALSE)
    selected <- rownames(values) %in% unique(result$cis[[id_column]])
    subset <- values[selected, , drop = FALSE]
    exported <- data.frame(trait_id = as.character(rownames(subset)), subset, check.names = FALSE)
    names(exported)[1L] <- id_column
    fwrite(exported, temporary[3L], sep = "\t", quote = FALSE, na = "NA")
    if (any(selected)) {
        run_observed_models(loaded$snps, loaded$positions, subset, covariates,
                            locations[selected, , drop = FALSE], parameters[1L], 1,
                            all_cis_file = temporary[4L], id_column = id_column)
    } else {
        table <- empty_associations(id_column)
        fwrite(table[, names(table) != "FDR"], temporary[4L], sep = "\t", quote = FALSE)
    }
    for (i in seq_along(outputs)) {
        if (!file.rename(temporary[i], outputs[i])) stop("could not write result: ", outputs[i])
    }
    cat(sprintf("significant cis pairs: %d; trans pairs: %d; cis %s: %d\n",
                nrow(result$cis), nrow(result$trans), settings[["trait"]], sum(selected)))
    cat("saved:", paste(outputs, collapse = "\n"), "\n")
}

if (sys.nframe() == 0L) main()
