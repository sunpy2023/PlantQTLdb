# map cis/trans-eQTLs from a QC VCF, complete normalized expression and covariates
# requires MatrixEQTL, data.table and bcftools; all sample IDs are retained as strings
# diploid biallelic SNPs and indels are read in chunks; GT counts the selected REF or ALT allele
# variant QC uses observed calls before imputation; IDs are CHROM:POS:REF:ALT
# gene/variant chromosomes must use the same reference and naming convention
# expression/covariate NA cause an error; invariant genes are removed without dropping samples
# MatrixEQTL fills remaining GT NA with each variant's observed mean
# significance threshold: 0.05 / retained variants, not the number of variant-gene pairs
# outputs: cis.significant.tsv, trans.significant.tsv, cis_egenes.expression.tsv,
# cis_egenes.all_cis.tsv; no cis-eGenes produces header-only eGene output tables

suppressPackageStartupMessages(library(MatrixEQTL))
suppressPackageStartupMessages(library(data.table))

missing_tokens <- c("", "NA", "NaN", "nan", "N/A", ".")

check_ids <- function(values, label) {
    if (anyNA(values) || any(values %in% missing_tokens) || anyDuplicated(values)
            || any(grepl("\\s", values, perl = TRUE))) {
        stop("invalid or duplicate ", label, call. = FALSE)
    }
}

read_table <- function(path, numeric_columns = FALSE) {
    line <- readLines(path, n = 1L, warn = FALSE)
    if (!length(line) || !nzchar(line) || endsWith(line, "\t")) stop("empty table/header: ", path)
    header <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    check_ids(header, paste("column names in", path))
    table <- withCallingHandlers(
        fread(path, sep = "\t", header = TRUE,
              colClasses = if (numeric_columns) list(character = 1L) else "character",
              na.strings = missing_tokens,
              check.names = FALSE, data.table = FALSE, showProgress = FALSE),
        warning = function(w) stop(conditionMessage(w), call. = FALSE))
    if (!identical(names(table), header) || nrow(table) < 1L) stop("invalid/empty table: ", path)
    table
}

read_matrix <- function(path, allow_missing = FALSE) {
    table <- read_table(path, numeric_columns = TRUE)
    if (ncol(table) < 2L) stop("matrix has no sample columns: ", path)
    check_ids(table[[1L]], paste("row IDs in", path))
    if (!all(vapply(table[-1L], function(column) is.numeric(column) || all(is.na(column)), logical(1)))) {
        stop("invalid numeric values: ", path)
    }
    values <- as.matrix(table[-1L])
    storage.mode(values) <- "double"
    if (any(is.infinite(values))) stop("non-finite numeric values: ", path)
    if (!allow_missing && anyNA(values)) stop("missing values in matrix: ", path)
    rownames(values) <- table[[1L]]
    values
}

check_design <- function(covariates) {
    design <- cbind(intercept = 1, t(covariates))
    scales <- apply(abs(design), 2L, max)
    scaled <- sweep(design, 2L, ifelse(scales == 0, 1, scales), "/")
    if (qr(scaled)$rank != ncol(design)) stop("covariates are collinear with an intercept")
    if (ncol(covariates) <= ncol(design) + 1L) stop("insufficient regression degrees of freedom")
}

read_locations <- function(path, ids) {
    locations <- read_table(path)
    if (!all(c("geneid", "chr", "s1", "s2") %in% names(locations))) {
        stop("phenotype locations require geneid, chr, s1 and s2 columns")
    }
    locations <- locations[, c("geneid", "chr", "s1", "s2")]
    check_ids(locations$geneid, "location IDs")
    check_ids(unique(locations$chr), "location chromosome names")
    if (!setequal(ids, locations$geneid)) stop("phenotype and location IDs differ")
    locations <- locations[match(ids, locations$geneid), ]
    for (column in c("s1", "s2")) {
        locations[[column]] <- suppressWarnings(as.numeric(locations[[column]]))
        if (any(!is.finite(locations[[column]]) | locations[[column]] < 1
                | locations[[column]] != floor(locations[[column]]))) stop("invalid phenotype coordinates")
    }
    if (any(locations$s1 > locations$s2)) stop("phenotype interval start exceeds end")
    locations
}

sliced_matrix <- function(values) {
    object <- SlicedData$new()
    object$CreateFromMatrix(values)
    object$ResliceCombined(2000L)
    object
}

run_bcftools <- function(args) {
    result <- system2(Sys.which("bcftools"), vapply(args, shQuote, character(1)), stdout = TRUE)
    status <- attr(result, "status")
    if (!is.null(status) && status != 0L) stop("bcftools failed with exit status ", status)
    result
}

gt_dosage <- function(gt) {
    gt <- gsub("|", "/", gt, fixed = TRUE)
    if (any(!grepl("^(?:[01.]/[01.]|\\.)$", gt, perl = TRUE))) {
        stop("expected diploid biallelic GT values, with phased/unphased or missing alleles")
    }
    unname(c("0/0" = 0L, "0/1" = 1L, "1/0" = 1L, "1/1" = 2L)[gt])
}

load_genotypes <- function(path, samples, min_maf, max_missing, chunk_size = 1000L,
                           counted_allele = "ALT") {
    if (!counted_allele %in% c("REF", "ALT")) stop("counted allele must be REF or ALT")
    vcf_samples <- run_bcftools(c("query", "-l", path))
    check_ids(vcf_samples, "VCF sample IDs")
    order <- match(samples, vcf_samples)
    if (anyNA(order)) stop("phenotype samples missing from VCF: ",
                          paste(samples[is.na(order)], collapse = ", "))
    format <- "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n"
    command <- paste(shQuote(Sys.which("bcftools")), "query -i", shQuote("N_ALT=1"),
                     "-f", shQuote(format), shQuote(path))
    connection <- pipe(command, "r")
    closed <- FALSE
    on.exit(if (!closed) try(close(connection), silent = TRUE), add = TRUE)
    object <- SlicedData$new()
    object$columnNames <- samples
    positions <- row_names <- list()
    queried <- retained_missing <- 0L
    repeat {
        lines <- readLines(connection, n = chunk_size, warn = FALSE)
        if (!length(lines)) break
        fields <- fread(text = paste0(paste(lines, collapse = "\n"), "\n"), sep = "\t",
                        header = FALSE, colClasses = "character", na.strings = NULL,
                        data.table = FALSE, showProgress = FALSE)
        if (ncol(fields) != length(vcf_samples) + 4L || anyNA(fields)) stop("malformed VCF query rows")
        queried <- queried + nrow(fields)
        gt <- as.matrix(fields[, order + 4L, drop = FALSE])
        dosage <- matrix(gt_dosage(as.vector(gt)), nrow(gt), ncol(gt))
        called <- rowSums(!is.na(dosage))
        alleles <- rowSums(dosage, na.rm = TRUE)
        maf <- pmin(alleles, 2 * called - alleles) / (2 * called)
        missing <- (ncol(dosage) - called) / ncol(dosage)
        variable <- called * rowSums(dosage^2, na.rm = TRUE) > alleles^2
        keep <- called > 0L & maf >= min_maf & missing <= max_missing & variable
        if (!any(keep)) next
        fields <- fields[keep, , drop = FALSE]
        position <- suppressWarnings(as.numeric(fields[[2L]]))
        if (any(!is.finite(position) | position < 1 | position != floor(position))) {
            stop("invalid VCF positions")
        }
        ids <- do.call(paste, c(fields[1:4], sep = ":"))
        i <- length(positions) + 1L
        retained <- dosage[keep, , drop = FALSE]
        if (counted_allele == "REF") retained <- 2L - retained
        object$setSlice(i, retained)
        row_names[[i]] <- ids
        positions[[i]] <- data.frame(snpid = ids, chr = fields[[1L]], pos = position)
        retained_missing <- retained_missing + sum(is.na(dosage[keep, , drop = FALSE]))
    }
    status <- close(connection)
    closed <- TRUE
    if (length(status) && status != 0L) stop("VCF query failed with exit status ", status)
    if (!length(positions)) stop("no variable biallelic variants pass MAF/missingness QC")
    object$rowNameSlices <- row_names
    positions <- as.data.frame(rbindlist(positions))
    check_ids(positions$snpid, "CHROM:POS:REF:ALT variant IDs")
    check_ids(unique(positions$chr), "variant chromosome names")
    cat(sprintf("VCF samples: %d; matched phenotype samples: %d; biallelic variants: %d; retained: %d\n",
                length(vcf_samples), length(samples), queried, nrow(positions)))
    cat(sprintf("variant QC: MAF >= %g; missingness <= %g; remaining GT NA: %d\n",
                min_maf, max_missing, retained_missing))
    list(snps = object, positions = positions)
}

main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    required <- c("--vcf", "--expression", "--gene-locations", "--covariates", "--out-dir")
    defaults <- c("--cis-distance" = "200000", "--min-maf" = "0.05", "--max-variant-missing" = "0.10",
                  "--counted-allele" = "ALT")
    usage <- paste("usage: Rscript run_eqtl_mapping.R --vcf genotype.qc.vcf.gz",
                   "--expression expression.normalized.tsv --gene-locations expression.locations.tsv",
                   "--covariates covariates.tsv --out-dir eqtl [--cis-distance 200000]",
                   "[--min-maf 0.05] [--max-variant-missing 0.10] [--counted-allele REF|ALT]")
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
                      "cis_egenes.expression.tsv", "cis_egenes.all_cis.tsv")
    outputs <- file.path(out_dir, output_names)
    if (any(vapply(outputs, normalizePath, character(1), mustWork = FALSE) %in% inputs)) {
        stop("output files must not overwrite inputs")
    }
    expression <- read_matrix(options[["--expression"]])
    covariates <- read_matrix(options[["--covariates"]])
    samples <- colnames(expression)
    if (!setequal(samples, colnames(covariates))) stop("expression and covariate sample sets differ")
    covariates <- covariates[, match(samples, colnames(covariates)), drop = FALSE]
    check_design(covariates)
    locations <- read_locations(options[["--gene-locations"]], rownames(expression))
    keep <- apply(expression, 1L, function(values) max(values) > min(values))
    cat(sprintf("expression genes: %d; removed invariant rows: %d\n",
                nrow(expression), sum(!keep)))
    expression <- expression[keep, , drop = FALSE]
    locations <- locations[keep, , drop = FALSE]
    if (!nrow(expression)) stop("no variable expression genes remain")
    loaded <- load_genotypes(options[["--vcf"]], samples, parameters[2L], parameters[3L],
                             counted_allele = counted_allele)
    threshold <- 0.05 / loaded$snps$nRows()
    cat(sprintf("MatrixEQTL %s; samples: %d; genes: %d; covariates: %d; cis distance: %g bp\n",
                packageVersion("MatrixEQTL"), length(samples), nrow(expression),
                nrow(covariates), parameters[1L]))
    cat(sprintf("cis/trans P threshold: 0.05 / %d = %.17g\n", loaded$snps$nRows(), threshold))
    cat("genotype dosage:", counted_allele, "allele count; beta is per counted allele\n")
    cat(sprintf("genes on chromosomes without retained variants: %d (trans tests only)\n",
                sum(!locations$chr %in% loaded$positions$chr)))
    working <- tempfile(".eqtl-", tmpdir = out_dir)
    dir.create(working)
    on.exit(unlink(working, recursive = TRUE), add = TRUE)
    temporary <- file.path(working, output_names)
    gene <- sliced_matrix(expression)
    cvrt <- sliced_matrix(covariates)
    result <- Matrix_eQTL_main(
        snps = loaded$snps, gene = gene, cvrt = cvrt, useModel = modelLINEAR,
        output_file_name = temporary[2L], pvOutputThreshold = threshold,
        output_file_name.cis = temporary[1L], pvOutputThreshold.cis = threshold,
        snpspos = loaded$positions, genepos = locations, cisDist = parameters[1L],
        pvalue.hist = FALSE, min.pv.by.genesnp = FALSE, noFDRsaveMemory = FALSE, verbose = TRUE)
    egenes <- unique(as.character(result$cis$eqtls$gene))
    selected <- rownames(expression) %in% egenes
    subset <- expression[selected, , drop = FALSE]
    fwrite(data.frame(gene_id = as.character(rownames(subset)), subset, check.names = FALSE),
           temporary[3L], sep = "\t", quote = FALSE, na = "NA")
    if (length(egenes)) {
        Matrix_eQTL_main(
            snps = loaded$snps, gene = sliced_matrix(subset), cvrt = cvrt, useModel = modelLINEAR,
            output_file_name = "", pvOutputThreshold = 0,
            output_file_name.cis = temporary[4L], pvOutputThreshold.cis = 1,
            snpspos = loaded$positions, genepos = locations[selected, , drop = FALSE],
            cisDist = parameters[1L], pvalue.hist = FALSE, min.pv.by.genesnp = FALSE,
            noFDRsaveMemory = TRUE, verbose = TRUE)
    } else {
        writeLines("SNP\tgene\tbeta\tt-stat\tp-value", temporary[4L])
    }
    for (i in seq_along(outputs)) {
        if (!file.rename(temporary[i], outputs[i])) stop("could not write result: ", outputs[i])
    }
    cat(sprintf("significant cis pairs: %d; trans pairs: %d; cis-eGenes: %d\n",
                result$cis$neqtls, result$trans$neqtls, length(egenes)))
    cat("saved:", paste(outputs, collapse = "\n"), "\n")
}

if (sys.nframe() == 0L) main()
