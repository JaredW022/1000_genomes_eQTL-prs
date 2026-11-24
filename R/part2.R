read_height_gwas <- function(file = "data/raw/Height.gwas.txt.gz") {
  gwas <- read.delim(file)

  # Keep SNPs with stronger evidence of association with height and reliable imputation quality
  keep <- gwas$P <= 0.001 & gwas$INFO > 0.8

  candidates <- gwas[
    which(keep),
    c("CHR", "BP", "SNP", "A1", "A2", "P", "OR")
  ]
  candidates$weight <- log(candidates$OR)

  candidates
}

read_height_bim <- function(file) {
  read.table(
    file,
    col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"),
    colClasses = c("integer", "character", "NULL", "integer",
                   "character", "character")
  )
}

match_height_positions <- function(gwas, reference) {
  positions <- list(gwas = gwas$BP, reference = reference$BP)
  unique_positions <- lapply(positions, function(x) {
    x[!duplicated(x) & !duplicated(x, fromLast = TRUE)]
  })
  shared <- Reduce(intersect, unique_positions)
  gwas_rows <- match(shared, positions$gwas)
  reference_rows <- match(shared, positions$reference)

  matches <- data.frame(
    gwas[gwas_rows, ],
    reference_row = reference_rows,
    reference_snp = reference$SNP[reference_rows],
    reference_a1 = reference$A1[reference_rows],
    reference_a2 = reference$A2[reference_rows]
  )
  reference_matches <-
    (matches$reference_a1 == matches$A1 & matches$reference_a2 == matches$A2) |
    (matches$reference_a1 == matches$A2 & matches$reference_a2 == matches$A1)
  matches[which(reference_matches), , drop = FALSE]
}

clump_height_chromosome <- function(variants, chromosome, reference_prefix,
                                    held_out_sample, output_dir, plink_executable) {
  output_prefix <- file.path(output_dir, paste0("chr", chromosome))
  clump_input <- paste0(output_prefix, "_clump_input.tsv")
  remove_file <- paste0(output_prefix, "_holdout.tsv")

  write.table(
    data.frame(ID = variants$reference_snp, P = variants$P),
    clump_input,
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    held_out_sample, remove_file, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE
  )

  status <- system2(plink_executable, c(
    "--bfile", shQuote(reference_prefix),
    "--remove", shQuote(remove_file),
    "--clump", shQuote(clump_input),
    "--clump-p1", "0.001",
    "--clump-r2", "0.1",
    "--clump-kb", "250",
    "--out", shQuote(output_prefix)
  ))

  if (status != 0) {
    stop("PLINK clumping failed on chromosome ", chromosome)
  }

  clump_file <- paste0(output_prefix, ".clumps")
  if (!file.exists(clump_file)) {
    return(variants[FALSE, , drop = FALSE])
  }
  clumps <- read.delim(clump_file)
  lead_rows <- match(clumps$ID, variants$reference_snp)
  variants[lead_rows, , drop = FALSE]
}

score_height_chromosome <- function(clumped, genotypes, sample_ids) {
  effect_copies <- genotypes[, clumped$reference_row, drop = FALSE]
  # BEDMatrix counts reference A1 copies.
  effect_is_a2 <- clumped$reference_a2 == clumped$A1
  effect_copies[, effect_is_a2] <- 2 - effect_copies[, effect_is_a2, drop = FALSE]
  setNames(as.numeric(effect_copies %*% clumped$weight), sample_ids)
}

run_height_prs <- function(output_dir, chromosomes = 1:22, sample_id = "HG00096",
                           project_root = ".") {
  if (dir.exists(output_dir)) {
    stop("Choose a new output directory for the clumping results.")
  }

  gwas <- read_height_gwas(file.path(project_root, "data/raw/Height.gwas.txt.gz"))
  plink_executable <- file.path(project_root, "bin/plink2.exe")
  totals <- NULL
  n_snps <- 0
  dir.create(output_dir, recursive = TRUE)

  for (chromosome in chromosomes) {
    reference_prefix <- file.path(
      project_root, "data/processed/LDREF", paste0("1000G.EUR.", chromosome)
    )
    samples <- read.table(
      paste0(reference_prefix, ".fam"),
      colClasses = c("character", "character", rep("NULL", 4))
    )
    sample_ids <- samples[[2]]
    sample_row <- which(sample_ids == sample_id)
    if (length(sample_row) != 1L) {
      stop("Expected one sample with ID ", sample_id, " on chromosome ", chromosome, ".")
    }

    reference <- read_height_bim(paste0(reference_prefix, ".bim"))
    matches <- match_height_positions(gwas[gwas$CHR == chromosome, ], reference)
    if (nrow(matches) == 0L) {
      next
    }

    genotypes <- BEDMatrix::BEDMatrix(paste0(reference_prefix, ".bed"), simple_names = TRUE)
    variants <- matches[!is.na(genotypes[sample_row, matches$reference_row]), , drop = FALSE]
    if (nrow(variants) == 0L) {
      next
    }
    clumped <- clump_height_chromosome(
      variants, chromosome, reference_prefix, samples[sample_row, ], output_dir,
      plink_executable
    )
    if (nrow(clumped) == 0L) {
      next
    }
    scores <- score_height_chromosome(clumped, genotypes, sample_ids)

    if (is.null(totals)) {
      totals <- scores
    } else {
      totals <- totals + scores[names(totals)]
    }
    n_snps <- n_snps + nrow(clumped)
  }

  if (n_snps == 0L) {
    stop("No variants were available for scoring.")
  }
  list(
    reference_scores = totals[names(totals) != sample_id],
    sample_score = totals[[sample_id]],
    n_snps = n_snps
  )
}
