#!/usr/bin/env Rscript

# Keep implicit layout devices off disk; explicit PNG/PDF devices are unchanged.
local({
  .plot_devices_before <- grDevices::dev.list()
  .plot_device_option <- getOption("device")
  .rplots_path <- file.path(getwd(), "Rplots.pdf")
  on.exit({
    for (.plot_device in setdiff(grDevices::dev.list(), .plot_devices_before)) {
      try(grDevices::dev.off(.plot_device), silent = TRUE)
    }
    options(device = .plot_device_option)
    # Remove only the default-device artifact in this execution directory.
    unlink(c(.rplots_path, file.path(dirname(.rplots_path), "._Rplots.pdf")))
  }, add = TRUE)
  options(device = function(...) grDevices::pdf(file = NULL))
  grDevices::pdf(file = NULL)


# Render regulatory-target pies and gene-body profiles.

options(stringsAsFactors = FALSE)
options(scipen = 999)

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])[1]
script_file <- gsub("~+~", " ", script_file, fixed = TRUE)
if (is.na(script_file)) stop("Cannot resolve pipeline script path", call. = FALSE)
DEFAULT_RUN_ROOT <- normalizePath(file.path(dirname(script_file), "..", "..", ".."), mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L || length(args) > 4L) stop("Usage: render_targets.R run_root mode output_root [publication_figures]", call. = FALSE)
RUN_ROOT <- normalizePath(args[[1]], mustWork = TRUE)
ANALYSIS_MODE <- args[[2]]
if (!ANALYSIS_MODE %in% c("promoters", "all_genes", "all_peaks")) stop("Unknown analysis mode: ", ANALYSIS_MODE, call. = FALSE)
PEAK_MODE <- ANALYSIS_MODE == "all_peaks"
ID_COLUMN <- if (PEAK_MODE) "peak_id" else "gene"
UNIT <- if (PEAK_MODE) "peaks" else "genes"
PROFILE_PREFIX <- "Metaprofile"
GROUP_BASE_DIR <- normalizePath(args[[3]], mustWork = FALSE)
PUBLICATION_DIR <- if (length(args) == 4L) normalizePath(args[[4]], mustWork = FALSE) else NULL
COVERAGE_AXIS_LABEL <- "Coverage"
REGTARGET_BASE_DIR <- GROUP_BASE_DIR
GROUP_VISUAL_DIR <- GROUP_BASE_DIR
REGTARGET_VISUAL_DIR <- GROUP_BASE_DIR
GROUP_DIR <- file.path(GROUP_BASE_DIR, "data")
REGTARGET_DIR <- file.path(REGTARGET_BASE_DIR, "data")
GTF <- file.path(RUN_ROOT, "reference", "gencode.vM25.annotation.gtf")
COMPUTE_MATRIX <- Sys.which("computeMatrix")
PLOT_PROFILE <- Sys.which("plotProfile")
if (!nzchar(COMPUTE_MATRIX)) COMPUTE_MATRIX <- "/opt/anaconda3/envs/cutrun_env/bin/computeMatrix"
if (!nzchar(PLOT_PROFILE)) PLOT_PROFILE <- "/opt/anaconda3/envs/cutrun_env/bin/plotProfile"
BIGWIGS <- c(
  MCM3 = file.path(RUN_ROOT, "cutrun_work", "03_bigwig", "IGV_representation", "MCM3_mean.bw"),
  NONO = file.path(RUN_ROOT, "cutrun_work", "03_bigwig", "IGV_representation", "NONO_mean.bw"),
  PSPC1 = file.path(RUN_ROOT, "cutrun_work", "03_bigwig", "IGV_representation", "PSPC1_mean.bw")
)

target_term <- if (ANALYSIS_MODE == "promoters") "direct targets" else "shared-peak DEGs"
GROUPS <- list(
  A = list(label = "MCM3 + NONO + PSPC1", region = "MCM3_NONO_PSPC1", tracks = c("MCM3", "NONO", "PSPC1"), title = paste0("MCM3-PSPC1-NONO ", target_term, " (n=%s)")),
  B = list(label = "MCM3 + PSPC1 only", region = "MCM3_PSPC1_only", tracks = c("MCM3", "PSPC1"), title = paste0("MCM3-PSPC1-only ", target_term, " (n=%s)")),
  C = list(label = "NONO + PSPC1 only", region = "NONO_PSPC1_only", tracks = c("NONO", "PSPC1"), title = paste0("NONO-PSPC1-only ", target_term, " (n=%s)")),
  D = list(label = "MCM3 + NONO only", region = "MCM3_NONO_only", tracks = c("MCM3", "NONO"), title = paste0("MCM3-NONO-only ", target_term, " (n=%s)"))
)
if (ANALYSIS_MODE == "all_genes") {
  for (name in names(GROUPS)) {
    GROUPS[[name]]$title <- paste0("Group ", name, ": shared-peak-associated DEGs (n=%s)")
  }
}
if (PEAK_MODE) {
  for (name in names(GROUPS)) GROUPS[[name]]$title <- paste0("Group ", name, ": DEG-linked shared peaks (n=%s)")
}
COLORS <- c(MCM3 = "#cf111f", NONO = "#f2a010", PSPC1 = "#4f7db1")
DIR_COLORS <- c(Up = "#D87070", Down = "#6FA37A", Mix = "#F3D37A")
BIN_SIZE <- 25L
BEFORE_BP <- 3000L
BODY_BP <- 3000L
AFTER_BP <- 3000L
DISPLAY_SCALE <- 1

for (pkg in c("data.table", "ggplot2", "grid", "png", "rtracklayer", "GenomicRanges", "GenomeInfoDb", "S4Vectors")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing R package: ", pkg, call. = FALSE)
}
dir.create(GROUP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REGTARGET_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(GROUP_VISUAL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REGTARGET_VISUAL_DIR, recursive = TRUE, showWarnings = FALSE)
if (!is.null(PUBLICATION_DIR)) dir.create(PUBLICATION_DIR, recursive = TRUE, showWarnings = FALSE)
for (path in c(GROUP_DIR, GTF, COMPUTE_MATRIX, PLOT_PROFILE, unname(BIGWIGS))) {
  if (!file.exists(path)) stop("Missing required input: ", path, call. = FALSE)
}
Sys.setenv(MPLCONFIGDIR = file.path(tempdir(), "mcm3_final_mplconfig"))

read_genes <- function(path) {
  tab <- data.table::fread(path, data.table = FALSE)
  if (!ID_COLUMN %in% names(tab)) stop("Missing identifier column: ", path, call. = FALSE)
  values <- tab[[ID_COLUMN]]
  sort(unique(trimws(as.character(values))[nzchar(trimws(as.character(values))) & !is.na(values)]))
}

write_genes <- function(path, genes) {
  data.table::fwrite(setNames(data.frame(sort(unique(genes))), ID_COLUMN), path, sep = "\t", quote = FALSE)
}

class_genes <- function(factor, class_name) {
  read_genes(file.path(GROUP_DIR, sprintf("promoter_vs_DEG_direction_%s_%s.tsv", factor, class_name)))
}

direction <- function(gene, factor) {
  if (PEAK_MODE) return(peak_directions[[factor]][match(gene, peak_directions$peak_id)])
  if (gene %in% class_sets[[factor]]$up) return("Up")
  if (gene %in% class_sets[[factor]]$down) return("Down")
  stop("Direct target has no assigned direction: ", gene, " / ", factor, call. = FALSE)
}

summarize_group <- function(genes, factors) {
  rows <- lapply(sort(unique(genes)), function(gene) {
    dirs <- vapply(factors, function(factor) direction(gene, factor), character(1))
    cls <- if (all(dirs == "Up")) "Up" else if (all(dirs == "Down")) "Down" else "Mix"
    out <- data.frame(gene = gene, class = cls, stringsAsFactors = FALSE)
    names(out)[1] <- ID_COLUMN
    for (factor in factors) out[[factor]] <- dirs[[factor]]
    out
  })
  if (!length(rows)) return(data.frame(gene = character(), class = character(), stringsAsFactors = FALSE))
  do.call(rbind, rows)
}

pie_plot <- function(title, tab, show_numbers = TRUE, show_labels = TRUE, levels = c("Up", "Down", "Mix"), radial_small_labels = FALSE) {
  dat <- as.data.frame(table(factor(tab$class, levels = levels)), stringsAsFactors = FALSE)
  names(dat) <- c("class", "n")
  dat <- dat[dat$n > 0, , drop = FALSE]
  dat$frac <- dat$n / sum(dat$n)
  dat$label <- if (!show_labels) "" else if (show_numbers) paste0(dat$class, "\n(n=", dat$n, ")") else as.character(dat$class)
  if (PEAK_MODE && show_labels && show_numbers) {
    small <- dat$frac < 0.08
    dat$label[small] <- paste0(dat$class[small], " (n=", dat$n[small], ")")
  }
  dat$angle <- 0
  dat$text_size <- 4
  dat$text_x <- 1
  if (radial_small_labels && show_labels) {
    # Match ggplot's reverse-alphabetical stacking; turn narrow-slice labels radially.
    order <- order(dat$class, decreasing = TRUE)
    mid <- numeric(nrow(dat))
    mid[order] <- cumsum(dat$frac[order]) - dat$frac[order] / 2
    dat$mid <- mid
    small <- dat$frac < 0.04
    dat$angle[small] <- (90 - 360 * mid[small] + 90) %% 180 - 90
    dat$text_size[small] <- 2.2
    dat$text_x[small] <- 1.25
  }
  ggplot2::ggplot(dat, ggplot2::aes(x = "", y = frac, fill = class)) +
    ggplot2::geom_col(width = 1, color = "#222222", linewidth = 0.6) +
    (if (radial_small_labels && show_labels)
      ggplot2::geom_text(ggplot2::aes(x = text_x, y = mid, label = label, angle = angle, size = text_size), fontface = "bold")
     else ggplot2::geom_text(ggplot2::aes(label = label), position = ggplot2::position_stack(vjust = 0.5), size = 4, fontface = "bold")) +
    (if (radial_small_labels && show_labels) ggplot2::scale_size_identity() else NULL) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::scale_fill_manual(values = c(DIR_COLORS, "No DEG" = "#BFBFBF")) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5), plot.margin = ggplot2::margin(8, 8, 8, 8))
}

binding_pie <- function(factor, all_targets, direct_targets, out, show_numbers = TRUE) {
  all_targets <- sort(unique(all_targets))
  direct_targets <- sort(unique(direct_targets))
  if (!all(direct_targets %in% all_targets)) {
    stop("CUT&RUN-associated targets are not a subset of all DEG targets for ", factor, call. = FALSE)
  }
  bound_label <- if (ANALYSIS_MODE == "promoters") "Targets w/ promoter binding" else "Shared-peak-associated DEGs"
  unbound_label <- if (ANALYSIS_MODE == "promoters") "Targets w/o promoter binding" else "Other DEGs"
  if (PEAK_MODE) {
    bound_label <- "Shared peaks linked to DEGs"
    unbound_label <- "Shared peaks without a DEG link"
  }
  dat <- data.frame(
    class = c(bound_label, unbound_label),
    n = c(length(direct_targets), length(all_targets) - length(direct_targets)),
    stringsAsFactors = FALSE
  )
  dat$frac <- dat$n / sum(dat$n)
  dat$label <- if (show_numbers) paste0(sprintf("%.1f", 100 * dat$frac), "%\n(n=", dat$n, ")") else ""
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = "", y = frac, fill = class)) +
    ggplot2::geom_col(width = 1, color = "#222222", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = label), position = ggplot2::position_stack(vjust = 0.5), size = 4, fontface = "bold") +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::scale_fill_manual(values = setNames(c("#BFD9F2", "#9AA3C7"), c(bound_label, unbound_label))) +
    ggplot2::theme_void(base_size = 12)
  if (show_numbers) {
    p <- p + ggplot2::labs(
      title = sprintf(if (PEAK_MODE) "%s KD: shared peaks (n=%s)" else "%s targets (n=%s)", factor, format(length(all_targets), big.mark = ",")),
      subtitle = if (PEAK_MODE) "At least one assigned protein-coding gene is a DEG" else if (ANALYSIS_MODE == "promoters") "Promoter binding from CUT&RUN (TSS +/- 1 kb)" else "Gene association with three-factor shared CUT&RUN peaks"
    ) + ggplot2::theme(
      legend.position = "bottom", legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 11),
      plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5)
    )
  } else {
    p <- p + ggplot2::theme(legend.position = "none", plot.margin = ggplot2::margin(8, 8, 8, 8))
  }
  ggplot2::ggsave(out, p, width = 6.5, height = 5.5, dpi = 300, bg = "white")
}

read_gtf_genes <- function(path) {
  gr <- rtracklayer::import(path)
  gr <- gr[S4Vectors::mcols(gr)$type == "gene"]
  gene_values <- if (ANALYSIS_MODE == "promoters") {
    as.character(S4Vectors::mcols(gr)$gene_name)
  } else {
    sub("\\..*$", "", as.character(S4Vectors::mcols(gr)$gene_id))
  }
  data.frame(
    chrom = as.character(GenomeInfoDb::seqnames(gr)), start = as.integer(GenomicRanges::start(gr)) - 1L,
    end = as.integer(GenomicRanges::end(gr)), strand = as.character(GenomicRanges::strand(gr)),
    gene = gene_values, stringsAsFactors = FALSE
  )
}

write_gene_body_bed <- function(gtf_genes, genes, path, as_peaks = PEAK_MODE) {
  if (as_peaks) {
    dat <- gtf_genes[gtf_genes$locus_id %in% genes, c("chrom", "start", "end", "locus_id")]
    if (!setequal(dat$locus_id, genes) || anyDuplicated(dat$locus_id)) stop("Peak BED membership mismatch", call. = FALSE)
    dat$score <- 0L
    dat$strand <- "."
    data.table::fwrite(dat, path, sep = "\t", col.names = FALSE, quote = FALSE)
    return(nrow(dat))
  }
  dat <- gtf_genes[gtf_genes$gene %in% genes, , drop = FALSE]
  missing <- setdiff(genes, unique(dat$gene))
  if (length(missing)) stop("Genes absent from GTF: ", paste(missing, collapse = ", "), call. = FALSE)
  dat <- dat[dat$end > dat$start & dat$strand %in% c("+", "-"), c("chrom", "start", "end", "gene", "strand")]
  dat$score <- 0L
  dat <- dat[, c("chrom", "start", "end", "gene", "score", "strand")]
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  write.table(dat, tmp, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  if (file.exists(path) && !file.remove(path)) stop("Cannot replace BED: ", path, call. = FALSE)
  if (!file.rename(tmp, path)) stop("Cannot install BED: ", path, call. = FALSE)
  nrow(dat)
}

run_command <- function(exe, args) {
  status <- system2(exe, args = vapply(as.character(args), shQuote, character(1)))
  if (!identical(status, 0L)) stop("Command failed: ", exe, call. = FALSE)
}

read_profile <- function(path, sample_order, peak_centered = PEAK_MODE) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines)) & !grepl("^\\s*#", lines)]
  lines <- lines[!(tolower(trimws(lines)) %in% c("bin labels", "bins"))]
  rows <- list()
  for (line in lines) {
    parts <- strsplit(trimws(line), "[[:space:]]+", perl = TRUE)[[1]]
    if (length(parts) < 3 || tolower(parts[[1]]) %in% c("bin", "bins")) next
    nums <- suppressWarnings(as.numeric(parts[-c(1, 2)]))
    nums <- nums[is.finite(nums)]
    if (!length(nums)) next
    rows[[length(rows) + 1L]] <- list(sample = parts[[1]], region = parts[[2]], values = nums)
  }
  if (!length(rows)) stop("Cannot parse plotProfile table: ", path, call. = FALSE)
  max_bins <- max(vapply(rows, function(row) length(row$values), integer(1)))
  parsed <- do.call(rbind, lapply(rows, function(row) {
    values <- c(row$values, rep(NA_real_, max_bins - length(row$values)))
    data.frame(sample = row$sample, region = row$region, t(values), check.names = FALSE)
  }))
  bin_cols <- names(parsed)[-(1:2)]
  out <- data.table::melt(data.table::as.data.table(parsed), id.vars = c("sample", "region"), measure.vars = bin_cols, variable.name = "bin", value.name = "signal")
  out[, bin_index := match(bin, bin_cols)]
  # Display the full -3 kb / scaled gene body / +3 kb span on a -3 to +3 axis.
  # TSS and TES are consequently internal landmarks at -1 and +1, respectively.
  out[, pos_kb := -3 + (bin_index - 1L) * (6 / (length(bin_cols) - 1L))]
  if (peak_centered) out[, pos_kb := (-BEFORE_BP + (bin_index - 0.5) * BIN_SIZE) / 1000]
  as.data.frame(out[, .(signal = mean(signal, na.rm = TRUE)), by = .(sample, region, pos_kb)])
}

make_group_profile <- function(group_name, genes, tracks, gtf_genes) {
  prefix <- file.path(REGTARGET_DIR, paste0(PROFILE_PREFIX, "_Group", group_name))
  bed <- paste0(prefix, "_regions.bed")
  n_bed <- write_gene_body_bed(gtf_genes, genes, bed)
  matrix <- paste0(prefix, ".matrix.gz")
  profile <- paste0(prefix, "_data.tsv")
  raw_plot <- tempfile(pattern = paste0("metaprofile_", group_name, "_"), fileext = ".png")
  on.exit(unlink(raw_plot), add = TRUE)
  if (PEAK_MODE) {
    run_command(COMPUTE_MATRIX, c("reference-point", "--referencePoint", "center", "-S", unname(BIGWIGS[tracks]), "-R", bed, "-b", BEFORE_BP, "-a", AFTER_BP, "--binSize", BIN_SIZE, "--missingDataAsZero", "--samplesLabel", tracks, "-o", matrix))
  } else {
    run_command(COMPUTE_MATRIX, c("scale-regions", "-S", unname(BIGWIGS[tracks]), "-R", bed, "--beforeRegionStartLength", BEFORE_BP, "--regionBodyLength", BODY_BP, "--afterRegionStartLength", AFTER_BP, "--binSize", BIN_SIZE, "--missingDataAsZero", "--samplesLabel", tracks, "-o", matrix))
  }
  run_command(PLOT_PROFILE, c("-m", matrix, "--plotType", "lines", "--legendLocation", "upper-right", "--plotTitle", sprintf(GROUPS[[group_name]]$title, n_bed), "--plotFileFormat", "png", "--plotHeight", "6", "--plotWidth", "9", "--outFileNameData", profile, "-out", raw_plot))
  long <- read_profile(profile, tracks)
  long$signal <- long$signal * DISPLAY_SCALE
  long$sample <- factor(long$sample, levels = tracks)
  yy <- long$signal[is.finite(long$signal)]
  ymax <- max(yy, na.rm = TRUE) * 1.05
  ymin <- 0
  p <- ggplot2::ggplot(long, ggplot2::aes(pos_kb, signal, color = sample)) +
    ggplot2::geom_line(linewidth = 1.8, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = if (PEAK_MODE) 0 else c(-1, 1), linetype = 2, linewidth = 1, color = "#222222") +
    ggplot2::scale_color_manual(values = COLORS[tracks], breaks = tracks, name = NULL) +
    ggplot2::scale_x_continuous(limits = c(-3, 3), breaks = if (PEAK_MODE) c(-3, 0, 3) else c(-3, -1, 1, 3), labels = if (PEAK_MODE) c("-3", "Peak center", "3") else c("-3", "TSS", "TES", "3")) +
    ggplot2::scale_y_continuous(limits = c(ymin, ymax)) +
    ggplot2::labs(title = sprintf(GROUPS[[group_name]]$title, n_bed), x = "Relative distance (kb)", y = COVERAGE_AXIS_LABEL) +
    ggplot2::theme_classic(base_size = 15) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 17, face = "bold", hjust = 0.5), axis.line = ggplot2::element_line(linewidth = 1.0), panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1.0), axis.ticks = ggplot2::element_line(linewidth = 1.0), axis.text = ggplot2::element_text(size = 14, face = "bold", color = "black"), axis.title = ggplot2::element_text(size = 16, face = "bold"), legend.position = "right", legend.text = ggplot2::element_text(size = 14, face = "bold"), plot.margin = ggplot2::margin(10, 12, 12, 12))
  out_png <- file.path(REGTARGET_VISUAL_DIR, paste0(PROFILE_PREFIX, "_Group", group_name, ".png"))
  # Match the single-panel Profile_MCM3_NONO aspect ratio (7.4:5.6).
  ggplot2::ggsave(out_png, p, width = 10.5, height = 10.5 * 5.6 / 7.4, dpi = 300, bg = "white")
  if (!is.null(PUBLICATION_DIR)) {
    p_no_text <- p + ggplot2::theme(
      plot.title = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
      legend.position = "none"
    )
    ggplot2::ggsave(file.path(PUBLICATION_DIR, paste0(PROFILE_PREFIX, "_Group", group_name, "_noTexts.png")), p_no_text,
                    width = 10.5, height = 10.5 * 5.6 / 7.4, dpi = 300, bg = "white")
  }
  for (path in c(out_png, profile, bed)) {
    target <- file.path(if (grepl("\\.png$", path)) GROUP_VISUAL_DIR else GROUP_DIR, basename(path))
    if (!identical(normalizePath(path, mustWork = FALSE), normalizePath(target, mustWork = FALSE))) file.copy(path, target, overwrite = TRUE)
  }
  summary <- data.frame(group = group_name, definition = GROUPS[[group_name]]$label, genes = length(genes), bed_rows = n_bed, tracks = paste(tracks, collapse = ","), stringsAsFactors = FALSE)
  if (PEAK_MODE) names(summary)[names(summary) == "genes"] <- "peaks"
  summary
}

make_mcm3_target_profile <- function(gtf_genes) {
  up <- class_sets$MCM3$up
  down <- class_sets$MCM3$down
  prefix <- file.path(REGTARGET_DIR, paste0(PROFILE_PREFIX, "_MCM3_target"))
  if (PEAK_MODE) {
    evidence <- data.table::fread(file.path(GROUP_DIR, "PeakGeneDEGEvidence.tsv"), data.table = FALSE)
    selected <- evidence[(evidence$peak_id %in% up & evidence$MCM3 == "Up") |
                         (evidence$peak_id %in% down & evidence$MCM3 == "Down"), , drop = FALSE]
    up <- sort(unique(selected$gene_id[selected$MCM3 == "Up"]))
    down <- sort(unique(selected$gene_id[selected$MCM3 == "Down"]))
    if (length(intersect(up, down))) stop("MCM3 gene direction sets overlap", call. = FALSE)
    data.table::fwrite(selected, paste0(prefix, "_gene_peak_assignments.tsv"), sep = "\t", quote = FALSE)
    data.table::fwrite(unique(selected[, c("gene_id", "gene", "gene_type", "MCM3")]),
                      paste0(prefix, "_genes.tsv"), sep = "\t", quote = FALSE)
    excluded <- data.frame(peak_id = class_sets$MCM3$mix,
                           reason = rep("Assigned MCM3 DEGs include both Up and Down", length(class_sets$MCM3$mix)))
    write.table(excluded, paste0(prefix, "_excluded_peaks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    gtf_genes <- read_gtf_genes(GTF)
  }
  up_bed <- file.path(REGTARGET_DIR, paste0(PROFILE_PREFIX, "_MCM3_target_UP_regions.bed"))
  down_bed <- file.path(REGTARGET_DIR, paste0(PROFILE_PREFIX, "_MCM3_target_DOWN_regions.bed"))
  n_up <- write_gene_body_bed(gtf_genes, up, up_bed, as_peaks = FALSE)
  n_down <- write_gene_body_bed(gtf_genes, down, down_bed, as_peaks = FALSE)
  matrix <- paste0(prefix, ".matrix.gz")
  profile <- paste0(prefix, "_data.tsv")
  raw_plot <- tempfile(pattern = "mcm3_target_metaprofile_", fileext = ".png")
  on.exit(unlink(raw_plot), add = TRUE)
  # This deepTools build has no computeMatrix --regionsLabel option. The two
  # region files retain their order and are relabelled explicitly below.
  beds <- c(up_bed, down_bed)
  n_regions <- c(n_up, n_down)
  direction_names <- c("MCM3 up-regulated targets", "MCM3 down-regulated targets")
  direction_colors <- c("#2855B2", "#74E37D")
  run_command(COMPUTE_MATRIX, c("scale-regions", "-S", BIGWIGS[["MCM3"]], "-R", beds, "-b", BEFORE_BP, "-a", AFTER_BP, "--binSize", BIN_SIZE, "--regionBodyLength", BODY_BP, "--missingDataAsZero", "--startLabel", "TSS", "--endLabel", "TES", "--samplesLabel", "MCM3", "-o", matrix))
  run_command(PLOT_PROFILE, c("-m", matrix, "--plotType", "lines", "--legendLocation", "upper-right", "--plotTitle", "MCM3 CUT&RUN across gene bodies (TSS to TES; scaled)", "--plotFileFormat", "png", "--plotHeight", "6", "--plotWidth", "9", "--outFileNameData", profile, "-out", raw_plot))
  long <- read_profile(profile, "MCM3", peak_centered = FALSE)
  long$signal <- long$signal * DISPLAY_SCALE
  region_levels <- unique(long$region)
  if (length(region_levels) != length(direction_names)) stop("Unexpected MCM3 target region classes.", call. = FALSE)
  long$direction <- direction_names[match(long$region, region_levels)]
  yy <- long$signal[is.finite(long$signal)]
  ymax <- max(yy, na.rm = TRUE) * 1.05
  ymin <- 0
  p <- ggplot2::ggplot(long, ggplot2::aes(pos_kb, signal, color = direction)) +
    ggplot2::geom_rect(
      data = data.frame(xmin = -1, xmax = 1, ymin = -Inf, ymax = Inf),
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = "#F2F2F2", color = NA
    ) +
    ggplot2::geom_line(linewidth = 1.8, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = 0.9, color = "#444444") +
    ggplot2::scale_color_manual(
      values = setNames(direction_colors, direction_names), breaks = direction_names, name = NULL
    ) +
    ggplot2::scale_x_continuous(limits = c(-3, 3), breaks = c(-3, -1, 1, 3), labels = c("-3.0", "TSS", "TES", "3.0")) +
    ggplot2::scale_y_continuous(limits = c(ymin, ymax)) +
    ggplot2::labs(title = "MCM3 CUT&RUN across gene bodies (TSS\u2192TES; scaled)", x = "Relative distance (kb)", y = COVERAGE_AXIS_LABEL) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 21, face = "bold", hjust = 0.5, margin = ggplot2::margin(b = 10)), axis.line = ggplot2::element_line(linewidth = 1), panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 1.0), axis.ticks = ggplot2::element_line(linewidth = 1), axis.text = ggplot2::element_text(size = 15, face = "bold", color = "black"), axis.title = ggplot2::element_text(size = 18, face = "bold"), legend.position = "top", legend.direction = "horizontal", legend.text = ggplot2::element_text(size = 15, face = "bold"), legend.key.width = grid::unit(1.1, "cm"), plot.margin = ggplot2::margin(18, 18, 18, 24))
  out_png <- file.path(REGTARGET_VISUAL_DIR, paste0(PROFILE_PREFIX, "_MCM3_target.png"))
  ggplot2::ggsave(out_png, p, width = 14.6667, height = 14.6667 * 5.6 / 7.4, dpi = 300, bg = "white")
  if (!is.null(PUBLICATION_DIR)) {
    p_no_text <- p + ggplot2::theme(
      plot.title = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
      legend.position = "none"
    )
    ggplot2::ggsave(file.path(PUBLICATION_DIR, paste0(PROFILE_PREFIX, "_MCM3_target_noTexts.png")), p_no_text,
                    width = 14.6667, height = 14.6667 * 5.6 / 7.4, dpi = 300, bg = "white")
  }
  target <- file.path(GROUP_VISUAL_DIR, basename(out_png))
  if (!identical(normalizePath(out_png, mustWork = FALSE), normalizePath(target, mustWork = FALSE))) file.copy(out_png, target, overwrite = TRUE)
  summary <- data.frame(direction = c("Up", "Down"), n_genes = n_regions)
  if (PEAK_MODE) summary$n_source_peaks <- c(length(class_sets$MCM3$up), length(class_sets$MCM3$down))
  write.table(summary, paste0(prefix, "_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

direct_sets <- setNames(lapply(c("MCM3", "NONO", "PSPC1"), function(factor) read_genes(file.path(GROUP_DIR, paste0("Venn_target_", factor, "_", UNIT, ".tsv")))), c("MCM3", "NONO", "PSPC1"))
if (PEAK_MODE) {
  peak_directions <- data.table::fread(file.path(GROUP_DIR, "Peak_DEG_directions.tsv"), data.table = FALSE)
  class_sets <- setNames(lapply(names(direct_sets), function(factor) {
    status <- peak_directions[[factor]]
    list(up = peak_directions$peak_id[status == "Up"], down = peak_directions$peak_id[status == "Down"],
         mix = peak_directions$peak_id[status == "Mix"], indirect = peak_directions$peak_id[status == "No DEG"])
  }), names(direct_sets))
  for (factor in names(direct_sets)) {
    tab <- data.frame(class = peak_directions[[factor]][peak_directions[[factor]] != "No DEG"])
    lev <- c("Up", "Down", "Mix")
    plot <- pie_plot(sprintf("%s KD: DEG-linked peaks (n=%s)", factor, nrow(tab)), tab,
                     levels = lev, radial_small_labels = TRUE) +
      ggplot2::labs(subtitle = "Direction of assigned-gene differential expression") +
      ggplot2::theme(plot.subtitle = ggplot2::element_text(face = "bold", hjust = .5, size = 11))
    ggplot2::ggsave(file.path(GROUP_VISUAL_DIR, paste0("Pie_DEGAssociation_", factor, ".png")), plot, width = 6.5, height = 5.5, dpi = 300, bg = "white")
    if (!is.null(PUBLICATION_DIR)) {
      bare <- pie_plot(NULL, tab, show_labels = FALSE, levels = lev)
      ggplot2::ggsave(file.path(PUBLICATION_DIR, paste0("Pie_DEGAssociation_", factor, "_noTexts.png")), bare, width = 6.5, height = 5.5, dpi = 300, bg = "white")
    }
  }
} else {
  class_sets <- setNames(lapply(c("MCM3", "NONO", "PSPC1"), function(factor) list(up = class_genes(factor, "A_bound_UP"), down = class_genes(factor, "B_bound_DOWN"), indirect = class_genes(factor, "D_DEG_only_indirect"))), c("MCM3", "NONO", "PSPC1"))
}

groups <- list(
  A = Reduce(intersect, direct_sets),
  B = setdiff(intersect(direct_sets$MCM3, direct_sets$PSPC1), direct_sets$NONO),
  C = setdiff(intersect(direct_sets$NONO, direct_sets$PSPC1), direct_sets$MCM3),
  D = setdiff(intersect(direct_sets$MCM3, direct_sets$NONO), direct_sets$PSPC1)
)
regions <- data.table::fread(file.path(GROUP_DIR, paste0("Venn_target_", UNIT, ".tsv")), data.table = FALSE)
summary_rows <- list()
pie_paths <- character()
pie_no_text_paths <- character()
for (name in names(GROUPS)) {
  expected <- sort(unique(regions[[ID_COLUMN]][regions$region == GROUPS[[name]]$region]))
  if (!setequal(groups[[name]], expected)) stop("Venn mismatch for group ", name, call. = FALSE)
  tab <- summarize_group(groups[[name]], GROUPS[[name]]$tracks)
  write_genes(file.path(GROUP_DIR, paste0("Group", name, "_", UNIT, ".tsv")), groups[[name]])
  write.table(tab, file.path(GROUP_DIR, paste0("Group", name, "_directions.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  counts <- table(factor(tab$class, levels = c("Up", "Down", "Mix")))
  summary_rows[[name]] <- data.frame(group = name, group_definition = GROUPS[[name]]$label, venn_region = GROUPS[[name]]$region, n_genes = nrow(tab), n_up = counts[["Up"]], n_down = counts[["Down"]], n_mix = counts[["Mix"]], stringsAsFactors = FALSE)
  if (PEAK_MODE) names(summary_rows[[name]])[4] <- "n_peaks"
  plot <- pie_plot(sprintf(if (PEAK_MODE) "Group %s: %s\nn = %s peaks" else "Group %s: %s\nn = %s", name, GROUPS[[name]]$label, nrow(tab)), tab)
  out <- file.path(GROUP_VISUAL_DIR, paste0("Pie_Group", name, ".png"))
  ggplot2::ggsave(out, plot, width = 5.2, height = 5.2, dpi = 300, bg = "white")
  pie_paths <- c(pie_paths, out)
  if (!is.null(PUBLICATION_DIR)) {
    plot_no_text <- pie_plot(NULL, tab, show_numbers = FALSE, show_labels = FALSE)
    out_no_text <- file.path(PUBLICATION_DIR, paste0("Pie_Group", name, "_noTexts.png"))
    ggplot2::ggsave(out_no_text, plot_no_text, width = 5.2, height = 5.2, dpi = 300, bg = "white")
    pie_no_text_paths <- c(pie_no_text_paths, out_no_text)
  }
}
write.table(do.call(rbind, summary_rows), file.path(GROUP_DIR, "Pie_Groups_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
definitions <- data.frame(group = names(GROUPS), group_definition = vapply(GROUPS, `[[`, character(1), "label"), venn_region = vapply(GROUPS, `[[`, character(1), "region"), n_genes = vapply(groups, length, integer(1)))
if (PEAK_MODE) names(definitions)[4] <- "n_peaks"
write.table(definitions, file.path(GROUP_DIR, "Pie_Groups_definitions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

grDevices::png(file.path(GROUP_VISUAL_DIR, "Pie_Groups.png"), width = 3150, height = 3300, res = 300, bg = "white")
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(3, 2, heights = grid::unit(c(0.6, 5.2, 5.2), "in"), widths = grid::unit(c(5.25, 5.25), "in"))))
grid::grid.text(if (PEAK_MODE) "DEG-linked peak groups (Up / Down / Mix)" else "Regulatory-target overlap groups (Up / Down / Mix)", vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2), gp = grid::gpar(fontsize = 18, fontface = "bold"))
for (i in seq_along(pie_paths)) {
  image <- png::readPNG(pie_paths[[i]])
  grid::pushViewport(grid::viewport(layout.pos.row = if (i <= 2) 2 else 3, layout.pos.col = if (i %% 2 == 1) 1 else 2))
  grid::grid.raster(image, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"), interpolate = TRUE)
  grid::popViewport()
}
grid::popViewport()
grDevices::dev.off()
if (!is.null(PUBLICATION_DIR)) {
  grDevices::png(file.path(PUBLICATION_DIR, "Pie_Groups_noTexts.png"), width = 3150, height = 3000, res = 300, bg = "white")
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2)))
  for (i in seq_along(pie_no_text_paths)) {
    image <- png::readPNG(pie_no_text_paths[[i]])
    grid::pushViewport(grid::viewport(layout.pos.row = if (i <= 2) 1 else 2, layout.pos.col = if (i %% 2 == 1) 1 else 2))
    grid::grid.raster(image, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"), interpolate = TRUE)
    grid::popViewport()
  }
  grid::popViewport()
  grDevices::dev.off()
}

input_parameters <- readLines(file.path(GROUP_DIR, "AnalysisParameters.tsv"), warn = FALSE)
writeLines(c(input_parameters, "mix_color\t#F3D37A (yellow)"), file.path(GROUP_DIR, "Pie_Groups_parameters.tsv"))

binding_summary <- list()
for (factor in c("NONO", "PSPC1")) {
  all_targets <- union(union(class_sets[[factor]]$up, class_sets[[factor]]$down), class_sets[[factor]]$indirect)
  if (PEAK_MODE) all_targets <- peak_directions$peak_id
  if (!all(direct_sets[[factor]] %in% all_targets)) {
    stop("Target-binding pie inputs are inconsistent for ", factor, call. = FALSE)
  }
  binding_summary[[factor]] <- data.frame(
    factor = factor,
    DEG_targets = length(unique(all_targets)),
    CUTRUN_associated_targets = length(unique(direct_sets[[factor]])),
    DEG_only_targets = length(setdiff(all_targets, direct_sets[[factor]])),
    CUTRUN_associated_percent = 100 * length(unique(direct_sets[[factor]])) / length(unique(all_targets)),
    stringsAsFactors = FALSE
  )
  if (PEAK_MODE) names(binding_summary[[factor]]) <- c("factor", "total_shared_peaks", "DEG_linked_peaks", "peaks_without_DEG_link", "DEG_linked_percent")
  binding_pie(factor, all_targets, direct_sets[[factor]], file.path(GROUP_VISUAL_DIR, paste0("Pie_TargetBinding_", factor, ".png")))
  if (!is.null(PUBLICATION_DIR)) binding_pie(factor, all_targets, direct_sets[[factor]],
    file.path(PUBLICATION_DIR, paste0("Pie_TargetBinding_", factor, "_noTexts.png")), show_numbers = FALSE)
}
write.table(do.call(rbind, binding_summary), file.path(GROUP_DIR, "Pie_TargetBinding_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

gtf_genes <- if (PEAK_MODE) data.table::fread(file.path(GROUP_DIR, "SharedPeakLoci.tsv"), data.table = FALSE) else read_gtf_genes(GTF)
profile_summary <- do.call(rbind, lapply(names(GROUPS), function(name) make_group_profile(name, groups[[name]], GROUPS[[name]]$tracks, gtf_genes)))
write.table(profile_summary, file.path(GROUP_DIR, paste0(PROFILE_PREFIX, "_Groups_summary.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
make_mcm3_target_profile(gtf_genes)
message("[DONE] Regulatory-target figures written to: ", GROUP_DIR)
})
