######################### Plotting Coverge and Insertion data for specific regions ----
# library(dplyr)
# library(tidyverse)
# library(plyranges)
# library(ggrepel)
# library(ggbio)
# library(OrganismDbi)
# library(doParallel)



##### supporting functions for plotRegion
# Modified Track Plots

#' Get gene promoter range string
#'
#' For a given gene, return the gene position from an ArchR project
#' @param proj An ArchR Project
#' @param gene The name of a gene
#' @param upstream The number of basepairs to extend the range upstream of the gene promoter
#' @param downstream The number of basepairs to extend the range downstream of the gene promoter
#' @param output_string Logical value, default TRUE. If TRUE outputs a range string (ie 'chr1: 1700000-1710000') else
#' outputs a granges object
#' @return A range string or granges object, depending on 'output_string' parameter.
get_promoter_range <- function(proj, gene, upstream = 50000, downstream = 50000, output_string = TRUE) {
  gene_pos <- ArchR::getGenes(proj, gene)
  pro_pos_gr <- GenomicRanges::promoters(gene_pos, upstream = upstream, downstream = downstream)
  if (output_string) {
    pro_pos_str <- sprintf(
      "%s: %s-%s",
      GenomicRanges::seqnames(pro_pos_gr),
      GenomicRanges::start(GenomicRanges::ranges(pro_pos_gr)),
      GenomicRanges::end(GenomicRanges::ranges(pro_pos_gr))
    )
    return(pro_pos_str)
  } else {
    return(pro_pos_gr)
  }
}

#' CountDF to Region GRanges
#'
#' Gets the range covered by counts in a count df and converts to granges
#' Used in `plotRegion()`
#'
#' @param countdf  A dataframe that comes from `getbpCounts()` or `getbpInserts()`. Expected columns "chr" and "Locus"
#' @return A granges object for the region defined in countdf
countdf_to_region <- function(countdf) {
  assertthat::assert_that("chr" %in% names(countdf))
  assertthat::assert_that("Locus" %in% names(countdf))

  chrom <- toString(unique(countdf$chr))
  startSite <- min(countdf$Locus)
  endSite <- max(countdf$Locus)
  regionGRanges <- GenomicRanges::GRanges(
    seqnames = chrom,
    ranges = IRanges::IRanges(start = startSite, end = endSite), strand = "*"
  )

  return(regionGRanges)
}

#' Get GeneBody GRanges for a given region
#'
#' Used in `plotRegion()`
#'
#' @param regionGRanges regionGRanges A region Granges object to retrieve gene bodies for. For example, the output of countdf_to_region.
#' @param TxDb A TxDb database
#' @param single.strand.genes.only Logical, default FALSE.
get_grange_genebody <- function(regionGRanges, TxDb, single.strand.genes.only = TRUE) {
  geneBody <- GenomicFeatures::genes(TxDb, single.strand.genes.only = single.strand.genes.only) %>%
    plyranges::join_overlap_intersect(regionGRanges)
  return(geneBody)
}


#' Default ggplot theme for counts plot
.counts_plot_default_theme <- list(
  panel.grid = ggplot2::element_blank(),
  plot.margin = grid::unit(c(0, 0, 0, 0), "cm"),
  legend.title = ggplot2::element_text(size = 12),
  legend.text = ggplot2::element_text(size = 10),
  axis.text.y = ggplot2::element_blank(),
  axis.ticks.y = ggplot2::element_blank(),
  panel.border = ggplot2::element_blank(),
  strip.background = ggplot2::element_blank(),
  strip.text.x = ggplot2::element_text(size = 10, angle = 0),
  strip.text.y = ggplot2::element_text(size = 10, angle = 0)
)


#' Plot normalized counts for each sample
#'
#' Used in `plotRegion()` for the counts tracks
#'
#' @param countdf  A dataframe that comes from `getbpCounts()` or `getbpInserts()`
#' @param plotType Options include 'overlaid','area', 'line', or 'RidgePlot'. default is 'area', which will plot a seperate track for each group.
#' Setting plotType to 'overlaid' will overlay count plot histograms across samples, instead of faceting out separately.
#' Setting plotType to 'RidgePlot' will generate a ridgeplot across all groups.
#' @param base_size Numeric, default 12. Global plot base text size parameter
#' @param counts_color Optional color palette. A named vector of color values where names are unique values in the `color_var` column
#' @param range_label_size Numeric value, default 4. Text size for the y-axis range label
#' @param legend.position Any acceptable `legend.position` argument to theme(). Default NULL will place legend for overlaid plots at (0.8,0.8),
#' or to the "right" for faceted plots.
#' @param facet_label_side Direction character value, default "top". Can also be "right", "left", or "bottom". Position of facet label.
#' @param counts_group_colors Optional named color vector. Values as colors, names are levels of `counts_color_var`. If provided, will color the plots specifically
#' using `scale_color_manual()`
#' @param counts_color_var Character value, default "Groups". Column name from countdf to use to color counts plots.
#' Only used if counts_group_colors provided
#' @param theme_ls A list of named theme arguments passed to theme(), defaults to `.counts_plot_default_theme`. For example, `list(axis.ticks = ggplot2::element_blank())`
#' @return A ggplot object of count histograms by sample.

counts_plot_samples <- function(countdf,
                                plotType = "area",
                                base_size = 12,
                                counts_color = NULL,
                                range_label_size = 4,
                                legend.position = NULL,
                                facet_label_side = "top",
                                counts_group_colors = NULL,
                                counts_color_var = "Groups",
                                theme_ls = .counts_plot_default_theme) {
  assertthat::assert_that("Counts" %in% names(countdf))
  assertthat::assert_that("Locus" %in% names(countdf))
  assertthat::assert_that("Groups" %in% names(countdf))
  assertthat::assert_that(counts_color_var %in% names(countdf))

  # Fill in theme any unspecified theme options with defaults
  default_theme <- .counts_plot_default_theme
  unspec_param <- setdiff(names(default_theme), names(theme_ls))
  if (length(unspec_param) > 0) {
    theme_ls <- c(theme_ls, default_theme[unspec_param])
  }

  # Range df for plotting
  df_range <- data.frame(
    x = max(countdf$Locus),
    y = max(countdf$Counts),
    label = paste0("Range:", 0, "-", round(max(countdf$Counts), digits = 2))
  )

  # Intialize plot
  p1 <- ggplot(data = countdf, aes(x = Locus, y = Counts)) +
    theme_bw(base_size = base_size)


  # Plots
  if (tolower(plotType) == "overlaid") {
    # Conditional Theme
    if (is.null(legend.position)) {
      legend.position <- c(0.8, 0.8)
    }
    theme_ls <- c(
      theme_ls,
      list(
        legend.position = legend.position
      )
    )

    # Base Plot
    p1 <- p1 +
      geom_line(aes(color = !!as.name(counts_color_var)), alpha = 0.75, size = 1.5) +
      ylab(NULL) +
      labs(Groups = "Groups") +
      coord_cartesian(clip = "off") +
      geom_text(
        data = df_range,
        aes(x = x, y = y, label = label),
        size = range_label_size,
        hjust = 1,
        vjust = 1
      ) +
      do.call(theme, theme_ls)

    # Conditional Plot Elements
    if (!is.null(counts_group_colors)) {
      # assertthat::assert_that(all(unique(countdf[[counts_color_var]]) %in% names(counts_group_colors)),
      #                        msg = "Must supply colors for all levels of color variable")
      p1 <- p1 + scale_color_manual(values = counts_group_colors, breaks = names(counts_group_colors))
    }
  } else if (tolower(plotType) == "area") {
    # Conditional Theme
    if (is.null(legend.position)) {
      legend.position <- "right"
    }
    theme_ls <- c(
      theme_ls,
      list(
        legend.position = legend.position
      )
    )

    # Base Plot, conditional
    if (is.null(counts_color)) {
      p1 <- p1 +
        geom_area(aes(fill = !!as.name(counts_color_var)), position = "identity")
    } else {
      p1 <- p1 +
        geom_area(fill = counts_color, position = "identity")
    }

    # Base Plot, common elements
    p1 <- p1 +
      ylab(NULL) +
      facet_wrap(vars(Groups), ncol = 1, strip.position = facet_label_side) +
      geom_text(
        data = df_range,
        aes(x = x, y = y, label = label),
        size = range_label_size,
        hjust = 1,
        vjust = 1
      ) +
      do.call(theme, theme_ls)

    # Conditional Plot Elements
    if (!is.null(counts_group_colors)) {
      # commenting this section out so we can supply specific colors to highlight just a subset if we want--other samples will just be gray
      # assertthat::assert_that(all(unique(countdf[[counts_color_var]]) %in% names(counts_group_colors)),
      #                        msg = "Must supply colors for all levels of color variable")
      p1 <- p1 +
        scale_fill_manual(values = counts_group_colors, breaks = names(counts_group_colors))
    }
  } else if (tolower(plotType) == "line") {
    # Conditional Theme
    if (is.null(legend.position)) {
      legend.position <- "right"
    }
    theme_ls <- c(
      theme_ls,
      list(
        legend.position = legend.position
      )
    )

    # Base Plot, conditional
    if (is.null(counts_color)) {
      p1 <- p1 +
        geom_line(aes(color = !!as.name(counts_color_var)), position = "identity")
    } else {
      p1 <- p1 +
        geom_line(color = counts_color, position = "identity")
    }

    # Base Plot, common elements
    p1 <- p1 +
      ylab(NULL) +
      facet_wrap(vars(Groups), ncol = 1, strip.position = facet_label_side) +
      geom_text(
        data = df_range,
        aes(x = x, y = y, label = label),
        size = range_label_size,
        hjust = 1,
        vjust = 1
      ) +
      do.call(theme, theme_ls)

    # Conditional Plot Elements
    if (!is.null(counts_group_colors)) {
      # commenting this section out so we can supply specific colors to highlight just a subset if we want--other samples will just be gray
      # assertthat::assert_that(all(unique(countdf[[counts_color_var]]) %in% names(counts_group_colors)),
      #                        msg = "Must supply colors for all levels of color variable")
      p1 <- p1 +
        scale_fill_manual(values = counts_group_colors, breaks = names(counts_group_colors))
    }
  } else if (tolower(plotType) == "ridgeplot") {


    ## Conditional elements needed to make RidgePlots work.
    countdf_tmp <- as.data.frame(table(countdf$Groups))
    countdf$Groups2 <- match(countdf$Groups, countdf_tmp$Var1)


    # Base plot, conditional
    p1 <- p1 +
      geom_ridgeline(
        data = as.data.frame(countdf),
        aes(
          x = Locus, y = Groups2, height = Counts,
          fill = Groups
        ),
        alpha = 0.25
      ) +
      ylab(NULL) +
      scale_y_continuous(
        breaks = c(1:length(tmp$Var1)),
        label = tmp$Var1
      ) +
      geom_text(
        data = data.frame(
          x = Inf, y = Inf,
          label = paste("Range:", 0, "-", round(max(countdf$Counts), digits = 2), sep = "")
        ),
        aes(x = x, y = y, label = label),
        size = 2, hjust = 0.9, vjust = 1.4
      ) +
      theme(legend.position = "none")
  } else {
    stop("Error: Plot type not recognized. Please check input for variable 'plotType'")
  }

  return(p1)
}

#' Get scaled values for custom breaks
#'
#' Get scaled values for custom breaks in a given value vector for use
#' defining breaks in scale_gradientn(), for example
#' @param vector of breaks
#' @param x vector of weights
breaks_to_scaledbreaks <- function(breaks, x) {
  rescaled_weights <- scales::rescale(x)
  rescaled_breaks <- quantile(rescaled_weights, probs = ecdf(x)(breaks))
  return(rescaled_breaks)
}

cleanup_breaks <- function(breaks, x) {
  breaks <- sort(breaks)
  n_lower <- sum(breaks < min(x))
  while (n_lower > 1) {
    breaks <- breaks[-1]
    n_lower <- sum(breaks < min(x))
  }
  n_higher <- sum(breaks > max(x))
  while (n_higher > 1) {
    breaks <- breaks[-length(breaks)]
    n_higher <- sum(breaks > max(x))
  }
  return(breaks)
}


#' Get GRanges of motif annotations within a region
#'
#' Helper to return the motif annotations within a specfic region outside of
#' the plotting function. Can supply either a counts df or a region string.
#'
#' @param ArchRProj The archr project containing motif annotations
#' @param motifSetName The name of the motif annotations in the ArchR project
#' @param countdf A counts data frame object from getPop
get_motifs_in_region <- function(motifsList, countdf = NULL, regionString = NULL, numCores = 1, metaColumn = NULL, cellSubsets = "ALL") {
  chrom <- toString(unique(countdf$chr))
  startSite <- min(countdf$Locus)
  endSite <- max(countdf$Locus)
  regionGRanges <- GRanges(seqnames = chrom, ranges = IRanges(start = startSite, end = endSite), strand = "*")

  specMotifs <- unlist(motifsList) %>%
    plyranges::mutate(name = gsub("_.*", "", names(.))) %>%
    plyranges::join_overlap_intersect(regionGRanges) %>%
    plyranges::mutate(type = "exon") %>%
    plyranges::mutate(index = seq(1, length(.), by = 1)) %>%
    plyranges::mutate(labels = paste(name, index, sep = "_"))

  return(specMotifs)
}

#' Plot Motifs on Counts Plot
#'
#' Add motif annotation labels on an existing counts plot. If providing weights for motif
#' text colors, ensure the input counts plot is formatted accordingly (ie counts_color = "gray90")
#' for best visibility.
#'
#' @param p1 The output of `counts_plot_samples()`
#' @param countdf  A dataframe that comes from `getbpCounts()` or `getbpInserts()`
#' @param ArchRProj An archR project containing motif annotations
#' @param motifSetName The name of the motif set in ArchRProj to use for annotation
#' @param motif_y_space_factor A factor for vertical spacing between motif labels. Default 4. Increase to make labels farther apart, decrease to make labels closer.
#' @param motif_stagger_labels_y = FALSE Logical value, default FALSE. If TRUE, will  stagger motif labels in adjacent columns in the vertical direction
#' @param motif_weights Optional numeric vector, default NULL. If provided will be used to color motif labels by the weighted values. Values must be uniquely named
#' with motif names, for example c(`KLF5`= 3.2, `STAT1 = 0.2`, `EOMES` = -1.4`). Weights can be anything relevant, for example if the peak/region is associated with
#' a specific group/sample then global motif enrichment results for that group: `-log10(FDR)*sign(change)`
#' @param motif_weight_name Character value, default "Motif Weight". Used to label the color legend.
#' @param weight_colors Named numeric vector. Names should be color values and breaks should be the corresponding values of motif_weights. Values outside the
#' highest and lowest value will appear as max or min defined color value.
#' @param motif_lab_size Numeric value, default 1. Size of motif labels.
#' @param motif_lab_alpha Numeric value, default 0.25. Alpha for motif labels.
#' @param motif_line_size Numeric value, default 1. Size of motif lines.
#' @param motif_line_alpha Numeric value, default 0.25. Alpha for motif lines.
#' @return The input ggplot object with motif labels overlaid
counts_plot_motif_overlay <- function(p1,
                                      countdf,
                                      motifsList,
                                      motif_y_space_factor = 4,
                                      motif_stagger_labels_y = FALSE,
                                      motif_weights = NULL,
                                      motif_weight_name = "Motif Weight",
                                      motif_weight_colors = c(darkblue = -10, gray = 0, darkred = 10),
                                      motif_lab_size = 1,
                                      motif_lab_alpha = 0.25,
                                      motif_line_size = 0.75,
                                      motif_line_alpha = 0.25) {

  # Retrieve annotations in region and format
  specMotifs <- get_motifs_in_region(
    countdf = countdf,
    motifsList = motifsList,
  )
  reduceMotifs <- plyranges::reduce_ranges(specMotifs, count = plyranges::n(), names = paste(labels, collapse = ","))
  splitMotifs <- strsplit(reduceMotifs$names, split = ",")

  # Calculate y spacing for motif labels
  minSize <- max(countdf$Counts) / (length(specMotifs) / motif_y_space_factor)
  y1 <- rep(0, length(specMotifs))
  names(y1) <- specMotifs$labels

  for (i in 1:length(splitMotifs)) {
    if (i %% 2 == 0 & motif_stagger_labels_y) { # stagger y spacing in neighboring columns
      start_val <- minSize / 2
    } else {
      start_val <- 0
    }

    y1[splitMotifs[[i]]] <- seq(start_val, start_val + minSize * (length(splitMotifs[[i]]) - 1), by = minSize)
  }

  # TF label coordinates and labels
  tmp_motifdf <- data.frame(
    x1 = start(specMotifs),
    x2 = start(specMotifs) + width(specMotifs),
    y = y1,
    name = specMotifs$labels
  ) %>%
    pivot_longer(cols = c("x1", "x2"), names_to = NULL, values_to = "x") %>%
    group_by(name) %>%
    dplyr::mutate(labels = ifelse(max(x) == x, gsub("_.*", "", name), NA))

  # Incorprate Weights
  if (!is.null(motif_weights)) {
    assertthat::assert_that(is.numeric(motif_weight_colors))

    if (length(intersect(names(motif_weights), tmp_motifdf$labels)) == 0) {
      warning(sprintf(
        "None of the supplied motif weight names match expected motif labels. Example motif label format: %s",
        paste(head(unique(na.omit(tmp_motifdf$label)), 3), collapse = ", ")
      ))
    }
    tmp_motifdf <- tmp_motifdf %>%
      dplyr::mutate(mweight = ifelse(labels %in% names(motif_weights), motif_weights[labels], NA)) %>%
      dplyr::mutate(mweight = ifelse(!all(is.na(mweight)), max(mweight, na.rm = T), NA)) %>% # still grouped by name. this ensures weights are applied to each row of same motif
      dplyr::mutate(mweight = case_when( # clip ends within color range
        mweight > max(motif_weight_colors) ~ max(motif_weight_colors),
        mweight < min(motif_weight_colors) ~ min(motif_weight_colors),
        is.na(mweight) ~ as.numeric(NA),
        TRUE ~ mweight
      ))
  } else {
    tmp_motifdf$mweight <- 1
  }

  # Plot
  p1 <- p1 +
    geom_line(
      data = tmp_motifdf,
      aes(x = x, y = y, group = name, color = mweight),
      alpha = motif_line_alpha,
      size = motif_line_size
    ) +
    # vertical labels
    ggrepel::geom_text_repel(
      data = tmp_motifdf[tmp_motifdf$y > 0 & !is.na(tmp_motifdf$labels), ], # removed the NA rows to prevent warnings for intentionally missing labels
      aes(x = x, y = y, label = labels, color = mweight),
      direction = "x",
      arrow = arrow(length = grid::unit(0.5, "mm")),
      alpha = motif_lab_alpha,
      size = motif_lab_size,
      segment.size = 0.25,
      max.overlaps = 50,
      min.segment.length = 0
    ) +
    # motif labels along x-axis
    ggrepel::geom_text_repel(
      data = tmp_motifdf[tmp_motifdf$y == 0 & !is.na(tmp_motifdf$labels), ], # removed the NA rows to prevent warnings for intentionally missing labels
      aes(x = x, y = y, label = labels, color = mweight),
      arrow = arrow(length = grid::unit(0.5, "mm")),
      alpha = motif_lab_alpha,
      size = motif_lab_size,
      direction = "y",
      segment.size = 0.25,
      min.segment.length = 0,
      nudge_y = -max(countdf$Counts) / 20
    ) +
    ylim(-max(countdf$Counts) / 10, max(countdf$Counts)) +
    xlim(min(countdf$Locus), max(countdf$Locus))

  # Color text by motif weights
  if (!is.null(motif_weights)) {
    scaled_breaks <- breaks_to_scaledbreaks(motif_weight_colors, na.omit(tmp_motifdf$mweight))

    p1 <- p1 + scale_color_gradientn(
      colors = names(motif_weight_colors), values = scaled_breaks,
      na.value = "gray", # may want to parameterize this
      name = motif_weight_name
    )
  } else {
    p1 <- p1 + scale_color_gradient(low = "black", high = "black", na.value = "gray", name = motif_weight_name, guide = "none")
  }

  return(p1)
}


#' Plot colors in a color palette
#'
#' Utility function to quickly visualize a color palette formatted for motif weight visualization
#'
#' @param pal A palette, named list of values where names are colors and values are the color breaks
#' @return Plots a base plot to visualize the input palette
#' @examples
#' motif_color_pal <- c(-20, -10, 0, 10, 20)
#' names(motif_color_pal) <- c("darkblue", "cyan", "gray", "orange", "darkred")
#' plot_pal(motif_color_pal)
plot_pal <- function(pal, cex = 4, cex.lab = 1) {
  npal <- length(pal)
  plotdim(npal * 1, 3)
  plot(1:npal, rep(1, npal), pch = 19, cex = cex, col = names(pal))
  text(1:npal, rep(1, npal), labels = pal, cex = cex.lab)
}

#' Get Gene Body Model
#'
#' Get Gene Body model for specific gene in region ranges
#'
#' Used in `plotRegion()` to preformat GRanges for selected genes
#'
#' @param whichGene Character value. Gene symbol for a single gene to plot.
#' @param regionGRanges Granges object for region of interest. Ie, an output of `countdf_to_region()`
#' @param countdf  A dataframe that comes from `getbpCounts()` or `getbpInserts()`
#' @param orgdb An organism database containing the gene
#' @param db_id_col Character value. Column in `orgdb` containing the output id. Default "REFSEQ".
#' @param verbose Logical value, default TRUE.
#' @return GRanges object for gene to plot
get_gene_body_model <- function(whichGene,
                                regionGRanges,
                                countdf,
                                orgdb,
                                db_id_col = "REFSEQ",
                                verbose = TRUE) {
  # Check for dependency RMariaDB needed for GenomicFeatures::makeTxDbFromUCSC
  if (!requireNamespace("RMariaDB", quietly=TRUE)) {
    stop("Couldn't load the package RMariaDB. RMariaDB must be installed to plot specific genes.")
  }
  
  # Get REFSEQ values for gene symbol
  txList <- tryCatch(
    unlist(AnnotationDbi::mapIds(
      x = orgdb,
      keys = whichGene,
      column = db_id_col,
      keytype = "SYMBOL",
      multiVals = "list"
    )),
    error = function(e) {
      stop(sprintf(
        "Error attempting to query symbol for 'whichGene' argument [%s] from specified database. Full error message: %s",
        paste(whichGene, collapse = ","), e
      ))
      return(NULL)
    }
  )

  # Pull specified database and table for the transcripts
  # This requires package "RMariaDB" available on CRAN and corresponding linux library "libmariadb-dev"
  newTxDB <- GenomicFeatures::makeTxDbFromUCSC(genome = "hg38", tablename = "refGene", transcript_ids = txList)
  if (verbose) {
    print(txList)
  }

  # Get granges object from database and intersect with region granges
  newTxDBgenes <- unlist(genes(newTxDB, single.strand.genes.only = FALSE)) %>%
    plyranges::filter(!grepl("_", seqnames))
  newGRanges <- biovizBase::crunch(newTxDB, which = newTxDBgenes) %>%
    plyranges::join_overlap_intersect(regionGRanges)
  colnames(values(newGRanges))[4] <- "model" # rename 'Type' column to models

  if (length(newGRanges) > 0) {
    newGRangesf <- as.data.frame(newGRanges)
    startSide <- newGRangesf[newGRangesf$start %in% start(regionGRanges) & newGRangesf$model == "gap", ]
    endSide <- newGRangesf[newGRangesf$end %in% end(regionGRanges) & newGRangesf$model == "gap", ]
    beginexon <- NULL
    endexon <- NULL

    if (dim(startSide)[1] > 0) {
      beginexon <- data.frame(
        seqnames = unique(countdf$chr),
        start = min(countdf$Locus) - 1,
        end = min(countdf$Locus) - 1,
        width = 1,
        strand = "*",
        tx_id = unique(startSide$tx_id),
        tx_name = unique(startSide$tx_name),
        gene_id = unique(startSide$gene_id),
        model = "exon"
      )
    }

    if (dim(endSide)[1] > 0) {
      endexon <- data.frame(
        seqnames = unique(countdf$chr),
        start = max(countdf$Locus) + 1,
        end = max(countdf$Locus) + 1,
        width = 1,
        strand = "*",
        tx_id = unique(endSide$tx_id),
        tx_name = unique(endSide$tx_name),
        gene_id = unique(endSide$gene_id),
        model = "exon"
      )
    }

    newModel <- rbind(newGRangesf, beginexon, endexon) %>%
      plyranges::mutate(nameList = paste(whichGene, tx_name, sep = ": ")) %>%
      makeGRangesListFromDataFrame(., split.field = "nameList", keep.extra.columns = TRUE)

    # check for transcripts that are identical within the window
    uniqueTranscripts <- stack(newModel)[!duplicated(ranges(stack(newModel)))] %>%
      plyranges::filter(model != "utr") %>%
      as.data.frame() %>%
      dplyr::select(tx_name) %>%
      unique(.)
    newModel <- newModel[grepl(paste(unlist(uniqueTranscripts), collapse = "|"), names(newModel))]

    names(newModel) <- rep(whichGene, length(newModel))

    return(newModel)
  } else {
    print("Warning: whichGene did not contain genes names for any genes in the window. Plotting all transcripts")
    return(NULL)
  }
}

#' Common theme for gene plots
.gene_plot_theme <- list(
  panel.grid = ggplot2::element_blank(),
  panel.border = ggplot2::element_blank(),
  axis.ticks = ggplot2::element_blank(),
  axis.text.y = ggplot2::element_blank(),
  plot.margin = grid::unit(c(.5, .5, .5, .5), "cm")
)


#' Generate a ggplot based on a Gene Model
#'
#' Creates a gene track ggplot based on gene model GRanges
#' Used conditinally in `plotRegion()` for gene plot.
#'
#' @param A gene granges model, ie generated by `get_gene_body_model()`
#' @param x_min Numeric value. Min x-axis limit. If NULL (default) will use min value from newModel granges.
#' @param x_max Numeric value. Max x-axis limit. If NULL (default) will use max value from newModel granges.
#' @param base_size Numeric, default 12. Global plot base text size parameter
#' @param theme_ls Named list of parameters passed to `theme()`. For defaults see `.gene_plot_theme`
plot_whichGene <- function(newModel,
                           x_lim = NULL,
                           base_size = 12,
                           theme_ls = .gene_plot_theme) {
  # if(is.null(x_min)){
  #     xmin = min(newModel$
  # }

  # Fill in theme any unspecified theme options with defaults
  default_theme <- .gene_plot_theme
  unspec_param <- setdiff(names(default_theme), names(theme_ls))
  if (length(unspec_param) > 0) {
    theme_ls <- c(theme_ls, default_theme[unspec_param])
  }

  p_gene <- ggbio::ggbio() +
    ggbio::geom_alignment(newModel,
      aes(type = model),
      cds.rect.h = 0.25,
      rect.height = 0.25 / 4
    )

  if (!is.null(x_lim)) {
    p_gene <- p_gene + xlim(x_lim[1], x_lim[2])
  }

  p_gene <- p_gene +
    scale_y_continuous(expand = c(0.3, 0.3)) +
    theme_bw(base_size = base_size) +
    coord_cartesian(clip = "off") +
    do.call(theme, theme_ls)

  return(p_gene)
}

#' Plot a genebody from a specific database
#'
#' Used in `plotRegion()` to plot gene track for granges input containing all genes in a region
#'
#' @param organismdb Database object, for example the output of `OrganismDbi::makeOrganismDbFromTxDb()`
#' @param geneBody_gr A GRanges object containing the gene, for example an output of `get_grange_genebody()`
#' @param xlim vector of (x_min, x_max)
#' @param base_size Numeric, default 12. Global plot base text size parameter
#' @param plot_theme Named list of `ggplot2::theme()` parameters.
#' @param collapseGenes Logical value, default FALSE. If TRUE will collapse all genes into one line.
#' sites into a single row
#' @return A ggplot object
plot_geneBody <- function(organismdb,
                          geneBody_gr,
                          x_lim = NULL,
                          base_size = 12,
                          theme_ls = .gene_plot_theme,
                          collapseGenes = FALSE) {
  # Fill in theme any unspecified theme options with defaults
  default_theme <- .gene_plot_theme
  unspec_param <- setdiff(names(default_theme), names(theme_ls))
  if (length(unspec_param) > 0) {
    theme_ls <- c(theme_ls, default_theme[unspec_param])
  }

  if (collapseGenes) {
    stat_val <- "reduce"
  } else {
    stat_val <- NULL
  }

  g_geneBody <- ggbio::autoplot(organismdb, wh = geneBody_gr) # , stat = stat_val)

  if (!is.null(x_lim)) { # had to move this here, otherwise theme would be overwritten...
    g_geneBody <- g_geneBody + xlim(x_lim[1], x_lim[2])
  }

  g_geneBody <- g_geneBody +
    scale_y_continuous(expand = c(0.3, 0.3)) +
    theme_bw(base_size = base_size) +
    coord_cartesian(clip = "off") +
    do.call(theme, theme_ls)


  return(g_geneBody)
}

##### simplifiedOrgDb: Generates a simplified OrganismDb object for plotting only the longest transcript for each gene.

#' Generate a simplified Organism Db object.
#'
#' Pulls out only the longest transcripts for each gene and repackages them into an organism database object for plotting.
#'
#' @param TxDb A TxDb object with transcript info.
#' @param orgdb A OrgDb object with gene name info.

simplifiedOrgDb <- function(TxDb = TxDb, orgdb = orgdb) {

  # retrieve transcript lengths
  txlen <- transcriptLengths(TxDb, with.utr5_len = TRUE, with.utr3_len = TRUE)
  setDT(txlen)
  txlen$len <- rowSums(as.matrix(txlen[, .(tx_len, utr5_len, utr3_len)]))

  setkey(txlen, gene_id, len, tx_id)

  # filter longesttranscript by gene_id
  ltx <- txlen[!is.na(gene_id)][, tail(.SD, 1), by = gene_id]$tx_id

  # filter txdb object
  txb <- as.list(TxDb)
  txb$transcripts <- txb$transcripts[txb$transcripts$tx_id %in% ltx, ]
  txb$splicings <- txb$splicings[txb$splicings$tx_id %in% ltx, ]
  txb$genes <- txb$genes[txb$genes$tx_id %in% ltx, ]
  chrominfo <- as.data.frame(seqinfo(TxDb)) %>%
    dplyr::mutate(chrom = rownames(.)) %>%
    dplyr::rename(length = seqlengths, is_circular = isCircular)
  txb2 <- GenomicFeatures::makeTxDb(txb$transcripts, txb$splicings, txb$genes, chrominfo,
    metadata =
      filter(
        metadata(TxDb),
        grepl("Genome|Organism|Data source|UCSC Table", name)
      )
  )
  Homo.sapiens.hg38 <- OrganismDbi::makeOrganismDbFromTxDb(txb2,
    orgdb = orgdb
  )

  return(Homo.sapiens.hg38)
}




get_link_plot <- function(regionGRanges, legend.position = NULL,
                          relativeHeights, linkdf) {
  linkdf2 <- dplyr::filter(linkdf, start + 250 > start(regionGRanges) & end - 250 < end(regionGRanges))

  ## Set Curvature to fit window
  if (is.null(legend.position) & relativeHeights["Links"] / 3 <= 0.5) {
    curveVal <- relativeHeights["Links"] / 3 - 0.05
  } else if (is.null(legend.position)) {
    curveVal <- 0.5
  } else if (relativeHeights["Links"] / 3 <= 0.5) {
    curveVal <- relativeHeights["Links"] / 3 - 0.2
  } else {
    curveVal <- 0.4
  }

  p5 <- ggplot() +
    geom_curve(aes(
      x = start + 250, xend = end - 250,
      y = y, yend = y, color = Correlation
    ),
    curvature = curveVal,
    data = cbind(linkdf2, y = rep(0, dim(linkdf2)[1]))
    ) +
    theme_minimal() +
    ylab(NULL) +
    xlab(NULL) +
    ylim(-1, 0) +
    scale_colour_viridis_c(breaks = c(
      ceiling(10 * min(linkdf2$Correlation)) / 10,
      0,
      floor(10 * max(linkdf2$Correlation)) / 10
    )) +
    coord_cartesian(clip = "off", ylim = c(-0.75, 0)) +
    theme(
      panel.grid = ggplot2::element_blank(), panel.border = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()
    )

  return(p5)
}
