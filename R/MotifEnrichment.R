#' @title \code{EnrichedRanges}
#'
#' @description Test for enrichment of motifs within Group1 against a background
#'   Group2.
#'
#' @param Group1 A GRanges object, such as a set of significant differential
#'   tiles.
#' @param Group2 A GRanges object containing background regions, non-overlapping
#'   with Group1
#' @param Category A GRanges representing a motif or other category of positions
#' @param type Optional, name of a metadata column in Group1 and Group2 to test
#'   for enrichment the number of unique entries in column given by 'type'.
#'   Default is NULL, which tests the number of number of Ranges.
#' @param returnTable Optional parameter. Will return the table used for the hypergeometric enrichment test instead of the p-value. 
#'
#' @return A data.frame of enrichment
#'
#' @noRd
#'
EnrichedRanges <- function(Group1,
                           Group2,
                           Category,
                           type = NULL,
                           returnTable = FALSE) {
  Group1Cat <- plyranges::filter_by_overlaps(Group1, Category)
  Group2Cat <- plyranges::filter_by_overlaps(Group2, Category)

  OnlyGroup1 <- plyranges::filter_by_non_overlaps(Group1, Category)
  OnlyGroup2 <- plyranges::filter_by_non_overlaps(Group2, Category)

  if (returnTable && is.null(type)) {
    dt_table <- data.frame(
      Group1 = c(length(Group1Cat), length(OnlyGroup1)),
      Group2 = c(length(Group2Cat), length(OnlyGroup2)),
      row.names = c("In Category", "Not in Category")
    )
    return(t(dt_table))
  } else if (returnTable && sum(c(
    colnames(GenomicRanges::mcols(Group1)),
    colnames(GenomicRanges::mcols(Group2))
  ) %in% type) == 2 && length(type) == 1) {
    dt_table <- data.frame(
      Group1 = c(
        length(unique(GenomicRanges::mcols(Group1Cat)[, type])),
        length(unique(GenomicRanges::mcols(OnlyGroup1)[, type]))
      ),
      Group2 = c(
        length(unique(GenomicRanges::mcols(Group2Cat)[, type])),
        length(unique(GenomicRanges::mcols(OnlyGroup2)[, type]))
      ),
      row.names = c("In Category", "Not in Category")
    )
    return(t(dt_table))
  } else if (returnTable) {
    stop("Incorrect method or column name. Please check input")
  }

  if (is.null(type)) {
    pVal <- phyper(
      q = length(Group1Cat),
      m = length(Group1),
      n = length(Group2),
      k = length(Group1Cat) + length(Group2Cat),
      lower.tail = FALSE
    )
    enrichment <- (
      (length(Group1Cat) / length(Group1)) /
        (length(Group2Cat) / length(Group2))
    )
  } else if (sum(c(
    colnames(GenomicRanges::mcols(Group1)),
    colnames(GenomicRanges::mcols(Group2))
  ) %in% type) == 2 && length(type) == 1) {
    pVal <- phyper(
      q = length(unique(GenomicRanges::mcols(Group1Cat)[, type])),
      m = length(unique(GenomicRanges::mcols(Group1)[, type])),
      n = length(unique(GenomicRanges::mcols(Group2)[, type])),
      k = (length(unique(GenomicRanges::mcols(Group1Cat)[, type]))
      + length(unique(GenomicRanges::mcols(Group2Cat)[, type]))),
      lower.tail = FALSE
    )

    enrichment <- (
      (length(unique(GenomicRanges::mcols(Group1Cat)[, type]))
      / length(unique(GenomicRanges::mcols(Group1)[, type]))) /
        (length(unique(GenomicRanges::mcols(Group2Cat)[, type]))
        / length(unique(GenomicRanges::mcols(Group2)[, type])))
    )
  } else {
    stop("Incorrect method or column name. Please check input")
  }

  return(data.frame(p_value = pVal, enrichment = enrichment))
}

#' @title \code{MotifEnrichment}
#'
#' @description Test for enrichment of motifs within Group1 against a background
#'   Group2.
#'
#' @param Group1 A GRanges object, such as a set of significant differential
#'   tiles.
#' @param Group2 A GRanges object containing background regions, non-overlapping
#'   with Group1
#' @param motifPosList A GRangesList of motifs and positions for each motif.
#'   Must be named for each motif.
#' @param type Optional, name of a metadata column in Group1 and Group2 to test
#'   for enrichment the number of unique entries in column given by 'type'.
#'   Default is NULL, which tests the number of number of Ranges.
#'
#' @return allEnrichmentList
#'
#' @export
#'
MotifEnrichment <- function(Group1, Group2, motifPosList, type = NULL) {
  allEnrichmentList <- pbapply::pblapply(motifPosList, function(x) {
    tmp_df <- EnrichedRanges(Group1, Group2, Category = x, type = type)
  }, cl = 1) # cl = 1 to disable multiprocessing as it is not much faster
  df_final <- do.call("rbind", allEnrichmentList)

  df_final$adjp_val <- p.adjust(df_final$p_value, method = "fdr")
  df_final$mlog10Padj <- -log10(df_final$adjp_val)

  return(df_final)
}

#' @title{Gene2Motif}
#'
#' @description
#'
#' @param TSS_Sites GRanges containing TSS sites of interest. Must include a
#'   column 'name' which has the associated gene name.
#' @param allTiles GRanges containing all tiles
#' @param TSS_Links a data.table object that record all the Tile-Tile links by
#'   co-accessibility. Must include columns named 'Tile1' and 'Tile2' which
#'   contain a string describing each Tile in the format 'chr1:100-2000' and
#'   must be identical to Tiles listed in allTiles
#' @param motifPosList A GRangesList of motifs and positions for each motif.
#'   Must be named for each motif.
#' @param numCores Optional, the number of cores to use with multiprocessing.
#'   Default is 1.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#'
#' @return
#'
#' @export
#'
Gene2Motif <- function(TSS_Sites, allTiles, TSS_Links, motifPosList,
                       numCores = 1, verbose = FALSE) {
  if (verbose) {
    message("Generating TSS-Tile Network.")
  }

  TSS_Network <- c(
    TSS_Links$Tile1, TSS_Links$Tile2,
    MOCHA::GRangesToString(TSS_Sites)
  ) %>%
    unique() %>%
    StringsToGRanges(.) %>%
    plyranges::filter_by_overlaps(allTiles, .)

  if (verbose) {
    message(
      "Finding all motifs related to each Tile within the ",
      "TSS-Tile Network."
    )
  }

  # Let's find all the motifs that overlap with each Tile
  # within the altTSS Network
  cl <- parallel::makeCluster(numCores)
  parallel::clusterExport(cl,
    varlist = c("motifPosList", "TSS_Network"),
    envir = environment()
  )

  tmpOverlap <- pbapply::pblapply(seq_along(motifPosList), function(x) {
    avgWidth <- mean(GenomicRanges::width(motifPosList[[x]]))
    ifelse(
      plyranges::count_overlaps(
        TSS_Network, motifPosList[[x]],
        minoverlap = avgWidth
      ) > 0,
      names(motifPosList)[x],
      NA
    )
  }, cl = cl)

  overlap_df <- do.call("cbind", tmpOverlap)
  colnames(overlap_df) <- names(motifPosList)
  rownames(overlap_df) <- MOCHA::GRangesToString(TSS_Network)

  parallel::clusterExport(cl,
    varlist = c("overlap_df"),
    envir = environment()
  )

  motifList <- pbapply::pblapply(c(seq_len(dim(overlap_df)[1])), function(x) {
    ifelse(any(!is.na(overlap_df[x, ])),
      list(overlap_df[x, which(!is.na(overlap_df[x, ]))]),
      NA
    )
  }, cl = cl)

  if (verbose) {
    message(
      "Finding all Tiles related to each gene within the ",
      "TSS-Tile Network."
    )
  }

  parallel::clusterExport(cl,
    varlist = c("TSS_Sites", "allTiles"),
    envir = environment()
  )

  name <- Tile1 <- Tile2
  Tile2Gene <- pbapply::pblapply(unique(TSS_Sites$name), function(x) {
    filtTSS <- plyranges::filter(TSS_Sites, name == x)
    geneTSS <- MOCHA::GRangesToString(
      plyranges::filter_by_overlaps(allTiles, filtTSS)
    )

    tmp <- TSS_Links[Tile1 %in% geneTSS | Tile2 %in% geneTSS, ]

    if (dim(tmp)[1] > 0) {
      unique(c(tmp$Tile1, tmp$Tile2, geneTSS))
    } else {
      geneTSS
    }
  }, cl = cl)
  names(Tile2Gene) <- unique(TSS_Sites$name)



  if (verbose) {
    message("Linking Motifs to each gene within the TSS-Tile Network")
  }
  ## Link all the genes to motifs via Tile2Gene and the motifList

  parallel::clusterExport(cl,
    varlist = c("TSS_Network", "motifList"),
    envir = environment()
  )

  Gene2Motif <- pbapply::pblapply(Tile2Gene, function(x) {
    # Find which indices of the AltTSS Network GRanges are linked to that Gene
    tmp <- GenomicRanges::findOverlaps(MOCHA::StringsToGRanges(x), TSS_Network)
    # Pull up and unlist all the motifs associated with those tiles.
    unlist(motifList[S4Vectors::subjectHits(tmp)])
  }, cl = cl)

  parallel::stopCluster(cl)

  return(Gene2Motif)
}
