#' @title \code{getAltTSS} Annotate Peaks falling in TSS Sites, and identify
#'   alternatively regulated TSSs for each gene.
#'
#' @description \code{getAltTSS} Pulls out all peaks that fall in TSS sites,
#'   annotates them with the name of gene, and identifies genes that have
#'   evidence for alternatively regulated TSSs, including both type i (only some
#'   of the open TSSs for a gene are significantly more (or less) accessible),
#'   and type ii (multiple TSSs are significant different, with some being more
#'   accessible and others less). Alternatively, this function  will return all
#'   open TSSs with differential measurements if the returnAllTSS flag is set to
#'   TRUE.
#'
#' @param completeDAPs GRanges object that contains the differential
#'   measurements across all peaks (unfiltered DAPs). Will also work with
#'   data.frame or data.table version of a GRanges object. If you want
#'   alternatively regulated TSSs, the object must include a column names 'FDR',
#'   and 'Log2FC_C', which is standard for MOCHA differentials.
#' @param returnAllTSS Flag to return all TSS sites with DAPs measurements,
#'   without filtering for alternative TSS usage. If multiple TSSs fall within
#'   the same tile, then that tile will be repeated for each TSS.
#' @param nuancedTSS True/False flag to determine if alternative TSS genes
#'   should be filtered out if all their differential TSS usage falls within too
#'   small of a range. Default is TRUE
#' @param nuancedTSSGAP Minimum distance betweeen TSSs needed for them to
#'   considered distinctly regulated TSSs. If two TSSs are too close, it is
#'   unclear and highly unlikely that ATAC data can distinguish between them.
#'   Default is 150 bp.
#' @param threshold FDR Threshold for determining significant vs non-significant
#'   changes in accessibility. Following MOCHA's standards, default is 0.2.
#'
#' @return tpeaks A GRanges containing annotated peaks falling in TSS sites.
#'
#' @export
#'
getAltTSS <- function(completeDAPs,
                      returnAllTSS = FALSE,
                      nuancedTSS = TRUE,
                      nuancedTSSGap = 150,
                      threshold = 0.2,
                      TxDb = TxDb.Hsapiens.UCSC.hg38.refGene,
                      Org = org.Hs.eg.db) {
  . <- exactTSS <- name <- FDR <- Log2FC_C <- strand <- seqnames <- NULL

  if (grepl("data.table|data.frame", class(completeDAPs)[1])) {
    DAP_GRanges <- GenomicRanges::makeGRangesFromDataFrame(
      as.data.frame(completeDAPs),
      keep.extra.columns = TRUE
    )
  } else if (class(completeDAPs)[1] == "GRanges") {
    DAP_GRanges <- completeDAPs
  } else {
    stop(
      "completeDAPs object type is not compatible. Please submit a GRanges",
      ", data.table, or SummarizedExperiment object."
    )
  }

  tss1 <- suppressWarnings(ensembldb::transcriptsBy(TxDb, by = ("gene")))

  names(tss1) <- AnnotationDbi::mapIds(Org, names(tss1), "SYMBOL", "ENTREZID")

  allT <- IRanges::stack(tss1) %>%
    GenomicRanges::trim(.) %>%
    GenomicRanges::promoters(., upstream = 0, downstream = 0) %>%
    plyranges::mutate(exactTSS = start(.)) %>%
    plyranges::filter(!duplicated(exactTSS)) %>%
    plyranges::anchor_3p(.) %>%
    plyranges::stretch(., extend = 125) %>%
    GenomicRanges::trim()

  tpeaks <- plyranges::join_overlap_intersect(allT, DAP_GRanges)

  if (returnAllTSS) {
    return(tpeaks)
  }

  if (!all(
    colnames(GenomicRanges::mcols(DAP_GRanges)) %in% c("FDR", "Log2FC_C")
  )) {
    stop(
      "GRanges input does not include columns for FDR and Log2FC_C, so ",
      "alternatively regulated TSSs cannot be identified. Please provide an",
      " object with those columns for filtering."
    )
  }

  # We want to select all the genes that have multiple TSSs (duplicated(name))
  # and at least one TSSs that has an FDR <= to the threshold. After that, we
  # want to filter for sites where only a subset of open TSSs change, or all
  # TSSs change, but it opposite directions.
  # may need to consider removing the exactTSS sites
  altTSS <- tpeaks %>%
    plyranges::filter(!duplicated(exactTSS)) %>%
    plyranges::group_by(name) %>%
    plyranges::filter(any(FDR <= threshold) & any(duplicated(name)))


  altTSS <- altTSS %>%
    plyranges::filter(
      ifelse(all(FDR <= threshold, na.rm = TRUE) & !any(is.na(FDR)),
        !all(Log2FC_C > 0) & !all(Log2FC_C < 0), TRUE
      )
    ) %>%
    plyranges::ungroup() %>%
    sort()

  if (nuancedTSS) {
    nuancedGenes <- split(altTSS, as.character(altTSS$name)) %>%
      parallel::mclapply(., function(x) {
        tmp <- IRanges::gaps(x) %>%
          filter(
            seqnames == GenomeInfoDb::seqnames(x) &
              strand == GenomicRanges::strand(x)
          ) %>%
          IRanges::width(.)
        ifelse(length(tmp) == 3 & tmp[2] < nuancedTSSGap, FALSE, TRUE)
      }) %>%
      unlist()

    altTSS <- altTSS %>% filter(name %in% names(nuancedGenes)[nuancedGenes])
  }
  return(altTSS)
}
