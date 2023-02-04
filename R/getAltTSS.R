#' @title \code{getStartSites} Annotate Peaks falling in TSS Sites
#'
#' @description \code{getStartSites} Pulls out all peaks that fall in TSS sites and annotates them with the name of gene
#'
#'
#' @param peakSet an ArchR Project, or a GRangesList of fragments
#' @param TxDb Transcript Db object for organism. Default is Hg38.
#' @param Org Organism Db for annotating gene names across databases. Must match organism of TxDb. Default is org.Hs.eg.db.
#'
#' @return tpeaks A GRanges containing annotated peaks falling in TSS sites.
#'
#' @export
#'
getStartSites <- function(peakSet,
                          TxDb = TxDb.Hsapiens.UCSC.hg38.refGene,
                          Org = org.Hs.eg.db) {
  . <- NULL
  if (grepl("data.table|SummarizedExperiment", class(peakSet)[1])) {
    peaksGRanges <- ExtractGR(peakSet)
  } else if (class(peakSet)[1] == "GRanges") {
    peaksGRanges <- peakSet
  }

  tss1 <- suppressWarnings(ensembldb::transcriptsBy(TxDb, by = ("gene")))

  names(tss1) <- suppressWarnings(mapIds(Org, names(tss1), "SYMBOL", "ENTREZID"))

  allT <- stack(tss1) %>%
    GenomicRanges::trim(.) %>%
    GenomicRanges::promoters(., upstream = 0, downstream = 0) %>%
    plyranges::mutate(exactTSS = start(.)) %>%
    plyranges::filter(!duplicated(exactTSS)) %>%
    GenomicRanges::trim()

  tpeaks <- plyranges::join_overlap_intersect(allT, peaksGRanges)

  return(tpeaks)
}


#' @title \code{getAltTSS} Annotate Peaks falling in TSS Sites
#'
#' @description \code{getAltTSS} Pulls out all peaks that fall in TSS sites and annotates them with the name of gene
#'
#' @param completeDAPs GRanges object that contains the differential measurements across all peaks (unfiltered DAPs)
#' @param returnAllTSS Flag to return all TSS sites with DAPs measurements, without filtering for alternative TSS usage
#' @param nuancedTSS True/False flag to determine if alternative TSS genes should be filtered out if all their differential TSS usage falls within too small of a range. Default is TRUE
#' @param threshold Basepair distance threshold for filtering out genes whose differential TSS sites all falls very close together. Default is 0.2.
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
  . <- exactTSS <- NULL
  if (grepl("data.table|SummarizedExperiment", class(completeDAPs)[1])) {
    DAP_GRanges <- ExtractGR(completeDAPs)
  } else if (class(completeDAPs)[1] == "GRanges") {
    DAP_GRanges <- completeDAPs
  } else {
    stop("completeDAPs object type is not compatible. Please submit a GRanges, data.table, or SummarizedExperiment object.")
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

  ## We want to select all the genes that have multiple TSSs (duplicated(name)) and at least one TSSs that has an FDR <= to the threshold.
  ## After that, we want to filter for sites where only a subset of open TSSs change, or all TSSs change, but it opposite directions.
  ###### may need to consider removing the exactTSS sites
  altTSS <- tpeaks %>%
    plyranges::filter(!duplicated(exactTSS)) %>%
    plyranges::group_by(name) %>%
    plyranges::filter(any(FDR <= threshold) & any(duplicated(name)))


  altTSS <- altTSS %>%
    plyranges::filter(ifelse(all(FDR <= threshold, na.rm = TRUE) & !any(is.na(FDR)),
      !all(Log2FC_C > 0) & !all(Log2FC_C < 0), TRUE
    )) %>%
    ungroup() %>%
    sort()

  if (nuancedTSS) {
    nuancedGenes <- split(altTSS, as.character(altTSS$name)) %>%
      parallel::mclapply(., function(x) {
        tmp <- gaps(x) %>%
          filter(seqnames == seqnames(x) & strand == strand(x)) %>%
          width(.)
        ifelse(length(tmp) == 3 & tmp[2] < nuancedTSSGap, FALSE, TRUE)
      }) %>%
      unlist()

    altTSS <- altTSS %>% filter(name %in% names(nuancedGenes)[nuancedGenes])
  }
  return(altTSS)
}
