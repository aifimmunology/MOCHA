#' @title \code{annotateTiles}
#'
#' @description \code{annotateTiles} annotates a set of sample-tile matrices
#'   given with gene annotations. Details on TxDb and Org annotation packages
#'   and available annotations can be found at Bioconductor:
#'   https://bioconductor.org/packages/3.15/data/annotation/
#'
#' @param SampleTileObj A RangedSummarizedExperment generated from getSampleTileMatrix,
#'   containing TxDb and Org in the metadata. This may also be a GRanges object.
#' @param TxDb The annotation package for TxDb object for your genome.
#'   Optional, only required if SampleTileObj is a GRanges.
#' @param Org The genome-wide annotation for your organism.
#'   Optional, only required if SampleTileObj is a GRanges.
#' @param promoterRegion Optional list containing the window size in basepairs
#' defining the promoter region. The format is (upstream, downstream).
#' Default is (2000, 100).
#'
#' @return SampleTileObj the input data structure with added gene annotations.
#'
#' @examples
#' \dontrun{
#' SampleTileMatricesAnnotated <- MOCHA::annotateTiles(
#'   SampleTileMatrices
#' )
#' }
#'
#' @export
annotateTiles <- function(SampleTileObj,
                          TxDb = NULL,
                          Org = NULL,
                          promoterRegion = c(2000, 100)) {
  if (class(SampleTileObj)[1] == "RangedSummarizedExperiment" & is.null(TxDb) & is.null(Org)) {
    if (!all(c("TxDb", "Org") %in% names(S4Vectors::metadata(SampleTileObj)))) {
      stop("Error: TxDb and/or Org are missing from SampleTileObj. SampleTileObj as a RangedSummarizedExperiment must contain a TxDb and Org in the metadata.")
    }
    tileGRanges <- SummarizedExperiment::rowRanges(SampleTileObj)
    TxDb <- AnnotationDbi::loadDb(S4Vectors::metadata(SampleTileObj)$TxDb)
    Org <- AnnotationDbi::loadDb(S4Vectors::metadata(SampleTileObj)$Org)
  } else if (class(tileList)[[1]] == "GRanges" & !is.null(TxDb) & !is.null(Org)) {
    tileGRanges <- SampleTileObj
  } else {
    stop("Error: Invalid inputs. Verify SampleTileObj is a RangedSummarizedExperiment. If SampleTileObj is a GRanges, TxDb and Org must be provided.")
  }

  txList <- suppressWarnings(GenomicFeatures::transcriptsBy(TxDb, by = ("gene")))
  names(txList) <- suppressWarnings(AnnotationDbi::mapIds(Org, names(txList), "SYMBOL", "ENTREZID"))
  exonList <- suppressWarnings(ensembldb::exonsBy(TxDb, by = "gene"))
  # Same TxDb, same gene names in same order.
  names(exonList) <- names(txList)

  txs <- stack(txList) %>%
    GenomicRanges::trim() %>%
    S4Vectors::unique(.)
  exonSet <- stack(exonList)

  promoterSet <- stack(txList) %>%
    GenomicRanges::trim(.) %>%
    S4Vectors::unique(.) %>%
    suppressWarnings(GenomicRanges::promoters(., upstream = promoterRegion[1], downstream = promoterRegion[2]))

  getOverlapNameList <- function(rowTiles, annotGR) {
    overlapGroup <- findOverlaps(rowTiles, annotGR) %>% as.data.frame()
    overlapGroup$Genes <- as.character(annotGR$name[overlapGroup$subjectHits])
    last <- overlapGroup %>%
      dplyr::group_by(subjectHits) %>%
      dplyr::summarize(Genes = paste(unique(Genes), collapse = ", "))
    return(last)
  }

  txs_overlaps <- getOverlapNameList(tileGRanges, txs)
  exon_overlaps <- getOverlapNameList(tileGRanges, exonSet)
  promo_overlaps <- getOverlapNameList(tileGRanges, promoterSet)

  tileType <- as.data.frame(mcols(tileGRanges)) %>%
    dplyr::mutate(Index = 1:nrow(.)) %>%
    dplyr::mutate(Type = dplyr::case_when(
      Index %in% promo_overlaps$subjectHits ~ "Promoter",
      Index %in% exon_overlaps$subjectHits ~ "Exonic",
      Index %in% txs_overlaps$subjectHits ~ "Intronic",
      TRUE ~ "Distal"
    )) %>%
    dplyr::left_join(promo_overlaps, by = c("Index" = "subjectHits")) %>%
    dplyr::rename("Promo" = Genes) %>%
    dplyr::left_join(exon_overlaps, by = c("Index" = "subjectHits")) %>%
    dplyr::rename("Exons" = Genes) %>%
    dplyr::left_join(txs_overlaps, by = c("Index" = "subjectHits")) %>%
    dplyr::rename("Txs" = Genes) %>%
    dplyr::mutate(Genes = ifelse(Type == "Promoter", Promo, NA)) %>%
    dplyr::mutate(Genes = ifelse(Type == "Exonic", Exons, Genes)) %>%
    dplyr::mutate(Genes = ifelse(Type == "Intronic", Txs, Genes))

  tileGRanges$tileType <- tileType$Type
  tileGRanges$Gene <- tileType$Genes

  # If input was as Ranged SE, then edit the rowRanges for the SE and return it.
  # Else, return the annotated tile GRanges SampleTileObject.
  if (class(SampleTileObj)[1] == "RangedSummarizedExperiment") {
    SummarizedExperiment::rowRanges(SampleTileObj) <- tileGRanges
    return(SampleTileObj)
  } else {
    return(tileGRanges)
  }
}
