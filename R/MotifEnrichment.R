#' @title Test for enrichment of motifs against a background
#'
#' @description Test for enrichment of motifs within Group1 against a background
#'   Group2 using a hypergeometric t-test.
#'
#' @param Group1 A GRanges object, such as a set of significant differential
#'   tiles.
#' @param Group2 A GRanges object containing background regions, non-overlapping
#'   with Group1
#' @param motifPosList A GRangesList of motifs and positions for each motif.
#'   Must be named for each motif.
#' @param type Optional, name of a metadata column in Group1 and Group2 to test
#'   for enrichment the number of unique entries in column given by 'type'.
#'   Default is NULL, which tests the number of Ranges.
#'
#' @return A data.frame containing enrichment for each group
#'
#' @export
#' @keywords downstream
MotifEnrichment <- function(Group1, Group2, motifPosList, type = NULL) {
  allEnrichmentList <- pbapply::pblapply(motifPosList, function(x) {
    EnrichedRanges(
      Group1,
      Group2,
      Category = x,
      type = type,
      returnTable = FALSE
    )
  }, cl = 1) # cl = 1 to disable multiprocessing as it is not much faster
  df_final <- do.call("rbind", allEnrichmentList)

  df_final$adjp_val <- stats::p.adjust(df_final$p_value, method = "fdr")
  df_final$mlog10Padj <- -log10(df_final$adjp_val)

  return(df_final)
}

#' @title Get motifset for a given cell type
#'
#' @description Extracts an existing motif set and filters it down to only motifs found within regions open in that celltype
#'                A necessary step for accurate motif enrichment analysis. 
#'
#' @param STM A Sample-Tile Object from MOCHA
#' @param cellPopulation The name of a cell type within MOCHA
#' @param MotifSetName The name of the total motif set, saved as metadata within the MOCHA object
#' @param specMotif Name of a specific motif or sets of motif that you wish to pull out. Default is NULL, with returns all motifs. 
#' @param asGRangesList Default is FALSE. If TRUE, then converts list of GRanges to GRangesList, which takes quite a while. 
#' 
#'
#' @return A GRangesList containing the motifs that are present in open regions for a given cell type
#'
#' @export
#' @keywords downstream
getCellTypeMotifs <- function(STM, cellPopulation, MotifSetName = 'Motifs', specMotif = NULL, 
                              asGRangesList = TRUE) {
  if(!MotifSetName %in% names(STM@metadata)){
    stop('MotifSetName not found within STM Object.')
  }
  allMotifs <- STM@metadata[[MotifSetName]]
    
  if(!is.null(specMotif)){
      if(all(specMotif %in% names(allMotifs))){
          
          allMotifs <- allMotifs[names(allMotifs) %in% specMotif]
          
      }else{
      
          stop('specMotif not found within motif names. Make sure you have an exact match to the motif name')
          
      }
    }
      
    
  #Pull out background tiles. 
  backgroundTiles = getCellTypeTiles(STM, cellPopulation)
  minWidth = min(unlist(lapply(allMotifs, GenomicRanges::width)))
  subMotifs <- lapply(allMotifs, function(MM){
              plyranges::filter_by_overlaps(MM, backgroundTiles, minoverlap = minWidth)
              # minimum motif size is 6. 
      })
  if(asGRangesList){
    subMotifs <- methods::as('GRangesList', subMotifs)
  }
  return(subMotifs)
}

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
#' @param returnTable If TRUE, return the table used for the hypergeometric
#'   enrichment test. Default is FALSE, return the p-value result from the
#'   hypergeometric test.
#'
#' @return A data.frame table of enrichment
#'
#' @noRd
#' @keywords internal
EnrichedRanges <- function(Group1,
                           Group2,
                           Category,
                           type = NULL,
                           returnTable = FALSE) {
  if (!methods::is(Group1, "GRanges") || !methods::is(Group2, "GRanges")) {
    stop(
      "Input Group1 and/or Group2 are not class 'GRanges'. Group1 and ",
      "Group2 must be GRanges objects."
    )
  }
  Group1Cat <- plyranges::filter_by_overlaps(Group1, Category)
  Group2Cat <- plyranges::filter_by_overlaps(Group2, Category)

  OnlyGroup1 <- plyranges::filter_by_non_overlaps(Group1, Category)
  OnlyGroup2 <- plyranges::filter_by_non_overlaps(Group2, Category)

  # Conduct hypergeometric test and return p-value+enrichment score
  if (is.null(type)) {
    if (returnTable) {
      dt_table <- data.frame(
        Group1 = c(length(Group1Cat), length(OnlyGroup1)),
        Group2 = c(length(Group2Cat), length(OnlyGroup2)),
        row.names = c("In Category", "Not in Category")
      )
      return(t(dt_table))
    }

    pVal <- stats::phyper(
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
    if (returnTable) {
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
    }

    pVal <- stats::phyper(
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
    stop(
      "Invalid 'type' or Group1 or Group2 input. Given 'type' must be ",
      "the name of a column in the metadata of both Group1 and Group2."
    )
  }

  return(data.frame(p_value = pVal, enrichment = enrichment))
}

