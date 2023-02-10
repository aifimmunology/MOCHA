#' @title \code{addMotifSet}
#'
#' @description \code{addMotifSet} Identify motifs within your peakset.
#'
#' @param SampleTileMatrix A SummarizedExperiment, specifically the output of
#'   getSampleTileMatrix
#' @param motifPWMs A pwms object for the motif database. Either PFMatrix,
#'   PFMatrixList, PWMatrix, or PWMatrixList
#' @param genome BSgenome object, DNAStringSet, or FaFile, or short string
#'   signifying genome build recognized by [BSgenome::getBSgenome()]
#' @param w Parameter for motifmatchr controlling size in basepairs of window for filtration.
#'   Default is 7.
#' @param returnSTM If TRUE, return the modified SampleTileMatrix with motif set
#'   added to metadata (default). If FALSE, return just the motifs from motifmatchr.
#' @param motifSetName Name to give motifList in the SampleTileMatrix's metadata
#'   if `returnSTM=TRUE`. Default is 'Motifs'.
#'
#' @return the modified SampleTileMatrix with motifs added to the metadata
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' # load a curated motif set from library(chromVARmotifs)
#' # included with ArchR installation
#' data(human_pwms_v2)
#' SE_with_motifs <- addMotifSet(
#'   SampleTileMatrix,
#'   motifPWMs = human_pwms_v2,
#'   returnSTM = TRUE, motifSetName = "Motifs", w = 7
#' )
#' }
#'
#' @export
#'
addMotifSet <- function(SampleTileMatrix, 
                        motifPWMs, 
                        genome = NULL,
                        w = 7, 
                        returnSTM = TRUE, 
                        motifSetName = "Motifs") {
  if (!requireNamespace("motifmatchr", quietly = TRUE)) {
    stop(
      "Package 'motifmatchr' is required for addMotifSet. ",
      "Please install 'motifmatchr' to proceed."
    )
  }
  TotalPeakSet <- SummarizedExperiment::rowRanges(SampleTileMatrix)
  
  if (is.null(genome)) {
    warning("Genome not provided - looking for genome in metadata of the ",
            "provided SampleTileMatrix")
    genome <- S4Vectors::metadata(SampleTileMatrix)$Genome
  }
  motif_ix <- motifmatchr::matchMotifs(
    pwms = motifPWMs,
    subject = TotalPeakSet,
    genome = genome,
    out = "positions", w = w)

  . <- NULL
  names(motif_ix) <- sub("_D_.*|_I_.*", "", names(motif_ix)) %>%
    sub("_I$|_D$", "", .) %>%
    sub(".*_LINE", "", .) %>%
    sub(".*_", "", .)

  motifList <- list(motif_ix)
  names(motifList) <- motifSetName

  if (returnSTM) {
    S4Vectors::metadata(SampleTileMatrix) <- append(
      S4Vectors::metadata(SampleTileMatrix), motifList
    )
    return(SampleTileMatrix)
  } else {
    return(motif_ix)
  }
}
