#' @title \code{addMotifSet}
#'
#' @description \code{addMotifSet}Identify motifs within peakset
#'
#' @param SE_Object your MOCHA SummarizedExperiment. Requires Genome AnnotationDbi object within the metadata added by getSampleTileMatrix
#' @param pwms a pwms object for the motif database. Either PFMatrix, PFMatrixList, PWMatrix, or PWMatrixLis'
#' @param w the width for motifmatchr
#' @param returnObj if TRUE, return the modified SE_Object with motif set added to metadata (default). If FALSE, return the motifs from motifmatchr.
#' @param motifSetName name of the motifList in the SE_object's metadata if returnObj=TRUE. Default is 'Motifs'.
#'
#' @return the modified SE_Object with motifs added to the metadata
#' @examples
#' \dontrun{
#' # load a curated motif set from library(chromVARmotifs) included with ArchR installation
#' data(human_pwms_v2)
#' SE_with_motifs <- addMotifSet(SE_Object, pwms = human_pwms_v2, returnObj = TRUE, motifSetName = "Motifs", w = 7)
#' }
#'
#' @export

addMotifSet <- function(SE_Object, pwms, w = 7, returnObj = TRUE, motifSetName = "Motifs") {
  TotalPeakSet <- rowRanges(SE_Object)
  genome <- metadata(SE_Object)$Genome
  motif_ix <- motifmatchr::matchMotifs(
    pwms = pwms,
    TotalPeakSet,
    genome = genome,
    out = "positions", w = w
  )
  names(motif_ix) <- sub("_D_.*|_I_.*", "", names(motif_ix)) %>%
    sub("_I$|_D$", "", .) %>%
    sub(".*_LINE", "", .) %>%
    sub(".*_", "", .)

  motifList <- list(motif_ix)
  names(motifList) <- motifSetName

  if (returnObj) {
    metadata(SE_Object) <- append(metadata(SE_Object), motifList)
    return(SE_Object)
  } else {
    return(motif_ix)
  }
}
