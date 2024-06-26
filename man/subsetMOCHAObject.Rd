% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subsetMOCHAObject.R
\name{subsetMOCHAObject}
\alias{subsetMOCHAObject}
\title{Subset a tileResults object by metadata}
\usage{
subsetMOCHAObject(
  Object,
  subsetBy,
  groupList,
  removeNA = TRUE,
  subsetPeaks = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{Object}{A MultiAssayExperiment or RangedSummarizedExperiment,}

\item{subsetBy}{The variable to subset by. Can either be 'celltype', or a
column from the sample metadata (see `colData(Object)`).}

\item{groupList}{the list of cell type names or sample-associated data that
should be used to subset the Object}

\item{removeNA}{If TRUE, removes groups in groupList that are NA. If FALSE,
keep groups that are NA.}

\item{subsetPeaks}{If `subsetBy` = 'celltype', subset the tile set to tiles
only called in those cell types. Default is TRUE.}

\item{verbose}{Set TRUE to display additional messages. Default is FALSE.}
}
\value{
Object the input Object, filtered down to either the cell type or
  samples desired.
}
\description{
\code{subsetMOCHAObject} subsets a tileResults-type object (from
  callOpenTiles), or a SummarizedExperiment-type object (from
  getSampleTileMatrix), either by cell type or sample metadata.
}
\keyword{utils}
