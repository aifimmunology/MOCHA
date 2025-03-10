% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/motifFootprint.R
\name{motifFootprint}
\alias{motifFootprint}
\title{Generate motif footprints}
\usage{
motifFootprint(
  SampleTileObj,
  motifName = "Motifs",
  specMotif = NULL,
  regions = NULL,
  cellPopulations = "ALL",
  windowSize = 500,
  normTn5 = TRUE,
  smoothTn5 = 10,
  groupColumn = NULL,
  subGroups = NULL,
  sampleSpecific = FALSE,
  numCores = 1,
  force = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{SampleTileObj}{A Sample-tile object from MOCHA's getSampleTileMatrix()}

\item{motifName}{The name of metadata entry with Motif location, added via addMotifSet()}

\item{specMotif}{An optional string specifying which motif to analyze. If blank, it will analyze all motifs.}

\item{regions}{An optional GRanges object or list of strings in the format chr1:100-200, specifiying which specific regions to look at when conducting motif footprinting.}

\item{cellPopulations}{A list of cell populations to conduct motif footprinting on.}

\item{windowSize}{A number, representing the window to analyze around each motif location. Default is 500 bp.}

\item{normTn5}{A boolean for whether to normalize by Tn5 insertion bias.}

\item{smoothTn5}{The window size for smoothen Tn5 insertions at each location. Ideal when looking at motif footprints over a smaller number of regions, rarer cell types, or sparse regions, where local noise can make it harder to see the overall pattern. Can be set to 0.}

\item{groupColumn}{A string, corresponding to a metadata column within the SampleTileObj, that describes the groups by which you want to summarize motif footprints. If sampleSpecific = FALSE, then motifFootprint will average insertions across samples within each group.}

\item{subGroups}{A list of subgroups, if you want to only look at specific groups within the groupColumn.}

\item{sampleSpecific}{A boolean for whether to generate average motif footprints within each group, or to return a data.frame for all samples.}

\item{numCores}{Number of cores to parallelize over}

\item{force}{Boolean. If FALSE, it will through an error if there's an empty sample or not overlap with regions and motifs. If TRUE will ignore these issues and continue (or return NULL)}

\item{verbose}{Boolean. Default is FALSE. Will print more messages if TRUE.}
}
\value{
A SummarizedExperiment containing motif footprinting data
}
\description{
Generate a plot of average normalized insertions around a motif center
}
\keyword{downstream}
