# Function to get sample-level metadata,
# from an ArchR project's colData
sampleDataFromCellColData <- function(cellColData, sampleLabel) {
  if (!(sampleLabel %in% colnames(cellColData))) {
    stop(paste(
      "`sampleLabel` must present in your ArchR Project's cellColData",
      "Check `names(getCellColData(ArchRProj)` for possible sample columns."
    ))
  }

  # Drop columns where all values are NA
  cellColDataNoNA <- BiocGenerics::Filter(function(x) {
    !all(is.na(x))
  }, cellColData)

  # Convert to data.table
  cellColDT <- data.table::as.data.table(cellColDataNoNA)

  BoolDT <- cellColDT[, lapply(.SD, function(x) {
    length(unique(x)) == 1
  }), by = c(sampleLabel)]
  trueCols <- apply(BoolDT, 2, all)
  trueCols[[sampleLabel]] <- TRUE
  cellColDF <- as.data.frame(cellColDT)

  sampleData <- dplyr::distinct(cellColDF[, names(which(trueCols)), drop = F])

  # Set sampleIDs as rownames
  rownames(sampleData) <- sampleData[[sampleLabel]]
  return(sampleData)
}


# Function to split the output of getPopFrags into a list
# of lists of GRanges, one named for each celltype.
# Resulting in a list containing each celltype, and each
# celltype has a list of GRanges name for each sample.
splitFragsByCellPop <- function(frags) {

  # Rename frags by cell population
  renamedFrags <- lapply(
    1:length(frags),
    function(y) {
      # Split out celltype and sample from the name
      x <- frags[y]
      celltype_sample <- names(x)
      splits <- unlist(stringr::str_split(celltype_sample, "#"))
      celltype <- splits[1]
      sample <- unlist(stringr::str_split(splits[2], "__"))[1]
      # Rename the fragments with just the sample
      names(x) <- sample
      # Return as a list named for celltype
      output <- list(x)
      names(output) <- celltype
      output
    }
  )

  # Group frags by cell population
  renamedFrags <- unlist(renamedFrags, recursive = FALSE)
  splitFrags <- split(renamedFrags, f = names(renamedFrags))
  return(splitFrags)
}

# Function to generate parallelization. 
makeMOCHACluster <- function(numCores = 1) {

  if (numCores > 1) {
    #if(.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(numCores)
    #} else {
   #   cl <- parallel::makeCluster(numCores, type = "FORK") # Use forking on unix
    #}
  } else {
    # Set numCores = 1 (or <= 1) for sequential
    # evaluation with pblapply
    cl <- NULL
  }

  return(cl)
}



# Tests if a string is a in the correct format to convert to GRanges
validRegionString <- function(regionString) {
  if (!is.character(regionString)) {
    return(FALSE)
  }

  pattern <- "([0-9]{1,2}|chr[0-9]{1,2}|chr[X-Y]{1,1}):[0-9]*-[0-9]*"
  matchedPattern <- stringr::str_extract(regionString, pattern)

  if (any(is.na(matchedPattern))) {
    return(FALSE)
  } else if (any(!matchedPattern == regionString)) {
    return(FALSE)
  }

  splits <- stringr::str_split(regionString, "[:-]")[[1]]
  start <- splits[2]
  end <- splits[3]
  if (any(start > end)) {
    return(FALSE)
  }

  # All conditions satisfied
  return(TRUE)
}

#' @title \code{StringsToGRanges}
#'
#' @description \code{StringsToGRanges} Turns a list of strings in the format chr1:100-200
#'   into a GRanges object
#'
#' @param regionString A string or list of strings each in the format chr1:100-200
#' @return a GRanges object with ranges representing the input string(s)
#'
#' @export
StringsToGRanges <- function(regionString) {
  # if (length(regionString)>1){
  #   boolList <- lapply(regionString, function(x){validRegionString(x)})
  #   if (!all(boolList)){
  #     stop("Some region strings are invalid. Given regions must all be strings matching format 'seqname:start-end', where start<end e.g. chr1:123000-123500")
  #   }
  # } else if(!validRegionString(regionString)) {
  #   stop("Region must be a string matching format 'seqname:start-end', where start<end e.g. chr1:123000-123500")
  # }
  . <- NULL
  chrom <- gsub(":.*", "", regionString)
  startSite <- gsub(".*:", "", regionString) %>%
    gsub("-.*", "", .) %>%
    as.numeric()
  endSite <- gsub(".*-", "", regionString) %>% as.numeric()

  if (any(startSite >= endSite)) {
    stop("Error in region string: Make sure the start of the genomic range occurs before the end")
  }
  regionGRanges <- GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::IRanges(start = startSite, end = endSite), strand = "*")
  return(regionGRanges)
}

#' @title \code{GRangesToString} Converts a GRanges object to a string in the format 'chr1:100-200'
#'
#' @description \code{GRangesToString} Turns a GRanges Object into
#'  a list of strings in the format chr1:100-200
#'
#' @param GR_obj the GRanges object to convert to a string
#' @return A string or list of strings in the format 'chr1:100-200' representing
#'  ranges in the input GRanges
#'
#' @export
GRangesToString <- function(GR_obj) {
  paste(GenomicRanges::seqnames(GR_obj), ":", GenomicRanges::start(GR_obj), "-", GenomicRanges::end(GR_obj), sep = "")
}

#' @title \code{differentialsToGRanges} Converts a data.frame matrix to a GRanges,
#'   preserving additional columns as GRanges metadata
#'
#' @param differentials a matrix/data.frame with a column tileColumn containing
#'   region strings in the format "chr:start-end"
#' @param tileColumn name of column containing region strings. Default is "Tile".
#'
#' @return a GRanges containing all original information
#' @export
differentialsToGRanges <- function(differentials, tileColumn = "Tile") {
  regions <- MOCHA::StringsToGRanges(differentials[[tileColumn]])
  GenomicRanges::mcols(regions) <- differentials
  regions
}
