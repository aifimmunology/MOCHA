#' Extract fragments by populations from an ArchR Project
#'
#' \code{getPopFrags} returns a list of fragments per cell subset as a GRanges.
#'
#' @param ArchRProj The ArchR Project.
#' @param metaColumn The name of metadata column that contains the populations
#'   of cells you want to merge and export.
#' @param cellSubsets Default is 'ALL'. If you want to export only some groups,
#'   then give it a list of group names. This needs to be unique - no duplicated
#'   names. This list of group names must be identical to names that appear in
#'   the metadata column of the ArchR Project (e.g. metaColumn).
#' @param region Optional parameter. Set this if you only want to extract
#'   fragments from particular regions of the genome. Format should be as a
#'   string (e.g. 'chr1:1000-2000'), or a GRanges object.
#' @param numCores Number of cores to use.
#' @param sampleSpecific Set to TRUE to further subset cells by sample
#' @param NormMethod Normalization method. Can be either "nFrags","nCells", or
#'   "Median".
#' @param blackList Blacklisted region to filter out. Default is to not filter
#'   out anything (i.e. NULL). Input should be provided as a GRanges object. Any
#'   fragments with more than a certain overlap will be thrown out.
#' @param overlapList The minimum overlap necessary for a fragment marked as
#'   overlapping with the blacklist region and thus thrown out.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#'
#' @return A list of GRanges containing fragments. Each GRanges corresponds to a
#'   population defined by cellSubsets (and sample, if
#'   \code{sampleSpecific=TRUE})
#'
#'
#' @export

getPopFrags <- function(ArchRProj,
                        metaColumn,
                        cellSubsets = "ALL",
                        region = NULL,
                        numCores = 1,
                        sampleSpecific = TRUE,
                        NormMethod = "nfrags",
                        blackList = NULL,
                        verbose = FALSE,
                        overlapList = 50) {
  nFrags <- NULL
  # Turn off ArchR logging messages
  ArchR::addArchRVerbose(verbose = FALSE)

  # Extract metadata
  metadf <- ArchR::getCellColData(ArchRProj)

  if (!any(colnames(metadf) %in% metaColumn)) {
    stop(paste(
      stringr::str_interp("Provided metaColumn ${metaColumn} does not exist in the cellColData of your ArchRProj."),
      stringr::str_interp("Available columns are: ${colnames(metadf)}")
    ))
  }

  # Get cell counts for all cell populations
  tmp <- metadf[, metaColumn]
  cellCounts <- table(tmp[!is.na(tmp)])

  if (all(cellSubsets == "ALL")) {
    cellPopulations <- names(cellCounts)
  } else {

    # Take our cell populations from given cellSubsets
    cellPopulations <- unique(cellSubsets[!is.na(cellSubsets)][order(cellSubsets)])

    # Filter cellCounts to selected cell subsets (populations)
    cellCounts <- cellCounts[cellPopulations]
    if (any(is.na(cellCounts))) {
      stop(paste(
        "Some cell populations have NA cell counts.",
        "Please verify that all given cellSubsets exist for the given metaColumn.",
        stringr::str_interp("cellSubsets with NA cell counts: ${cellPopulations[is.na(names(cellCounts))]}")
      ))
    }
  }

  if (length(cellCounts) == 0) {
    stop("No cells were found for the given cellSubsets and/or metaColumn.")
  }

  cellNames <- rownames(metadf)[metadf[, metaColumn] %in% cellPopulations]

  #### Up to here, the above is only specific to sampleSpecific=False
  allArrows <- ArchR::getArrowFiles(ArchRProj)
  # Get sample names from our populations of interest
  samplesToExtract <- paste(
    unique(metadf$Sample[which(metadf[, metaColumn] %in% cellPopulations)]),
    collapse = "|"
  )
  # Get arrows containing the samples we want
  arrows <- allArrows[grepl(samplesToExtract, allArrows)]

  if (length(arrows) == 0) {
    stop(paste(
      "Found no arrows containing samples from your selected cell subset(s).",
      "Check cellSubsets or the input ArchR Project."
    ))
  }

  if (is.null(region)) {

    # Extract fragments from all available regions
    arrowList <- lapply(seq_along(arrows), function(x){ list(arrows[x], cellNames)})
    if (verbose) {
      message(stringr::str_interp("Extracting fragments from arrow files"))        
    }

    fragsList <- pbapply::pblapply(arrowList, simplifiedFragments, cl = numCores)

  } else {

    # Extract fragments from a given region only
    # Validate region and interpret as string or GRanges
    if (validRegionString(region) & tolower(NormMethod) == "raw") {
      regionGRanges <- StringsToGRanges(region)
    } else if (class(region)[1] == "GRanges" & tolower(NormMethod) == "raw") {
      regionGRanges <- region
    } else if (tolower(NormMethod) != "raw") {
      stop(paste(
        "Wrong NormMethod set. Please set NormMethod = 'Raw'",
        "if you wish to extract fragments from a specific region.",
        "\n Alternatively, you can extract the entire genome and then subset."
      ))
    } else {
      stop(paste(
        "Invalid region input.",
        "Region must either be a string matching the format",
        "'seqname:start-end', or a GRanges object."
      ))
    }

    chrom <- regionGRanges %>%
      as.data.frame() %>%
      dplyr::select(.data$seqnames)

    arrowList <- lapply(seq_along(arrows), function(x){ list(arrows[x], cellNames,as.character(chrom[, 1]), regionGRanges)})

    if (verbose) {
        message(stringr::str_interp("Extracting fragments from arrow files.}"))
    }

    fragsList <- pbapply::pblapply(arrowList, simplifiedFragments, cl = numCores) 
  
    fragsList <- lapply(fragsList, function(x){  # Filter according to provided blacklist
      if (is.null(blackList)) {
        fragsGRanges
      } else if (class(blackList)[1] == "GRanges") {
        fragsGRanges %>% plyranges::filter_by_non_overlaps(blackList, minoverlap = overlapList)
      } else {
        stop("Error: Wrong format for blackList region. Please provide GRanges")
      }
    })
  }

  # From MOCHA - sorts cell barcodes by population
  barcodesByCellPop <- lapply(cellPopulations, function(x) {
    row.names(metadf)[which(metadf[, metaColumn] == x)]
  })

  # Add normalization factor.
  if (tolower(NormMethod) == "raw") {
    names(barcodesByCellPop) <- paste(gsub("#", "_", gsub(" |_", "_", cellPopulations)), 1, sep = "__")
  } else if (tolower(NormMethod) == "ncells") {
    names(barcodesByCellPop) <- paste(gsub("#", "_", gsub(" |_", "_", cellPopulations)), cellCounts / 1000, sep = "__")
  } else if (tolower(NormMethod) == "nfrags") {
    # Calculate the total nFrags for our cell populations
    nFragsNorm <- as.data.frame(metadf[, c(metaColumn, "nFrags")]) %>%
      dplyr::filter(!is.na(get(metaColumn))) %>%
      dplyr::filter(get(metaColumn) %in% cellPopulations) %>%
      dplyr::group_by(get(metaColumn)) %>%
      dplyr::summarize(nFrags = sum(nFrags))

    # Verify expected cellPopulations have corresponding nFrags
    if (!all(nFragsNorm[, 1] == cellPopulations)) {
      stop("Names of nFrags don't match the cell populations.", nFragsNorm[, 1])
    }
    names(barcodesByCellPop) <- paste(gsub("#", "_", gsub(" |_", "_", cellPopulations)), nFragsNorm$nFrags / 10^6, sep = "__")
  } else if (tolower(NormMethod) == "median") {
    names(barcodesByCellPop) <- paste(gsub("#", "_", gsub(" |_", "_", cellPopulations)), "Median", sep = "__")
  } else if (tolower(NormMethod) == "medmax") {
    names(barcodesByCellPop) <- paste(gsub("#", "_", gsub(" |_", "_", cellPopulations)), "MedianMax", sep = "__")
  } else if (tolower(NormMethod) == "TotalSampleByNCells") {
    # This calculates total fragments per sample.
    totalsampleFrags <- unlist(fragsList, length)
    names(barcodesByCellPop) <- paste(gsub("#", "_", gsub(" |_", "_", cellPopulations)), totalsampleFrags / 10^6 * cellCounts / 1000, sep = "__")
  } else {
    stop("Error: Incorrect NormMethod given.")
  }

  ## Identify which subset of arrows the population can be found it.
  ## Speeds up sample-specific population extraction
  fragsListIndex <- lapply(cellPopulations, function(x) {
    names(arrows) %in% unique(metadf$Sample[which(metadf[, metaColumn] == x)])
  })

  # Sort fragments into a list by cell population

  popFrags <- lapply(seq_along(barcodesByCellPop), function(x) {
    if (verbose) {
      message("Extracting fragments for cellPopulation__normalization: ", names(barcodesByCellPop)[x])
    }
    if (sum(fragsListIndex[[x]]) > 1) {
      
      fragIterList <- lapply(which(fragsListIndex[[x]]), function(x){
        list(barcodesByCellPop[[x]], fragsList[[y]])
      })

      tmp <- pbapply::pblapply(fragIterList, subset_Frag, cl = numCores)

    } else {
      tmp <- list(subset_Frag(list(barcodesByCellPop[[x]], fragsList[[which(fragsListIndex[[x]])]])))
    }

    # For this population, get sample-specific normalization factors
    # and rename the fragment GRanges
    if (sampleSpecific) {
      subSet_tmp <- gsub("__.*", "", names(barcodesByCellPop)[x])

      if (tolower(NormMethod) == "ncells") {
        tmp_cellCount <- unlist(lapply(tmp, function(y) {
          length(unique(y$RG))
        }))

        names(tmp) <- paste(
          subSet_tmp,
          "#",
          names(arrows)[unlist(fragsListIndex[[x]])],
          "__",
          tmp_cellCount / 1000,
          sep = ""
        )
      } else if (tolower(NormMethod) == "nfrags") {
        names(tmp) <- paste(
          subSet_tmp,
          "#",
          names(arrows)[unlist(fragsListIndex[[x]])],
          "__",
          unlist(lapply(tmp, length)) / 10^6,
          sep = ""
        )
      } else {
        subSet_tmp <- gsub("__.*", "", names(barcodesByCellPop)[x])
        norm_tmp <- gsub(".*__", "", names(barcodesByCellPop)[x])

        names(tmp) <- paste(
          subSet_tmp,
          "#",
          names(arrows)[unlist(fragsListIndex[[x]])], "__", norm_tmp,
          sep = ""
        )
      }

      tmp
    } else {
      IRanges::stack(methods::as(tmp, "GRangesList"))
    }
  })

  if (sampleSpecific) {
    popFrags <- unlist(popFrags, recursive = FALSE)
  } else {
    names(popFrags) <- names(barcodesByCellPop)
  }

  # Turn off ArchR logging messages
  ArchR::addArchRVerbose(verbose = TRUE)

  return(popFrags)
}

#' Extract fragments from an arrow file based on one variable
#'
#' \code{simplifiedFragments} returns a list of fragments for a given set of cell names
#'
#' @param ref a list where the first index if the name of the arrow and the second index is a vector of strings describing cell names
#'
#' @return A Granges object for all fragments from a set of cells within a given arrow file. 
#'
#'
#' @noRd

simplifiedFragments <- function(ref){
   arrows <- ref[[1]]
   cellNames <- ref[[2]]

   if(length(ref) > 2){
    regionGRanges = ref[[4]]
    chrom = ref[[3]]
    frags <- ArchR::getFragmentsFromArrow(
        ArrowFile = arrows,
        cellNames = cellNames,
        chr = chrom,
        verbose = FALSE
      ) 
    frags <- plyranges::filter_by_overlaps(frags, regionGRanges)
    return(frags)

   }else{

    frags <- ArchR::getFragmentsFromArrow(
        ArrowFile = arrows,
        cellNames = cellNames,
        verbose = FALSE
      )
    return(frags)
   }

}

#' subsets fragments out by cellnames. 
#'
#' \code{subset_Frag} returns a sorted set of fragments by populations based on cell barcode lists
#'
#' @param ref a list where the first index is the cellnames for that population, and the second index is a GRanges of fragments
#'
#' @return A Granges object for all fragments from a set of cells within a given arrow file. 
#'
#'
#' @noRd

# From MOCHA - Function to sort fragments by populations based on cell barcode lists
subset_Frag <- function(ref) {
  cellNames = ref[[1]]
  fragsGRanges = ref[[2]]
  fragsTable <- as.data.table(fragsGRanges)
  idx <- which(fragsTable$RG %in% cellNames)
  return(fragsGRanges[idx])
}
