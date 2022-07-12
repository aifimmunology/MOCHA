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
#'
#' @return A list of GRanges containing fragments. Each GRanges corresponds to a
#'   population defined by cellSubsets (and sample, if
#'   \code{sampleSpecific=TRUE})
#'
getPopFrags <- function(ArchRProj, metaColumn, cellSubsets = "ALL", region = NULL, numCores = 1, sampleSpecific = TRUE,
                        NormMethod = "nfrags", blackList = NULL, overlapList = 50) {
  # Extract metadata
  metadf <- getCellColData(ArchRProj)

  if (!any(colnames(metadf) %in% metaColumn)) {
    stop(paste(
      str_interp("Provided metaColumn ${metaColumn} does not exist in the cellColData of your ArchRProj."),
      str_interp("Available columns are: ${colnames(metadf)}")
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
        str_interp("cellSubsets with NA cell counts: ${cellPopulations[is.na(names(cellCounts))]}")
      ))
    }
  }

  if (length(cellCounts) == 0) {
    stop("No cells were found for the given cellSubsets and/or metaColumn.")
  }


  #### Up to here, the above is only specific to sampleSpecific=False
  allArrows <- getArrowFiles(ArchRProj)
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
    fragsList <- parallel::mclapply(seq_along(arrows), function(x) {
      getFragmentsFromArrow(
        ArrowFile = arrows[x],
        cellNames = NULL,
        verbose = FALSE
      )
    }, mc.cores = numCores)
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
      dplyr::select(seqnames)
    fragsList <- parallel::mclapply(seq_along(arrows), function(x) {
      fragsGRanges <- getFragmentsFromArrow(
        ArrowFile = arrows[x],
        cellNames = NULL,
        chr = as.character(chrom[, 1]),
        verbose = FALSE
      ) %>% plyranges::join_overlap_intersect(regionGRanges)

      # Filter according to provided blacklist
      if (is.null(blackList)) {
        fragsGRanges
      } else if (class(blackList)[1] == "GRanges") {
        fragsGRanges %>% plyranges::filter_by_non_overlaps(blackList, minoverlap = overlapList)
      } else {
        stop("Error: Wrong format for blackList region. Please provide GRanges")
      }
    }, mc.cores = numCores)
  }

  # From scMACS - sorts cell barcodes by population
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
      print(nFragsNorm[, 1])
      stop("Names of nFrags don't match the cell populations.")
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

  # From scMACS - Function to sort fragments by populations based on cell barcode lists
  subset_Frag <- function(cellNames, fragsGRanges) {
    fragsTable <- as.data.table(fragsGRanges)
    idx <- which(fragsTable$RG %in% cellNames)
    fragsGRanges[idx]
  }

  ## Identify which subset of arrows the population can be found it.
  ## Speeds up sample-specific population extraction
  fragsListIndex <- lapply(cellPopulations, function(x) {
    names(arrows) %in% unique(metadf$Sample[which(metadf[, metaColumn] == x)])
  })

  # Sort fragments into a list by cell population

  popFrags <- lapply(seq_along(barcodesByCellPop), function(x) {
    print(paste("Sorting ", names(barcodesByCellPop)[x], sep = ""))
    if (sum(fragsListIndex[[x]]) > 1) {
      tmp <- parallel::mclapply(which(fragsListIndex[[x]]), function(y) {
        subset_Frag(barcodesByCellPop[[x]], fragsList[[y]])
      }, mc.cores = 20)
    } else {
      tmp <- list(subset_Frag(
        barcodesByCellPop[[x]],
        fragsList[[which(fragsListIndex[[x]])]]
      ))
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
      stack(as(tmp, "GRangesList"))
    }
  })

  if (sampleSpecific) {
    popFrags <- unlist(popFrags, recursive = FALSE)
  } else {
    names(popFrags) <- names(barcodesByCellPop)
  }
  return(popFrags)
}