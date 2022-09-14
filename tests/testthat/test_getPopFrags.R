# Check for existence of ArchR test data first:
capture.output(ArchR::addArchRThreads(threads = 10), type = "message")
withr::local_options(timeout = 600) # 10 minute timeout for download
capture.output(ArchR::addArchRVerbose(verbose = FALSE), type = "message")


# Download (first time only) and load ArchR project
capture.output(testProj <- ArchR::getTestProject(), type = "message") # Requires an internet connection
metaColumn <- "Clusters" # Column with cell population groups

# # Test getPopFrags with full genome and all normalization methods
# # Includes snapshot test for output
# for (sampleSpecific in c(TRUE, FALSE)) {
#   for (NormMethod in c("raw", "nFrags", "nCells", "Median", NULL)) {
#     test_that(
#       stringr::str_interp(
#         "getPopFrags works for NormMethod='${NormMethod}', sampleSpecific=${sampleSpecific}"
#       ),
#       {
#         testthat::local_edition(3)
#         capture.output(
#           popFrags <- getPopFrags(
#             testProj,
#             metaColumn = metaColumn,
#             numCores = 1,
#             NormMethod = NormMethod,
#             sampleSpecific = sampleSpecific
#           ),
#           type = "message"
#         )
#         
#         # Test population+sample names are as expected
#         populations <- sort(unique(getCellColData(testProj)[[metaColumn]]))
#         
#         if (sampleSpecific) {
#           samples <- unique(testProj$Sample)
#           combinations <- purrr::cross2(populations, samples) %>% purrr::map(purrr::lift(paste))
#           expected_names <- gsub(" ", "#", combinations)
#         } else {
#           expected_names <- populations
#         }
#         sampleNames <- names(popFrags)
#         actual_names <- sort(gsub("__.*", "", sampleNames))
#         
#         expect_equal(actual_names, expected_names)
#         
#         # Remove population+sample names and check that the content is equal
#         names(popFrags) <- NULL
#         
#         snapshotVariant <- str_interp("${NormMethod}")
#         if (sampleSpecific) {
#           snapshotVariant <- paste(snapshotVariant, "sampleSpecific", sep = "_")
#         }
#         expect_snapshot_output(
#           popFrags,
#           variant = snapshotVariant
#         )
#       }
#     )
#   }
# }


# Test getPopFrags with specific region
test_that(
  "getPopFrags works with a specific region when NormMethod='Raw'",
  {
    capture.output(
      popFrags <- getPopFrags(
        testProj,
        metaColumn = metaColumn,
        numCores = 1,
        region = "chr2:1-187350807",
        NormMethod = "Raw",
        sampleSpecific = FALSE
      ),
      type = "message"
    )

    # Test population names are as expected
    actual_names <- sort(gsub("__.*", "", names(popFrags)))
    expected_names <- sort(unique(getCellColData(testProj)[[metaColumn]]))
    expect_equal(actual_names, expected_names)
  }
)

# Test getPopFrags with the wrong NormMethod when specifying a region
for (NormMethod in c("nFrags", "nCells", "Median", NULL)) {
  test_that(
    str_interp("getPopFrags throws an error if the wrong NormMethod (${NormMethod}) is set when asking for a specific region"),
    {
      expect_error(
        capture.output(popFrags <- getPopFrags(
          testProj,
          metaColumn = metaColumn,
          numCores = 1,
          region = "chr2:1-187350807",
          NormMethod = NormMethod,
          sampleSpecific = FALSE
        ), type = "message"),
        "Wrong NormMethod"
      )
    }
  )
}

# Test validRegionString
test_that(
  "validRegionString works on edge cases",
  {
    expect_true(MOCHA:::validRegionString("chr1:1-123412"))
    expect_true(MOCHA:::validRegionString("chr20:12300-12400"))
    expect_true(MOCHA:::validRegionString("4:145-146"))
    expect_false(MOCHA:::validRegionString("chr20_12300-12400"))
    expect_false(MOCHA:::validRegionString("chr20:12300-12200"))
  }
)

# Test getPopFrags with an incorrectly formatted region string
test_that(
  str_interp("getPopFrags throws an error for an incorrectly formatted region string"),
  {
    expect_error(
      capture.output(
        popFrags <- getPopFrags(
          testProj,
          metaColumn = metaColumn,
          numCores = 1,
          region = "chr2_1-187350807",
          NormMethod = "Raw",
          sampleSpecific = FALSE
        ),
        type = "message"
      ),
      "Invalid region input."
    )
    # A character vector or region strings will work, but a list will not.
    expect_error(
      capture.output(
        popFrags <- getPopFrags(
          testProj,
          metaColumn = metaColumn,
          numCores = 1,
          region = list("chr1:1-2", "chr2:1-2"),
          NormMethod = "Raw",
          sampleSpecific = FALSE
        ),
        type = "message"
      ),
      "Invalid region input."
    )
  }
)


# Test getPopFrags with a metaColumn not in cellColData
test_that(
  str_interp("getPopFrags throws an error for a metaColumn not in cellColData"),
  {
    expect_error(
      capture.output(
        popFrags <- getPopFrags(
          testProj,
          metaColumn = "IDoNotExist",
          numCores = 1,
          sampleSpecific = FALSE
        ),
        type = "message"
      ),
      "does not exist in the cellColData of your ArchRProj"
    )
  }
)

# Test getPopFrags with nonexistent cell subsets
test_that(
  str_interp("getPopFrags throws an error for missing cellSubsets"),
  {
    expect_error(
      capture.output(
        popFrags <- getPopFrags(
          testProj,
          metaColumn = "Clusters",
          cellSubsets = c("C1", "IDoNotExist", "C3", "C5"),
          numCores = 1,
          sampleSpecific = FALSE
        ),
        type = "message"
      ),
      "cellSubsets with NA cell counts: IDoNotExist"
    )
  }
)
