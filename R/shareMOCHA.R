#' @title \code{packMOCHA}
#'
#' @description \code{packMOCHA} combines a MOCHA object (Sample-Tile
#'   Matrix or tileResults) with its saved coverage tracks into a single zip
#'   archive. This allows MOCHA objects and the necessary coverage files for
#'   plotting to be shared to other file systems. See also:
#'   \link[MOCHA]{unpackMOCHA}
#'
#' @param MOCHAObj A MultiAssayExperiment or RangedSummarizedExperiment, from
#'   MOCHA
#' @param zipFile Filename and path of the zip archive.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#'
#' @zipFile Path to zip archive.
#'
#' @export
#'
packMOCHA <- function(MOCHAObj,
                      zipFile,
                      verbose = FALSE) {
  if (class(zipFile) != "character" || !grepl("\\.zip$", zipFile)) {
    stop("`zipFile` must be a character string ending in .zip")
  }

  saveRDS(MOCHAObj, paste(MOCHAObj@metadata$Directory, "/MOCHA_Object.RDS", sep = ""))

  zip::zipr(zipFile, MOCHAObj@metadata$Directory)
  if (verbose) {
    message("Object & Coverage files saved: ", zipFile)
  }

  return(zipFile)
}


#' @title \code{unpackMOCHA}
#'
#' @description \code{unpackMOCHA} will unpack a zip archive created by
#'   \link[MOCHA]{unpackMOCHA}, setting the stored MOCHA object's stored
#'   directory path to the new location. See also: \link[MOCHA]{packMOCHA}
#'
#' @param zipFile Filepath to the packed MOCHA object.
#' @param outDir The path to the external directory where you want to unpack the MOCHA object.
#'
#' @return MOCHAObj the MOCHA object (tileResults or Sample-Tile Matrix)
#'
#' @export
#'
unpackMOCHA <- function(zipFile,
                        outDir,
                        verbose = FALSE) {
  if (!grepl("\\.zip$", zipFile)) {
    stop("`zipFile` must be a character string ending in .zip")
  }

  zip::unzip(zipfile = zipFile, outDir = outDir)

  # Extract directory name
  newDirectory <- sub("\\/$", "", zip::zip_list(zipFile)[1, 1], )

  # load MOCHA Object into memory
  MOCHAObj <- readRDS(paste(outDir, newDirectory, "MOCHA_Object.RDS", sep = "/"))

  # Change directory path for the MOCHA object to the new directory
  MOCHAObj@metadata$Directory <- paste(outDir, newDirectory, sep = "/")

  # Over-write object so that it has the correct directory paths saved.
  saveRDS(MOCHAObj, paste(outDir, newDirectory, "MOCHA_Object.RDS", sep = "/"))

  return(MOCHAObj)
}
