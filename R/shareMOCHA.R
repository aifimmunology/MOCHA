#' @title \code{packMOCHA}
#'
#' @description \code{packMOCHA} compresses a MOCHA sample-tile object or tile results object, along with the coverage tracks necessary for plotting, 
#'      and exports 
#'
#' @param Object A MultiAssayExperiment or RangedSummarizedExperiment, from MOCHA
#' @param zipName zipName Name of the compressed file that will be exported. 
#'
#' @return path to zipped object. 
#'
#'
#' @export
#' 
#' 

shareMOCHA <- function(MOCHA_Obj, zipName = 'ZippedMOCHA.zip'){

    saveRDS(MOCHA_Obj, paste(MOCHA_Obj@metadata$Directory, '/MOCHA_Object.RDS',sep=''))

    zip::zip(zipName, MOCHA_Obj@metadata$Directory)

    return(paste(getwd(), zipName, sep ="/"))
}


#' @title \code{unpackMOCHA}
#'
#' @description \code{unpackMOCHA} Take a packed MOCHA object, and unpacks it locally, while resetting the directory path to the new location. 
#' 
#' @packedObject a string. The path to the packed MOCHA object.  
#' @exdir a string. The path to the external directory where you want to unpack the MOCHA object. 
#' 
#' @return the MOCHA object (tile results or Sample tile matrix)
#'
#'
#' @export
#' 
#' 

shareMOCHA <- function(packedObject = 'ZippedMOCHA.zip', exdir = getwd()){

    zip::unzip(zipfile = packedObject, exdir = exdir)

    MOCHA_Obj = readRDS(paste(exdir, '/MOCHA/MOCHA_Object.RDS',sep=''))

    MOCHA_Obj@metadata$Directory = paste(exdir, '/MOCHA',sep='')

    return(MOCHA_Obj)
}