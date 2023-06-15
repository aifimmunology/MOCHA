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

    if(class(zipName) != 'character'){
        stop('Please provide a character string ending in .zip for the variable zipName.')
    }else if(!grepl('\\.zip$',zipName)){
        stop('zipName does not include .zip. Please provide a character string ending in .zip for the variable zipName.')
    }

    saveRDS(MOCHA_Obj, paste(MOCHA_Obj@metadata$Directory, '/MOCHA_Object.RDS',sep=''))

    zip::zipr(zipName, MOCHA_Obj@metadata$Directory)

    return(paste('Object & Coverage files saved: ', getwd(), zipName, sep ="/"))
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

unpackMOCHA <- function(packedObject = NULL, exdir = getwd()){

    if(is.null(packedObject)){
        stop('Please provide the name of the zipped MOCHA object. Should be in the format Name.zip')
    }else if(!grepl('\\.zip$',packedObject)){
        stop('File name does not include .zip. Please provide the name of the zipped MOCHA object. Should be in the format Name.zip')
    }

    zip::unzip(zipfile = packedObject, exdir = exdir)

    #Extract directory name
    newDirectory <- sub('\\/$','',zip::zip_list(packedObject)[1,1], )

    #load MOCHA Object into memory
    MOCHA_Obj = readRDS(paste(exdir, newDirectory, 'MOCHA_Object.RDS',sep='/'))

    #Change directory path for the MOCHA object to the new directory
    MOCHA_Obj@metadata$Directory = paste(exdir, newDirectory, sep='/')

    #Over-write object so that it has the correct directory paths saved.
    saveRDS(MOCHA_Obj, paste(exdir, newDirectory, 'MOCHA_Object.RDS',sep='/'))

    return(MOCHA_Obj)
}