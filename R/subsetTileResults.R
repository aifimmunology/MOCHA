#' @title \code{subsetTileResults}
#'
#' @description \code{subsetTileResults} subsets a tileResults-type object (from callOpenTiles),
#' 										  either by cell type or sample metadata.
#' @param tileObject A MultiAssayExperiment object from callOpenTiles, 
#' @param subsetBy the variable to subset by. Can either be 'celltype', or a column from the sample metadata (see colData(tileObject) )
#' @param groupList the list of cell type names or sample-associated data that should be used to subset the tileObject 
#' @param na.rm removes groups that are NA if set to true. If set to false, then you filter for everything in the groupList and also NA values. 
#'
#' @return tileObject the input tileObject, filtered down to either the cell type or samples desired.
#'
#'
#' @export


subsetTileResults <- function(tileObject,
						  subsetBy,
                          groupList,
						  na.rm = TRUE
						  ) {
						  
		
  sampleData <- MultiAssayExperiment::colData(tileObject)
   
  if (!(subsetBy %in% colnames(sampleData)) & tolower(subsetBy) != 'celltypes') {
    stop("Error: subsetBy must either be a column name within colData(tileObjects), or 'celltype'.")
  } 
  
  if (subsetBy %in% colnames(sampleData)) {
    if(!all(groupList %in% unique(sampleData[[subsetBy]]))){ 
		stop("Error: groupList includes names not found within the object sample data. Please check groupList.")
	}
	
	
	if(na.rm){
		keep <- rownames(sampleData)[which(sampleData[[subsetBy]] %in% groupList)]
		return(tileObject[,keep])
	}else{
	
		keep <- rownames(sampleData)[which(sampleData[[subsetBy]] %in% groupList | is.na(sampleData[[subsetBy]]))]
		return(tileObject[,keep])
	}
	
  } 
  
  
  # To subset by cell type, first we have to verify that all cell type names were found within the tile object.
  # then we simply do a simple subsetting process, like you would with a list.
  if (tolower(subsetBy) == 'celltypes') {
    if(!all(groupList %in% names(tileObject))){ 
		stop("Error: groupList includes celltypes not found within tileObject.")
	}
	
	keep <- which(names(tileObject) %in% groupList)
	return(tileObject[keep])
	
  } 
 
   
}