#' @title \code{GeneTileAssociations}
#'
#' @description \code{GeneTileAssociations} leverages zero-inflated generalized linear mixed effect modeling to associate accessibility
#                       and gene expression.
#' @param rnaSE SummarizedExperiment object from getPseudobulkRNA
#' @param atacSE A SampleTileMatrix object from getSampleTileMatrix or combinedSampleTileMatrix
#' @param cellPopulation - name of the assay (cell population) to extract and use from rnaSE and atacSE
#' @param sampleColumn - A string: the name of the column with Sample names from the metadata. 
#'                  Must be the same from atacSE and rnaSE. This is used for aligning atacSE and rnaSE
#' @param continuousFormula : Formula used for modeling. Must be in the form exp1 ~ exp2 + other factors + (1|RandomEffect). 
#'          Exp1 will be the dummy veriable for accessibility, while exp2 will be the dummy variable for gene expression. 
#' @param ziFormula list of  to tests from atacSE to test. 
#' @param zi_threshold Zero-inflated threshold ( range = 0-1), representing the fraction of samples with zeros. At or above this threshold, the zero-inflated modeling kicks in.
#'            We don't recommend changing this unless you have fully tested out the model and understand what it is doing. 
#' @param geneList a list of genes from rnaSE to test. 
#' @param tileList optional list of tiles to test. Default is NULL, which will test all tiles within the provided distance of the TSS. 
#' @param distance A number limiting the distance between tiles and genes TSS at which point associations will not be tested. 
#' @param initialSampling An integer. Default is 5. Represents the number of initial models to test when generating a null set (which is used when models fails). If your model fails often, increase this number or choose a better formula.
#' @param numCores Optional, the number of cores to use with multiprocessing. Default is 1.
#' 
#' @return A summarized experiment summarizing the output. 
#'
#' @details
#'
#'
#' @export

GeneTile_Associations <- function(atacSE, rnaSE, cellPopulation, sampleColumn, 
                                    continuousFormula, 
                                    ziFormula, 
                                    zi_threshold = 0,
                                    geneList, 
                                    tileList = NULL,
                                    distance = 10^6,
                                    initialSampling = 5, numCores = 2) {
  
 
  if(!sampleColumn %in% colnames(SummarizedExperiment::colData(atacSE)) |
          !sampleColumn %in% colnames(SummarizedExperiment::colData(rnaSE))){

    stop('sampleColumn missing from atacSE and/or rnaSE.')
  }

  #Check whether samples align.     
  if(!all(atacSE@colData[,sampleColumn] %in% rnaSE@colData[,sampleColumn])) {
      stop('sampleColumns are not the same. Please ensure that sample names match in atacSE and rnaSE')
  }else if(!all(atacSE@colData[,sampleColumn] == rnaSE@colData[,sampleColumn])){
      warning('Reording rnaSE sample names to match atacSE.')
      rnaSE <- rnaSE[,match(colnames(atacSE), colnames(rnaSE))]
  }

  if(any(c(colnames(SummarizedExperiment::colData(atacSE)), 
      colnames(SummarizedExperiment::colData(rnaSE))) %in% c('exp1','exp2'))){

    stop('metadata of atacSE and/or rnaSE contains a column that contains the name exp1 or exp2.',
    'exp1 and exp2 are hardcoded to represent the data from atacSE and rnaSE, not the metadata.
    Please remove these and try again.')

  }

   #Test whether tileList and geneList are found in their assays.
  if(!all(geneList %in% rownames(rnaSE))){
     stop('geneList not found within rnaSE. Please read check that your genes of interest are present in rnaSE.')
  }
    
  if (class(continuousFormula) == 'character') {
    continuousFormula <- as.formula(continuousFormula)
  }

  if (class(ziFormula) == 'character') {
    ziFormula <- as.formula(ziFormula)
  }

  #Subset down to one cell type
  if (length(cellPopulation) > 1) {
  stop(
    "More than one cell population was provided. ",
    "cellPopulation must be length 1. To run over multiple cell types, ",
    "run combineSampleTileMatrix() to produce a new combined TSAM_Object and set ",
    "cellPopulation = 'counts'."
  )
  } else if (
    (!cellPopulation %in% names(SummarizedExperiment::assays(atacSE)) |
      !cellPopulation %in% names(SummarizedExperiment::assays(rnaSE)))
  ) {
    stop("cellPopulation was not found within atacSE and/or rnaSE.")
  } else if(cellPopulation == 'counts'){
    newatacSE <- atacSE
    newrnaSE <- rnaSE
    
  }else{
    newatacSE <- combineSampleTileMatrix(subsetMOCHAObject(atacSE, subsetBy = 'celltype', groupList = cellPopulation, subsetPeaks = TRUE))
    newrnaSE <- rnaSE
    SummarizedExperiment::assays(newrnaSE) <- SummarizedExperiment::assays(newrnaSE)[cellPopulation]
  }

  if(!all(tileList %in% rownames(newatacSE)) & !is.null(tileList)){
     stop('tileList not found within atacSE.  Please read check that your tiles of interest are present in atacSE.')
  }else if(is.null(tileList)){
    tileList = rownames(newatacSE)
  }


  #Test whether the continuousFormula is in the right format
  if(!(all.vars(as.formula(continuousFormula))[1] =='exp1' & any(all.vars(as.formula(continuousFormula)) %in% 'exp2'))){
    stop('continuousFormula is not in the correct format. Data from rnaSE will be used to predict atacSE, because atacSE is more complicated to model. Please format continuousFormula as exp1 ~ exp2 + OtherFixedEffects + RandomEffect.')
  }

  if(!( all(all.vars(as.formula(ziFormula))[1] %in% 'exp1'))){
    stop('ziFormula is not in the correct format. Data from rnaSE will be used to predict atacSE, because atacSE is more complicated to model. Please format ziFormula as  ~ Factor1 + Factor2 + exp2. ziFormula is optional, but should ideally include at least FragNumber.')
  }


  if (zi_threshold < 0 | zi_threshold > 1 | ! is.numeric(zi_threshold)) {
    stop("zi_threshold must be between 0 and 1.")
  }

  #Subset down each SummarizedExperiment to the rows of interest, and assays of interest. 
  newatacSE <- newatacSE[rownames(newatacSE)  %in% tileList, ]
  newrnaSE <- newrnaSE[rownames(newrnaSE)  %in% geneList, ]

  # initialize parallel processing
  cl <- parallel::makeCluster(numCores)

  #Find all combinations of tileList and geneList to test. 
  allCombo <- as.data.frame(expand.grid(as.character(geneList),as.character(tileList)), stringsAsFactors = FALSE)
  browser()

  geneTileAssociations <- .multiModalModeling(newatacSE, newrnaSE, allCombo = allCombo, continuousFormula =  continuousFormula,
                        ziFormula = ziFormula, zi_threshold = zi_threshold, initialSampling = initialSampling, family = stats::gaussian(),
                        modality = 'GeneTile',
                        numCores = numCores)

  return(geneTileAssociations)

}


#' @title \code{scATAC_Associations}
#'
#' @description \code{scATAC_Associations} leverages zero-inflated generalized linear mixed effect modeling to associate accessibility
#                       and gene expression.
#' @param atacSE SummarizedExperiment object from getPseudobulkRNA
#' @param cellPopulation - name of the assay (cell population) to extract and use from atacSE
#' @param generalSE A SummarizedExperiment containing data from another modality, like ChromVAR, 
#' @param generalAssay Name of the assay from generalSE that you wish to use for modeling.
#' @param sampleColumn - A string: the name of the column with Sample names from the metadata. 
#'                  Must be the same from atacSE and generalSE. This is used for aligning atacSE and generalSE.
#' @param continuousFormula : Formula used for modeling. Must be in the form exp1 ~ exp2 + other factors + (1|RandomEffect). 
#'          Exp1 will be the dummy veriable for accessibility, while exp2 will be the dummy variable for gene expression. 
#' @param ziFormula list of  to tests from atacSE to test. 
#' @param zi_threshold Zero-inflated threshold ( range = 0-1), representing the fraction of samples with zeros. At or above this threshold, the zero-inflated modeling kicks in.
#'            We don't recommend changing this unless you have fully tested out the model and understand what it is doing. 
#' @param tileList a list of tiles to test. 
#' @param generalList A list of which rownames of generalAssay should be used for modeling associations. 
#' @param distance A number limiting the distance between tiles and genes TSS at which point associations will not be tested. 
#' @param initialSampling An integer. Default is 5. Represents the number of initial models to test when generating a null set (which is used when models fails). If your model fails often, increase this number or choose a better formula.
#' @param numCores Optional, the number of cores to use with multiprocessing. Default is 1.
#' 
#' @return A summarized experiment summarizing the output. 
#'
#' @details
#'
#'
#' @export

scATAC_Associations <- function(atacSE, cellPopulation, 
                                    generalSE, generalAssay,
                                    sampleColumn, 
                                    continuousFormula, 
                                    ziFormula, 
                                    zi_threshold = 0,
                                    tileList,
                                    generalList, 
                                    initialSampling = 5, numCores = 2) {
  
 
  if(!sampleColumn %in% colnames(SummarizedExperiment::colData(atacSE)) |
          !sampleColumn %in% colnames(SummarizedExperiment::colData(generalSE))){

    stop('sampleColumn missing from atacSE and/or generalSE.')
  }

  #Check whether samples align.     
  if(!all(atacSE@colData[,sampleColumn] %in% generalSE@colData[,sampleColumn])) {
      stop('sampleColumns are not the same. Please ensure that sample names match in atacSE and generalSE')
  }else if(!all(atacSE@colData[,sampleColumn] == generalSE@colData[,sampleColumn])){
      warning('Reording generalSE sample names to match atacSE.')
      generalSE <- generalSE[,match(colnames(atacSE), colnames(generalSE))]
  }

  if(any(c(colnames(SummarizedExperiment::colData(atacSE)), 
      colnames(SummarizedExperiment::colData(generalSE))) %in% c('exp1','exp2'))){

    stop('metadata of atacSE and/or generalSE contains a column that contains the name exp1 or exp2.',
    'exp1 and exp2 are hardcoded to represent the data from atacSE and generalSE, not the metadata.
    Please remove these and try again.')

  }

   #Test whether tileList and geneList are found in their assays.
  if(!all(geneList %in% rownames(generalSE))){
     stop('geneList not found within generalSE. Please read check that your genes of interest are present in generalSE.')
  }
    
  if (class(continuousFormula) == 'character') {
    continuousFormula <- as.formula(continuousFormula)
  }

  if (class(ziFormula) == 'character') {
    ziFormula <- as.formula(ziFormula)
  }

  #Subset down to one cell type
  if (length(cellPopulation) > 1) {
  stop(
    "More than one cell population was provided. ",
    "cellPopulation must be length 1. To run over multiple cell types, ",
    "run combineSampleTileMatrix() to produce a new combined TSAM_Object and set ",
    "cellPopulation = 'counts'."
  )
  } else if (
    (!cellPopulation %in% names(SummarizedExperiment::assays(atacSE)) |
      !cellPopulation %in% names(SummarizedExperiment::assays(generalSE)))
  ) {
    stop("cellPopulation was not found within atacSE and/or generalSE.")
  } else if(cellPopulation == 'counts'){
    newatacSE <- atacSE
    
  }else{
    newatacSE <- combineSampleTileMatrix(subsetMOCHAObject(atacSE, subsetBy = 'celltype', groupList = cellPopulation, subsetPeaks = TRUE))
   
  }
  newGeneral <- generalSE
  SummarizedExperiment::assays(newGeneral) <- SummarizedExperiment::assays(newGeneral)[generalAssay]

  if(!all(tileList %in% rownames(newatacSE))){
     stop('tileList not found within atacSE.  Please read check that your tiles of interest are present in atacSE.')
  }


  #Test whether the continuousFormula is in the right format
  if(!(all.vars(as.formula(continuousFormula))[1] =='exp1' & any(all.vars(as.formula(continuousFormula)) %in% 'exp2'))){
    stop('continuousFormula is not in the correct format. Data from generalSE will be used to predict atacSE, because atacSE is more complicated to model. Please format continuousFormula as exp1 ~ exp2 + OtherFixedEffects + RandomEffect.')
  }

  if(!( all(all.vars(as.formula(ziFormula))[1] %in% 'exp1'))){
    stop('ziFormula is not in the correct format. Data from generalSE will be used to predict atacSE, because atacSE is more complicated to model. Please format ziFormula as  ~ Factor1 + Factor2 + exp2. ziFormula is optional, but should ideally include at least FragNumber.')
  }


  if (zi_threshold < 0 | zi_threshold > 1 | ! is.numeric(zi_threshold)) {
    stop("zi_threshold must be between 0 and 1.")
  }

  #Subset down each SummarizedExperiment to the rows of interest, and assays of interest. 
  newatacSE <- newatacSE[rownames(newatacSE)  %in% tileList, ]
  newGeneral <- newGeneral[rownames(newGeneral)  %in% generalList, ]

  # initialize parallel processing
  cl <- parallel::makeCluster(numCores)

  #Find all combinations of tileList and geneList to test. 
  allCombo <- as.data.frame(expand.grid(as.character(generalList),as.character(tileList)), stringsAsFactors = FALSE)

  geneTileAssociations <- .multiModalModeling(newatacSE, newGeneral, allCombo = allCombo, continuousFormula =  continuousFormula,
                        ziFormula = ziFormula, zi_threshold = zi_threshold, initialSampling = initialSampling, family = stats::gaussian(), modality = 'GeneTile',
                        numCores = numCores)

  return(geneTileAssociations)

}



#' @title \code{scRNA_Associations}
#'
#' @description \code{scRNA_Associations} leverages generalized linear mixed effect modeling to associate gene expression to another modality. 
#' @param rnaSE SummarizedExperiment object from getPseudobulkRNA
#' @param cellPopulation - name of the assay (cell population) to extract and use from rnaSE and atacSE
#' @param generalSE A SummarizedExperiment containing data from another modality, like ChromVAR, 
#' @param generalAssay Name of the assay from generalSE that you wish to use for modeling.
#' @param sampleColumn - A string: the name of the column with Sample names from the metadata. 
#'                  Must be the same from atacSE and rnaSE. This is used for aligning atacSE and rnaSE
#' @param formula : Formula used for modeling. Must be in the form exp1 ~ exp2 + other factors + (1|RandomEffect), where 
#'          exp1 is the dummy variable for gene expression and exp2 is the dummy variable for the rownames from generalSE. 
#' @param geneList a list of genes from rnaSE to test. 
#' @param generalList A list of which rownames of generalAssay should be used for modeling associations. 
#' @param initialSampling An integer. Default is 5. Represents the number of initial models to test when generating a null set (which is used when models fails). If your model fails often, increase this number or choose a better formula.
#' @param numCores Optional, the number of cores to use with multiprocessing. Default is 1.
#' 
#' @return A summarized experiment summarizing the output. 
#'
#' @details
#'
#'
#' @export

scRNA_Associations <- function(rnaSE, cellPopulation, 
                                    generalSE, generalAssay,
                                    sampleColumn, 
                                    formula, 
                                    geneList, 
                                    generalList = NULL,
                                    initialSampling = 5, numCores = 2) {
  
 
  if(!sampleColumn %in% colnames(SummarizedExperiment::colData(generalSE)) |
          !sampleColumn %in% colnames(SummarizedExperiment::colData(rnaSE))){

    stop('sampleColumn missing from generalSE and/or rnaSE.')
  }

  #Check whether samples align.     
  if(!all(generalSE@colData[,sampleColumn] %in% rnaSE@colData[,sampleColumn])) {
      stop('sampleColumns are not the same. Please ensure that sample names match in generalSE and rnaSE')
  }else if(!all(generalSE@colData[,sampleColumn] == rnaSE@colData[,sampleColumn])){
      warning('Reording rnaSE sample names to match generalSE.')
      rnaSE <- rnaSE[,match(colnames(generalSE), colnames(rnaSE))]
  }

  if(any(c(colnames(SummarizedExperiment::colData(generalSE)), 
      colnames(SummarizedExperiment::colData(rnaSE))) %in% c('exp1','exp2'))){

    stop('metadata of generalSE and/or rnaSE contains a column that contains the name exp1 or exp2.',
    'exp1 and exp2 are hardcoded to represent the data from generalSE and rnaSE, not the metadata.
    Please remove these and try again.')

  }

   #Test whether tileList and geneList are found in their assays.
  if(!all(geneList %in% rownames(rnaSE))){
     stop('geneList not found within rnaSE. Please read check that your genes of interest are present in rnaSE.')
  }
    
  if (class(continuousFormula) == 'character') {
    continuousFormula <- as.formula(continuousFormula)
  }

  if (class(ziFormula) == 'character') {
    ziFormula <- as.formula(ziFormula)
  }

  #Subset down to one cell type
  if (length(cellPopulation) > 1) {
  stop(
    "More than one cell population was provided. ",
    "cellPopulation must be length 1. To run over multiple cell types, ",
    "run combineSampleTileMatrix() to produce a new combined TSAM_Object and set ",
    "cellPopulation = 'counts'."
  )
  } else if (
    (!cellPopulation %in% names(SummarizedExperiment::assays(generalSE)) |
      !cellPopulation %in% names(SummarizedExperiment::assays(rnaSE)))
  ) {
    stop("cellPopulation was not found within rnaSE and/or generalAssay within generalSE.")
  } else if(cellPopulation == 'counts'){
    newrnaSE <- rnaSE
    
    
  }else{
    newrnaSE <- rnaSE
    SummarizedExperiment::assays(newrnaSE) <- SummarizedExperiment::assays(newrnaSE)[cellPopulation]
  }

  newGeneral <- generalSE
  SummarizedExperiment::assays(newGeneral) <- SummarizedExperiment::assays(newGeneral)[generalAssay]

  if(!all(generalList %in% rownames(newGeneral))){
     stop('generalList not found within generalSE.  Please read check that your tiles of interest are present in generalSE.')
  }


  #Test whether the continuousFormula is in the right format
  if(!(all.vars(as.formula(formula))[1] =='exp1' & any(all.vars(as.formula(formula)) %in% 'exp2'))){
    stop('continuousFormula is not in the correct format. Data from rnaSE will be used to predict atacSE, because atacSE is more complicated to model. Please format continuousFormula as exp1 ~ exp2 + OtherFixedEffects + RandomEffect.')
  }

  #Subset down each SummarizedExperiment to the rows of interest, and assays of interest. 
  newGeneral <- newGeneral[rownames(newGeneral)  %in% generalList, ]
  newrnaSE <- newrnaSE[rownames(newrnaSE)  %in% geneList, ]

  # initialize parallel processing
  cl <- parallel::makeCluster(numCores)

  #Find all combinations of tileList and geneList to test. 
  allCombo <- as.data.frame(expand.grid(as.character(generalList),as.character(geneList)), stringsAsFactors = FALSE)
  browser()

  geneAssociations <- .multiModalModeling(newrnaSE, newGeneral, allCombo = allCombo, continuousFormula = formula,
                        ziFormula = ~0, zi_threshold = 0, initialSampling = initialSampling, family = stats::poisson(),
                        modality = 'GeneralGene',
                        numCores = numCores)

  return(geneAssociations)

}







#' @title \code{general_associations}
#'
#' @description \code{general_associations} allows for generalized linear mixed effect modeling across multiple modalities,
#'  as long as samples are aligned. This is only designed for integrating datatypes with normal distributions.
#' @param SE1 SummarizedExperiment object for measurement-type 1. 
#' @param SE2 SummarizedExperiment object for measurement-type 2. 
#' @param assay1 - name of the assay to extract and use from SE1
#' @param assay2 - name of the assay to extract and use from SE2
#' @param sampleColumn - A string: the name of the column with Sample names from the metadata. 
#'                  Must be the same from SE1 and SE2. This is used for aligning SE1 and SE2
#' @param formula : Formula used for modeling. Must be in the form exp1 ~ exp2 + other factors + (1|RandomEffect)
#' @param sig1 A list of rows from SE1 to test. 
#' @param sig2 A list of rows from SE2 to test. 
#' @param initialSampling An integer. Default is 5. Represents the number of initial models to test when generating a null set (which is used when models fails). If your model fails often, increase this number or choose a better formula.
#' @param numCores Optional, the number of cores to use with multiprocessing. Default is 1.
#' 
#' @return A summarized experiment summarizing the output. 
#'
#' @details
#'
#'
#' @export

general_associations <- function(SE1, SE2, assay1, assay2, sampleColumn, formula, sig1, sig2, family = stats::gaussian(), initialSampling = 5, numCores) {
  
   #Check whether samples align.     
  if(!all(SE1@colData[,sampleColumn] %in% SE2@colData[,sampleColumn])) {
      stop('sampleColumns are not the same. Please ensure that sample names match in SE1 and SE2')
  }else if(!all(SE1@colData[,sampleColumn] == SE2@colData[,sampleColumn])){
      warning('Reording SE1 sample names to match SE2.')
      SE2 <- SE2[,match(colnames(SE1), colnames(SE2))]
  }

  if(any(c(colnames(SummarizedExperiment::colData(SE1)), 
      colnames(SummarizedExperiment::colData(SE2))) %in% c('exp1','exp2'))){

    stop('metadata of SE1 and/or SE2 contains a column that contains the name exp1 or exp2.',
    'exp1 and exp2 are hardcoded to represent the data from SE1 and SE2, not the metadata.
    Please remove these and try again.')

  }
 
  if(!all(sig2 %in% rownames(SE2))){
     stop('sig2 not found within SE2. Please read documentation.')
  }
    
  #Test whether sig1 and sig2 are found in their assays.
  if(!all(sig1 %in% rownames(SE1))){
     stop('sig1 not found within SE1. Please read documentation.')
  }

  #Test whether the formula is in the right format
  if(!(all.vars(as.formula(formula))[1] =='exp1' & any(all.vars(as.formula(formula)) %in% 'exp2'))){
    stop('formula is not in the correct format. Data from SE2 will be used to predict SE1. Please format formula as exp1 ~ exp2 + OtherFixedEffects + RandomEffect.')
  }

  #Subset down each SummarizedExperiment to the rows of interest, and assays of interest. 
  newSE1 <- SE1[rownames(SE1)  %in% sig1, ]
  newSE2 <- SE2[rownames(SE2)  %in% sig2, ]
  SummarizedExperiment::assays(newSE1) <- SummarizedExperiment::assays(newSE1)[assay1]
  SummarizedExperiment::assays(newSE2) <- SummarizedExperiment::assays(newSE2)[assay2]

  # initialize parallel processing
  cl <- parallel::makeCluster(numCores)

  #Find all combinations of sig1 and sig2 to test. 
  allCombo <- as.data.frame(expand.grid(as.character(sig2),as.character(sig1)), stringsAsFactors = FALSE)

  generalAssociations <- .multiModalModeling(newSE1, newSE2, allCombo = allCombo, continuousFormula = formula,
                        ziFormula = ~0, zi_threshold =0, initialSampling, family = family, modality = 'GeneralAssociation',
                        numCores = numCores)

  return(generalAssociations)

}




#' @title Internal function to run generating null results, in case of model failure when associations multiple modalities
#'
#' @description \code{generateNULL_Associations} Runs null results, either for glmmTMB or lmerTest/lme4
#' @param modelingData the data used for modeling. 
#' @param metaData the metadata associated with the data
#' @param continuousFormula The formula for the continuous data that should be used within glmmTMB. It should be in the
#'   format (exp ~ factors). All factors must be found in column names of the SE_Object colData (i.e. sample metadata).
#'   See \link[glmmTMB]{glmmTMB}.
#' @param ziFormula The formula for the zero-inflated data that should be used within glmmTMB. It should be in the
#'   format ( ~ factors). All factors must be found in column names of metaData
#' @param family family for data distribution to be used in the model
#' @param initialSampling Size of data to use for pilot
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' 
#' @return output_vector A linear model
#'
#' @noRd

generateNULL_Associations <- function(mat1, mat2, metaData, allCombo, continuousFormula, ziFormula, family, initialSampling){

  #Pull out random pairs to test
  pilotIndices <- sample(x = 1:length(allCombo), size = initialSampling, replace = FALSE)

  #Get the full name (we need to reverse the order to in order to have the first index be mat1, and the second related to mat2.)
  pilotIndices2 <- allCombo[pilotIndices]

  #Generate a list of models
  modelList <- lapply(pilotIndices2, function(x) {
     df <- data.frame(exp1 = unlist(mat1[x[1],]), exp2 = unlist(mat2[x[2],]), metaData, stringsAsFactors= FALSE)

    tryCatch(
      {
        glmmTMB::glmmTMB(as.formula(continuousFormula),
            ziformula = as.formula(ziFormula),
            data = df,
            family = family,
            REML = TRUE)
       
      },
      error = function(e) {
        NA
      }
    )
  })

  names(modelList) <- unlist(lapply(pilotIndices2, function(x) { paste(x[1],x[2],sep = '_')}))

  NAList <- unlist(lapply(modelList, function(x){
                    ## First identify all the pilot models that ran at all without errors. 
                   tmp1 <- all(is.na(x))
                   if(tmp1){ 
                    return(FALSE)
                   }else{
                     ## Then identify all the models that converged. 
                    tryCatch(
                      {
                        tmp1 <- summary(x)$coefficients
                        return(TRUE)
                      }, error = function(e){ return(FALSE)})
                   }
            }))
   
  if (all(!NAList)) {
    stop("For the initial sampling, every test model failed. Reconsider modelFormula or increase initial sampling size.")
  } else {
    idx <- which(NAList)
    # Did any of the models have both zero-inflated and continous portions?
    bothZI_Cont <- unlist(lapply(idx, function(x){
      ziDF <- summary(modelList[[x]])$coefficients$zi
      !is.null(ziDF) | sum(dim(ziDF)) != 0

    }))

    if(all(!bothZI_Cont) & !ziFormula %in% c('~ 0', '~0')){
      warning('No working models using the zero-inflated formula. Do you need to modify the zero-inflated formula?')
      modelRes <- modelList[[idx[1]]]
    }else if(ziFormula %in% c('~ 0', '~0')){
      modelRes <- modelList[[idx[1]]]
    }else {
      modelRes <- modelList[[idx[which(bothZI_Cont)[1]]]]
    }

    #Extract the first representative model. And use it to create a null template. 
    coeff2 <- summary(modelRes)$coefficients
    coeff <- lapply(coeff2, as.data.frame)
    coeff$cond[!is.na(coeff$cond)] <- NA
    rownames(coeff$cond)[grepl('(Intercept)',rownames(coeff$cond))] = 'Intercept'
    if(all(!bothZI_Cont)){
      combinedCoeff <- coeff$cond
    }else{
      coeff$zi[!is.na(coeff$zi)] <- NA
      rownames(coeff$zi)[grepl('(Intercept)',rownames(coeff$zi))] = 'Intercept'
      rownames(coeff$zi) <- paste('ZI', rownames(coeff$zi), sep ='_')
      combinedCoeff <- rbind(coeff$cond, coeff$zi)
    }
    
    Resid <- stats::resid(modelRes)
    
    if(!all(metaData$Sample %in% names(Resid))){
      NA_samples <- rep(NA, sum(!metaData$Sample %in% names(Resid)))
      names(NA_samples) <- metaData$Sample[!metaData$Sample %in% names(Resid)]
      Resid <- c(Resid, NA_samples)
    }
    Resid <- Resid[match(names(Resid),metaData$Sample)]
    Resid[!is.na(Resid)] <- NA
    varCorrObj <- glmmTMB::VarCorr(modelRes)
    if(is.null(varCorrObj$cond) & all(!bothZI_Cont)){
      varCorrObj <- list('cond' = list('NoRandomEffects'= NA), 'zi' = NULL)
      residual = c('Residual' = NA)
    } else if(is.null(varCorrObj$cond)){
      varCorrObj <- list('cond' = list('NoRandomEffects'= NA), 'zi' = list('None'= NA))
      residual = c('Residual' = NA)
    }else{
      residual = as.vector(attr(varCorrObj$cond, "sc")^2)
      names(residual) = 'Residual'
    }
    cond_other = unlist(varCorrObj$cond)
    names(cond_other) = paste('Cond', names(cond_other), sep = "_")
    
   
    if(!is.null(varCorrObj$zi)){
      zi_other = unlist(varCorrObj$zi)
      names(zi_other) = paste('ZI', names(zi_other), sep = "_")
      varcor_df <- c(cond_other, zi_other,residual)
    }else{
      varcor_df <- c(cond_other, residual)
    }
    varcor_df[!is.na(varcor_df)] = NA

    nullDFList <- list('Coeff' = combinedCoeff, 'Resid' = Resid, 'VCov'= varcor_df)

    rm(modelList)
  }
  return(nullDFList)
}

#' @title \code{individualAssociations}
#'
#' @description \code{mindividualAssociations} allows for generalized linear mixed effect modeling across multiple modalities,
#'  as long as samples are aligned. This is only designed for integrating datatypes with normal distributions.
#' @param x combination to test for associations. First index is the rowname from mat1, and second index is rowname from mat2 to test.
#'            mat1, mat2, and metadata1, as well as nullDFList should all be exported to the core beforehand.  
#' 
#' @return A summarized experiment summarizing the output. 
#'
#' @details
#'
#' @noRd


individualAssociations <- function(x){
  
  df <- data.frame(exp1 = unlist(mat1[x[[1]],]), exp2 = unlist(mat2[x[[2]],]), metaData)


  output_vector <- tryCatch(
    {
     
      if(sum(df$exp1 == 0, na.rm = TRUE)/length(df$exp1) == 0){
        modelRes <- glmmTMB::glmmTMB(as.formula(continuousFormula),
          ziformula = ~ 0,
          data = df,
          family = family,
          REML = TRUE,
          control = glmmTMB::glmmTMBControl(parallel = 1)
        )
       
      }else if(sum(df$exp1 == 0, na.rm = TRUE)/length(df$exp1) <= zi_threshold){
        df$exp1[df$exp1 == 0] = NA
        modelRes <- glmmTMB::glmmTMB(as.formula(continuousFormula),
          ziformula = ~ 0,
          data = df,
          family = family,
          REML = TRUE,
          control = glmmTMB::glmmTMBControl(parallel = 1)
        )
       
      }else {
        modelRes <- glmmTMB::glmmTMB(as.formula(continuousFormula),
          ziformula = as.formula(ziFormula),
          data = df,
          family = family,
          REML = TRUE,
          control = glmmTMB::glmmTMBControl(parallel = 1)
        )
      }

      Coeff <- lapply(summary(modelRes)$coefficients, as.data.frame)
      rownames(Coeff[['cond']])[grepl('(Intercept)',rownames(Coeff[['cond']]))] = 'Intercept'
      if(sum(dim(Coeff$zi)) != 0){
        rownames(Coeff$zi)[grepl('(Intercept)',rownames(Coeff$zi))] = 'Intercept'
        rownames(Coeff$zi) <- paste('ZI', rownames(Coeff$zi), sep ='_')
      }else{
        Coeff$zi = nullDFList$Coeff[grepl('ZI',rownames(nullDFList$Coeff)),]
      }
      combinedCoeff <- rbind(Coeff$cond, Coeff$zi)


      Resid <- stats::resid(modelRes)
      
      #Clean up residuals. Data that is NA will be removed from the residuals, so we add them back in as NAs for the sake of completion. 
      if(!all(metaData$Sample %in% names(Resid))){
        NA_samples <- rep(NA, sum(!metaData$Sample %in% names(Resid)))
        names(NA_samples) <- metaData$Sample[!metaData$Sample %in% names(Resid)]
        Resid <- c(Resid, NA_samples)
      }
      Resid <- Resid[match(names(Resid),metaData$Sample)]

      #Now process the variance from random effects. 
      varcor_df <- tryCatch({
        varCorrObj <- glmmTMB::VarCorr(modelRes)
        cond_other = unlist(varCorrObj$cond)
        names(cond_other) = paste('Cond', names(cond_other), sep = "_")
        residual = as.vector(attr(varCorrObj$cond, "sc")^2)
        names(residual) = 'Residual'

        #Process variance
        if(!is.null(varCorrObj$zi)){
          zi_other = unlist(varCorrObj$zi)
          names(zi_other) = paste('ZI', names(zi_other), sep = "_")
          varcor_df <- c(cond_other, zi_other,residual)
        }else if(all(df$exp1 !=0, na.rm = TRUE)){
          subNull = nullDFList[grepl('ZI_', names(nullDFList))]
          zi_other = rep(0, length(subNull))
          names(zi_other) = names(subNull)
          varcor_df <- c(cond_other, zi_other,residual)
        }else {
          varcor_df <- c(cond_other, residual)
        }
        varcor_df
      }, error = function(e){
        nullDFList$VCov
      })
      
      list('Coeff' = combinedCoeff, 'Resid' = Resid, 'VCov'= varcor_df)
    },
    error = function(e) {
      nullDFList
    }
  )
  return(output_vector)
}




#' @title \code{multiModalModeling}
#'
#' @description \code{multiModalModeling} allows for generalized linear mixed effect modeling across multiple modalities,
#'  as long as samples are aligned. This is only designed for integrating datatypes with normal distributions.
#' @param mat1 SummarizedExperiment object for measurement-type 1. 
#' @param mat2 SummarizedExperiment object for measurement-type 2. 
#' @param allCombos all combinations of rows from mat1 and mat2 to test. 
#' @param metaData - metadata from samples
#' @param distance 
#' @param formula : Formula used for modeling. Must be in the form exp1 ~ exp2 + other factors + (1|RandomEffect)
#' @param ziformula1 : Formula used for modeling the zero-inflated component. Must be in the form  ~ exp2 + other factors + (1|RandomEffect)
#' @param family distribution family to be provided to glmmTMB
#' @param numCores Optional, the number of cores to use with multiprocessing. Default is 1.
#' 
#' @return A summarized experiment summarizing the output. 
#'
#' @details
#'
#'
#' @noRd

.multiModalModeling <- function(SE1, SE2, allCombo, continuousFormula, ziFormula, zi_threshold, initialSampling, family, modality, numCores) {
  
  ## Extract matrices and metadata
  mat1 <- SummarizedExperiment::assays(SE1)[[1]]
  mat2 <- SummarizedExperiment::assays(SE2)[[1]]

  if(any(modality %in% c('GeneTile','General_scATAC'))){
    mat1 <- log2(mat1 + 1)
  }

  metaData <- SummarizedExperiment::colData(SE1)

  if(numCores <= 1){
    stop('numCores must be greater than 1. This method is meant to be parallelized.')
  }

  if (any(c('exp1','exp2') %in% colnames(metaData))) {
    stop("Variables for exp1 and exp2 are within the metaData.")
  }

  if (c(class(continuousFormula)) != "formula") {
    stop("continuousFormula was not provided as a formula.")
  }
  if (c(class(ziFormula)) != "formula" & !is.null(ziFormula) ) {
    stop("ziFormula was not provided as a formula.")
  }

  #Transform in character strings for multithreading. 
  continuousFormula <- deparse(continuousFormula)
  ziFormula <- deparse(ziFormula)

  if (!all(all.vars(continuousFormula) %in% c("exp1", "exp2", colnames(metaData)))) {
    stop("Model formula is not in the correct format (exp1 ~ factors) or model factors are not found in column names of metadata.")
  }

  if (!all(all.vars(ziFormula) %in% c( "exp2", colnames(metaData)))) {
    stop("factors from the ziFormula were not found in the metadata.")
  }

  # convert allCombo into a list to iterate ove (need to reverse order for modeling, since expand grid tends to flip it.)
  allComboList <- lapply(c(1:dim(allCombo)[1]), function(x){
                   rev(as.character(unlist(allCombo[x,])))
      })


  #Generate the nullDFList, incase of model failure. 
  nullDFList <- generateNULL_Associations(mat1=mat1, mat2 = mat2, metaData = metaData,
                       allCombo = allComboList, continuousFormula = continuousFormula, ziFormula = ziFormula,
                       family = family, initialSampling = initialSampling)

  cl <- parallel::makeCluster(numCores)
  #export glmmTMB to each cores
  parallel::clusterEvalQ(cl, {
        library('glmmTMB')
      })

  #Export data to each core
  parallel::clusterExport(cl, c("mat1", "mat2", "metaData", "continuousFormula", "ziFormula", "nullDFList", "family",
  "zi_threshold", "individualAssociations"), envir = environment())
    
  #Test for each individual combination. 
  allRes <- pbapply::pblapply(cl = cl, X =allComboList, individualAssociations)

  # stop parallel processing
  parallel::stopCluster(cl)
  
  rownamesList <- unlist(lapply(allComboList, function(x){ paste(x[1],  x[2], sep = "_")}))

  outputList <- processModelOutputs(modelOutputList=allRes, nullDFList = nullDFList, 
                      rownamesList = rownamesList,
                                  SummarizedExperimentObj = SE1, ranged = FALSE, returnList = TRUE) 


  # Process rowData. 
  rowData1 <- as.data.frame(SummarizedExperiment::rowData(SE1))
  rowData2 <- as.data.frame(SummarizedExperiment::rowData(SE2))
  if(any(dim(rowData1) ==0) & any(dim(rowData1) == 0)){
    #No row data available
    allRowData = NULL
  }else if(any(dim(rowData1) ==0)){
    #Only rowData from SE2
    allRowData <- rowData1[allCombo$Var2,]
    allRowData$Obj2 = allCombo$Var2
    rownames(allRowData) = paste(allCombo$Var2,  allCombo$Var1, sep = "_")
      
  }else if(any(dim(rowData2) ==0)){
    #Only rowData from SE1
    allRowData <- rowData1[allCombo$Var1,]
    allRowData$Obj1 = allCombo$Var1
    rownames(allRowData) = paste(allCombo$Var2,  allCombo$Var1, sep = "_")

  }else{

    allRowData1 <- rowData1[allCombo$Var1,]
    allRowData1$Obj1 = allCombo$Var1
    rownames(allRowData1) = paste(allCombo$Var2,  allCombo$Var1, sep = "_")
    allRowData2 <- rowData2[allCombo$Var2,]
    allRowData2$Obj2 = allCombo$Var2
    rownames(allRowData2) = paste(allCombo$Var2,  allCombo$Var1, sep = "_")
    allRowData <- cbind(allRowData1, allRowData2)
  }
  residualDF <- outputList[[2]]
  colnames(residualDF)[is.na(colnames(residualDF))] = rownames(SummarizedExperiment::colData(SE1))[
      !rownames(SummarizedExperiment::colData(SE1)) %in% colnames(residualDF)]
  #Create ResidualObject, and metadata list.
  ResidualSE <- SummarizedExperiment::SummarizedExperiment(
                      list('Residual' =  residualDF),
                      colData = SummarizedExperiment::colData(SE1),
                      rowData = allRowData,
                      metadata = S4Vectors::metadata(SE1)
              )
              
  #Package up the metadata list. 
  metaDataList <- list('Residuals' = ResidualSE,
                      'RandomEffectVariance' = outputList[[3]])

  #Bind results into one data.frame.
  resDF <- SummarizedExperiment::SummarizedExperiment(outputList[[1]],
                    rowData = allRowData,
                    metadata= metaDataList)
  resDF@metadata = append(resDF@metadata, list('Type' = modality))

  # concatenate results and return
  return(resDF)
}