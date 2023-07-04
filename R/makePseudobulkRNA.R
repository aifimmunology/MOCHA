#' @title \code{makePseudobulkRNA}
#'
#' @description \code{makePseuduoulkRNA} pseodubulks a Seurat object by sample and cell type into a SummarizedExperiment, similar to MOCHA

#'
#' @param SO Seurat Object 
#' @param cellTypeColumn The column of the Seurat object with cell type information
#' @param Sample The column of the Seurat object with sample information
#' @param normalizeCounts Boolean, to determine whether to return normalized (via DESeq2) or raw counts in a SummarizedExperiment-type object. Default is TRUE and returns normalized counts.
#' @param filterByRowData A boolean flag that determines whether to export a summarizedExperiment with Genomic location for each gene (and thus filter the object to genes that match the database)
#' @param TxDb Transcript database to be used for identifying gene locations. 
#' @param OrgDb Organism database to match up gene names and locations.
#' @return A SummarizedExperiment carrying pseudobulked average expression per 1000 cells for each cell type. 
#'
#'
#' @export
makePseudobulkRNA <- function(SO, cellTypeColumn, sampleColumn = "Sample", 
                                 cellPopulations = 'All',
                                 addGenomicLocations = TRUE,
                                 Seurat_format = 'SYMBOL',
                                 TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene", 
                                 OrgDb = "org.Hs.eg.db"){

    if(any(!c(cellTypeColumn, sampleColumn) %in% colnames(SO@meta.data))){

        stop('Sample or cell type columns are missing from Seurat object. Please verify that the provided column names are correct.')

    }

    #Generate sample and cell type column in metadata for pseudobulking (after filtering down)
    metadata <- as.data.frame(SO@meta.data)
    counts <- SO@assays$RNA@counts
    metadata$CellTypeColumn = factor(unlist(metadata[,cellTypeColumn]))

    if(tolower(cellPopulations) == 'all'){
        cellPopulations <- unique(metadata$CellTypeColumn)
        metadata <- dplyr::filter(metadata, !!as.name(cellTypeColumn) %in% cellPopulations)
        counts <- counts[,rownames(metadata)]
    }

    cellTypeList <- unique(metadata$CellTypeColumn)
    sampleList <- unique(SO@meta.data[,sampleColumn])
    emptySample = counts[,1]
    emptySample[TRUE] = 0
    counts_ls <- lapply(1:length(cellTypeList), function(i) {
        
        do.call('cbind',lapply(sampleList, function(z){
            subSamples <- rownames(dplyr::filter(metadata, 
                            CellTypeColumn == cellTypeList[i] & !!as.name(sampleColumn) == z))
            if(length(subSamples) > 0){
                unlist(rowSums(as.matrix(counts[,subSamples])))
            }else{
                emptySample
            }

        }))
    })

    names(counts_ls) <- cellTypeList

    ##Sort sample-level metadata
    cellColDataNoNA <- BiocGenerics::Filter(function(x) {
        !all(is.na(x))
    }, SO@meta.data)   
    cellColDT <- data.table::as.data.table(cellColDataNoNA)
    BoolDT <- cellColDT[, lapply(.SD, function(x) {
             length(unique(x)) == 1
    }), by = c(sampleColumn)]  

    trueCols <- apply(BoolDT[,-1], 2, all)
    trueCols[[sampleColumn]] <- TRUE
    cellColDF <- as.data.frame(cellColDT)

    sampleData <- dplyr::distinct(cellColDF[, names(which(trueCols)), drop = F])

    # Set sampleIDs as rownames
    rownames(sampleData) <- sampleData[[sampleColumn]]

    #Get cell numbers
    cellCounts <- as.data.frame(table(metadata[, "Sample"], metadata[, 'CellTypeColumn']))
    names(cellCounts) <- c("Sample", "CellPop", "CellCount")
    cellCounts <- tidyr::pivot_wider(
        cellCounts,
        id_cols = "CellPop",
        names_from = "Sample",
        values_from = "CellCount"
    )
    allCellCounts <- as.data.frame(cellCounts[, -1])
    rownames(allCellCounts) <- cellCounts$CellPop

    #Process the meta data and extract the total counts and features for each population, as well as cell counts. 
    cellColDataCopy <- data.frame(metadata)

    cellColDataCopy[] <- lapply(cellColDataCopy, function(x) {
        utils::type.convert(as.character(x), as.is = TRUE)
    })

    # Assume all numeric columns are to be saved as additionalCellData
    isNumericCol <- unlist(lapply(cellColDataCopy, function(x) is.numeric(x)))
    additionalCellData <- colnames(cellColDataCopy)[isNumericCol]

    # Group by Sample (rows) and cellPop (columns)
    summarizedData <- dplyr::group_by(metadata, CellTypeColumn, !!as.name(sampleColumn))
    if (!is.null(additionalCellData)) {
        
        additionalMetaData <- lapply(additionalCellData, function(x) {

        suppressMessages(
            summarizedData2 <- dplyr::summarize(summarizedData, meanValues = mean(!!as.name(x), .groups = "drop")) 
        )
       
        summarizedData2 <- tidyr::pivot_wider(summarizedData2,
                id_cols = CellTypeColumn,
                names_from = !!as.name(sampleColumn),
                values_from = meanValues
            )

        summarizedData2 <- as.data.frame(summarizedData2)
        rownames(summarizedData2) <- summarizedData2[['CellTypeColumn']]
        summarizedData2 <- summarizedData2[, -1, drop = FALSE]

        # Filter to specific cellPopulations
        summarizedData2 <- summarizedData2[
            rownames(summarizedData2) %in% cellPopulations, ,
            drop = FALSE
        ]

        summarizedData2
        })
        names(additionalMetaData) <- additionalCellData

    }else if (is.null(additionalCellData)) {

        additionalMetaData <- NULL

    }
    remove(cellColDataCopy)

    summarizedData <- SummarizedExperiment::SummarizedExperiment(
            append(
            list(
                "CellCounts" = allCellCounts
            ),
            additionalMetaData
            ),
            colData = sampleData
        )

    rnaSE <-  SummarizedExperiment::SummarizedExperiment(counts_ls, colData = sampleData,
         metadata = list('summarizedData' = summarizedData,
            History = paste("makePseudobulkRNA", utils::packageVersion("MOCHA"))))

    #Decide if you are interweaving this data with genomic databases/locations. 
    if(addGenomicLocations){
        
        #Pull in Transcript and Organism databases. 
        TxDb <- MOCHA:::getAnnotationDbFromInstalledPkgname(dbName=TxDb, type="TxDb")
        OrgDb <- MOCHA:::getAnnotationDbFromInstalledPkgname(dbName=OrgDb, type="OrgDb")

        #Extact the GenomicRanges for all genes, while filtering out non-standard chromosomes that mess things up. 
        txList <-suppressWarnings(GenomicFeatures::genes(TxDb, single.strand.genes.only = TRUE))
        txList <- GenomeInfoDb::keepStandardChromosomes(sort(txList), species='Homo_sapiens',
                                      pruning.mode = 'coarse')
        #Reduce ranges merges transcripts for the same gene into one, so that we can one ball-park stop and end. 
        txList <- plyranges::reduce_ranges(plyranges::group_by(txList, gene_id))
        txList$GeneSymbol <- suppressWarnings(AnnotationDbi::mapIds(OrgDb, as.character(txList$gene_id), Seurat_format, "ENTREZID"))
        txList <- plyranges::filter(txList, !is.na(GeneSymbol) & GeneSymbol %in% rownames(rnaSE))
        txList <- plyranges::reduce_ranges(plyranges::group_by(txList, GeneSymbol))
        names(txList) <- txList$GeneSymbol

        #Subset down the read count matrix to just the transcripts we ahve in this database.
        newSE <- rnaSE[rownames(rnaSE)  %in% txList$GeneSymbol, ]
        newSE <- rnaSE[match(txList$GeneSymbol,rownames(newSE)), ]

        # Repackage the gene-sample matrices, the sampleData, transcript Genomic Ranges, and associated metadata (countInfo) into one SummarizedExperiment object. 
        rowRanges(newSE) = txList
       
    }else{
        newSE <- rnaSE

    }
    return(newSE)
}



#' @title \code{normalizePseudobulk }
#'
#' @description \code{normalizePseudobulk} Takes the output of makePseudobulkRNA and normalizes it. 

#'
#' @param rnaSE  the output of makePseudobulkRNA a
#' @param sampleColumn The column of the Seurat object with sample information
#' @return A SummarizedExperiment with normalized average expression
#'
#'
#' @export

normalizePseudobulk <- function(rnaSE, sampleColumn = 'Sample'){

    allMat <- lapply(names(SummarizedExperiment::assays(rnaSE)), function(x){
                old_mat <- SummarizedExperiment::assays(rnaSE)[[x]]
                dds <- DESeqDataSetFromMatrix(countData = old_mat[,colSums(old_mat) > 0],
                              colData = SummarizedExperiment::colData(rnaSE)[colSums(old_mat) > 0,],
                                design = as.formula(paste( '~', sampleColumn)))

                dds <- estimateSizeFactors(dds)
                new_mat <- DESeq2::counts(dds, normalize = TRUE)
                if(any(colSums(old_mat) == 0)){
                    filled_data <- do.call('cbind', lapply(which(colSums(old_mat) == 0), function(x){
                            rep(0, dim(new_mat)[1])
                    }))
                    new_mat <- cbind(new_mat, filled_data)
                }
                new_mat[,sort(colnames(new_mat))]
        })
    names(allMat) <- names(SummarizedExperiment::assays(rnaSE))
    se <- SummarizedExperiment::SummarizedExperiment(allMat, 
                        colData = SummarizedExperiment::colData(rnaSE),
                        rowRanges = SummarizedExperiment::rowRanges(rnaSE),
                        metadata = rnaSE@metadata)

    se@metadata$History <- append(se@metadata$History, paste("normalizePseudobulk", utils::packageVersion("MOCHA")))
    return(se)
    
    }

#' @title \code{subsetPseudobulkRNA }
#'
#' @description \code{subsetPseudobulkRNA} Takes the output of makePseudobulkRNA or normalizedPseudobulk

#'
#' @param rnaSE  the output of makePseudobulkRNA or normalizePseudobulk
#' @param sampleColumn The column of the Seurat object with sample information
#' @return A SummarizedExperiment with normalized average expression
#'
#'
#' @export
#' 
#

subsetPseudobulk <- function(rnaSE,
                              subsetBy,
                              groupList,
                              verbose = FALSE) {

                      
    if (!grepl("SummarizedExperiment",class(rnaSE)[1] )) {
        stop('The variable rnaSE that you provided is not a Summarized Experiment.')
    }

    sampleData <- SummarizedExperiment::colData(rnaSE)


    if (grepl("celltype", tolower(subsetBy))) {
        if(all(groupList %in% names(SummarizedExperiment::assays(rnaSE)))){
            newSE <- rnaSE
            SummarizedExperiment::assays(newSE) <- SummarizedExperiment::assays(newSE)[groupList]
        }else{
            stop('You asked to subset by cell type, but the provided list of cell types (groupList) includes some or all that are not found within rnaSE.')
        }

    }else if (subsetBy %in% colnames(sampleData)) {
        if(all(groupList) %in% unique(sampleData[,subsetBy])){
            newSE <- newSE[,sampleData[[subsetBy]] %in% groupList]
            newSE@metadata$summarizedData[,]
        }else{
            stop('Some indices of groupList were not found within the metadata of column provided to subsetBy.')
        }
    }
    return(newSE)
}

#' @title \code{combinePseudobulkRNA }
#'
#' @description \code{combinePseudobulkRNA} Takes the output of combinePseudobulkRNA and normalizes it. 

#'
#' @param rnaSE  the output of makePseudobulkRNA or normalizePseudobulk
#' @param cellPopulations Cell populations to keep through the combining process. 
#' @return A SummarizedExperiment with one assay named 'counts' which contains all Sample and Cell type gene expression data.
#'              It also moves all the summarizedData into the colData slot of the SummarizedExperiment.
#'
#' @export
#' 
#

combinePseudobulk <- function(rnaSE,
                              cellPopulations = 'All') {

    if (!grepl("SummarizedExperiment",class(rnaSE)[1] )) {
        stop('The variable rnaSE that you provided is not a Summarized Experiment.')
    }

    if(all(tolower(cellPopulations) == 'all')){
        subObject = rnaSE
    }else{
        subObject <- subsetPseudobulk(rnaSE, subsetBy = 'celltype', groupList = cellPopulations)
    }


    summarizedData <- S4Vectors::metadata(subObject)$summarizedData
    sampleData <- SummarizedExperiment::colData(subObject)

    cellNames <- SummarizedExperiment::assays(subObject)
    # Let's generate a new assay, that will contain the
    # the intensity for a given cell, as well as the
    # median intensity per sample-tile for all other cell types (i.e. the background)

    newAssays <- list(do.call("cbind", methods::as(cellNames, "list")))
    newSamplesNames <- unlist(lapply(names(cellNames), function(x) {
        gsub(" ", "_", paste(x, colnames(rnaSE), sep = "__"))
    }))

    names(newAssays) <- "counts"
    colnames(newAssays[[1]]) <- newSamplesNames

    # Combine Sample and Cell type for the new columns.
    # This takes the colData and repeats it across cell types for a given sample
    # Rows are now CellType__Sample. So for example 'CD16 Mono' and 'Sample1' becomes "CD16_Mono__Sample1"
    allSampleData <- as.data.frame(do.call("rbind", lapply(names(cellNames), function(x) {
        tmp_meta <- sampleData
        tmp_meta$Sample <- gsub(" ", "_", paste(x, tmp_meta$Sample, sep = "__"))
        tmp_meta$CellType <- rep(x, dim(tmp_meta)[1])
        rownames(tmp_meta) <- tmp_meta$Sample
        tmp_meta
    })))

    summarizedData <- S4Vectors::metadata(rnaSE)$summarizedData

    addMetaData <- as.data.frame(do.call('cbind',lapply(SummarizedExperiment::assays(summarizedData), function(x){
        tmpMeta <- x
        tmpMeta$CellType = rownames(x)
        tmpMeta2 <- tidyr::pivot_longer(
            tmpMeta,
            cols = colnames(x),
            names_to = "Sample",
            values_to = "Freq"
        )
        tmpMeta2 <- dplyr::mutate(
            tmpMeta2,
            Sample = gsub(" ", "_", paste(CellType, Sample, sep = "__"))
        )
        output <- unlist(tmpMeta2$Freq )
        names(output) <- tmpMeta2$Sample
        output
    })))
    addMetaData$Sample = rownames(addMetaData)
    fullMeta <- dplyr::full_join( as.data.frame(allSampleData), as.data.frame(addMetaData), by = 'Sample')

    History <- append(rnaSE@metadata$History, paste("combinePseudobulk", utils::packageVersion("MOCHA")))

    newSE <- SummarizedExperiment::SummarizedExperiment(newAssays, 
                    colData = fullMeta,
                    rowRanges = SummarizedExperiment::rowRanges(rnaSE),
                    metadata = list('History' = History))
    return(newSE)

}