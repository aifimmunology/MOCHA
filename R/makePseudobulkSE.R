
## work flow


makePseudobulkSE <- function(SO, cellTypeColumn, sampleColumn, 
                                 filterByRowData = TRUE,
                                 TxDb = TxDb.Hsapiens.UCSC.hg38.refGene, 
                                 OrgDb = org.Hs.eg.db){

    if(c(cellTypeColumn, sampleColumn) %in% colnames(SO@meta.data)){

        stop('Sample or cell type columns are missing from Seurat object. Please verify that the provided column names are correct.')

    }

    SO@meta.data$sample_celltype = paste(SO@meta.data[,cellTypeColumn], SO@meta.data[,sampleColumn], sep = '__')

    mat <- Seurat::AverageExpression(SO, assays = 'RNA', slot = 'counts', group.by = 'sample_celltype', return.seurat = FALSE)
    mat_df <- as.data.frame(mat$RNA)
    colnames(mat_df) <- gsub(" ","_", colnames(mat_df))
    mat_df$AllGenes <- rownames(mat_df)
    ## Clean up isoforms, if any. 
    mat_df$Genes <- gsub("\\.*","", mat_df$AllGenes) 

    summarizeMat <- dplyr::group_by(mat_df, Genes) 
    summarizeMat <- dplyr::select(summarizeMat, !AllGenes)
    summarizeMat <- dplyr::summarise_all(summarizeMat, mean)
    newSumMat <- as.data.frame(summarizeMat[,-1])
    rownames(newSumMat) <- summarizeMat$Genes 

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

    fullMeta <- dplyr::group_by_at(as.data.frame(SO@meta.data), c(sampleColumn, cellTypeColumn))

    fullMeta[,percentMito] = as.numeric(unlist(fullMeta[,percentMito] ))
    countInfo <- dplyr::summarise(fullMeta, nCount = sum(nCount_RNA), nFeature = sum(nFeature_RNA),
                                      CellCount = dplyr::n())

    if(filterByRowData){
        txList <-suppressWarnings(GenomicFeatures::genes(TxDb, single.strand.genes.only = TRUE))
        txList <- GenomeInfoDb::keepStandardChromosomes(sort(txList), species='Homo_sapiens',
                                      pruning.mode = 'coarse')
        txList <- plyranges::reduce_ranges(plyranges::group_by(txList, gene_id))
        txList$GeneSymbol <- suppressWarnings(AnnotationDbi::mapIds(OrgDb, as.character(txList$gene_id), "SYMBOL", "ENTREZID"))
        txList <- plyranges::filter(txList, !is.na(GeneSymbol) & GeneSymbol %in% rownames(newSumMat))
        txList <- plyranges::reduce_ranges(plyranges::group_by(txList, GeneSymbol))
        names(txList) <- txList$GeneSymbol

        subsetMat <- newSumMat[rownames(newSumMat) %in% txList$GeneSymbol, ]
        subsetMat <- subsetMat[match(txList$GeneSymbol,rownames(subsetMat)), ]

        cellTypes <- gsub("__.*", "", colnames(subsetMat))
        cellMatList <- lapply(unique(cellTypes), function(x){
            
            cellMat <- subsetMat[, x == cellTypes]
            colnames(cellMat) <- gsub(".*__","", colnames(cellMat))
            cellMat
            })
        names(cellMatList) = unique(cellTypes)

        rnaSE <- SummarizedExperiment::SummarizedExperiment(cellMatList, colData = sampleData, rowRanges = txList, metadata = countInfo)

    }else{

        cellTypes <- gsub("__.*", "", colnames(newSumMat))
        cellMatList <- lapply(unique(cellTypes), function(x){
    
            cellMat <- subsetMat[, x == cellTypes]
            colnames(cellMat) <- gsub(".*__","", colnames(cellMat))
            cellMat
        })
        names(cellMatList) = cellTypes

        rnaSE <-  SummarizedExperiment::SummarizedExperiment(mat, colData = sampleData, metadata = countInfo)
    }
    
    return(rnaSE)
}