#' @title \code{makePseudobulkRNA}
#'
#' @description \code{makePseuduoulkRNA} pseodubulks a Seurat object by sample and cell type into a SummarizedExperiment, similar to MOCHA

#'
#' @param SO Seurat Object 
#' @param cellTypeColumn The column of the Seurat object with cell type information
#' @param Sample The column of the Seurat object with sample information
#' @param filterByRowData A boolean flag that determines whether to export a summarizedExperiment with Genomic location for each gene (and thus filter the object to genes that match the database)
#' @param TxDb Transcript database to be used for identifying gene locations. 
#' @param OrgDb Organism database to match up gene names and locations.
#' @return A SummarizedExperiment carrying pseudobulked average expression per 1000 cells for each cell type. 
#'
#'
#' @export
makePseudobulkRNA <- function(SO, cellTypeColumn, sampleColumn = "Sample", 
                                 addGenomicLocations = TRUE,
                                 Seurat_format = 'SYMBOL',
                                 TxDb = "TxDb.Hsapiens.UCSC.hg38.refGene", 
                                 OrgDb = "org.Hs.eg.db"){

    if(any(!c(cellTypeColumn, sampleColumn) %in% colnames(SO@meta.data))){

        stop('Sample or cell type columns are missing from Seurat object. Please verify that the provided column names are correct.')

    }

    #Generate sample and cell type column in metadata for pseudobulking.
    counts <- SO@assays$RNA@counts
    metadata <- seurat@meta.data
    metadata$CellTypeColumn = factor(unlist(metadata[,cellTypeColumn]))
    sce <- SingleCellExperiment(assays = list(counts = counts), 
                        colData = metadata)
    groups <- colData(sce)[, c('cellTypeColumn', sampleColumn)]
    aggr_counts <- aggregate.Matrix(t(counts(sce)), 
                            groupings = groups, fun = "sum") 
    cellTypeList <- unique(metadata$CellTypeColumn)
    counts_ls <- lapply(1:length(cellTypeList), function(i) {
        ## Extract indexes of columns in the global matrix that match a given cluster
        column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == cellTypeList[i])
        ## Store corresponding sub-matrix as one element of a list
        aggr_counts[, column_idx]
    })
    names(counts_ls) <- cellTypeList

    metadata_ls <- lapply(1:length(counts_ls), function(i) {
        ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
        df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
        ## Use tstrsplit() to separate cluster (cell type) and sample IDs
        df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
        df$sample_id  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
        ## Retrieve cell count information for this cluster from global cell count table
        idx <- which(colnames(t) == unique(df$cluster_id))
        cell_counts <- t[, idx]
        ## Match order of cell_counts and sample_ids
        sample_order <- match(df$sample_id, names(cell_counts))
        cell_counts <- cell_counts[sample_order]
        ## Append cell_counts to data frame
        df$cell_count <- cell_counts
        ## Join data frame (capturing metadata specific to cluster) to generic metadata
        df <- plyr::join(df, metadata, 
                        by = intersect(names(df), names(metadata)))
        ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
        rownames(df) <- df$cluster_sample_id
        df
    })

    names(metadata_ls) <- cellTypeList

    mat_df <- as.data.frame(mat$RNA)
    colnames(mat_df) <- gsub(" ","_", colnames(mat_df))
    mat_df$AllGenes <- rownames(mat_df)
    ## Clean up isoforms, if any (duplicate gene names. More likely an issue with ensembl IDs.)
    mat_df$Genes <- gsub("\\.*","", mat_df$AllGenes) 

    #Multiply by 1000 to get the average gene counts (normalized, presumably) per 1000 cells. 
    newSumMat <- as.data.frame(summarizeMat[,-1])*1000
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

    #Process the meta data and extract the total counts and features for each population, as well as cell counts. 
    fullMeta <- dplyr::group_by_at(as.data.frame(SO@meta.data), c(sampleColumn, cellTypeColumn))

    countInfo <- dplyr::summarise(fullMeta, nCount =  mean(nCount_RNA, na.rm = TRUE), nFeature = mean(nFeature_RNA, na.rm = TRUE),
                                      CellCount = dplyr::n())

    sampleList <- unique(SO@meta.data[,sampleColumn])

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
        txList <- plyranges::filter(txList, !is.na(GeneSymbol) & GeneSymbol %in% rownames(newSumMat))
        txList <- plyranges::reduce_ranges(plyranges::group_by(txList, GeneSymbol))
        names(txList) <- txList$GeneSymbol

        #Subset down the read count matrix to just the transcripts we ahve in this database.
        subsetMat <- newSumMat[rownames(newSumMat) %in% txList$GeneSymbol, ]
        subsetMat <- subsetMat[match(txList$GeneSymbol,rownames(subsetMat)), ]

        #now iterate over each cell type, and generate a sample-gene matrix (while adding in empty columns for samples that didn't have any of a given cell type detected.)
        #Failing to add in empty samples will break the SummarizedExperiment object at the end. 
        cellTypes <- gsub("__.*", "", colnames(subsetMat))
        cellMatList <- lapply(unique(cellTypes), function(x){
            
            cellMat <- subsetMat[, x == cellTypes]
            colnames(cellMat) <- gsub(".*__","", colnames(cellMat))
            if(!all(sampleList %in% colnames(cellMat))){
                    emptyData = rep(0, dim(cellMat)[1])
                    emptyDF = do.call('cbind',lapply(1:sum(!sampleList %in% colnames(cellMat)), function(x) emptyData))
                    colnames(emptyDF) = sampleList[!sampleList %in% colnames(cellMat)]
                    cellMat <- cbind(cellMat, emptyDF)
            }
            cellMat[,match(unlist(sampleData[,sampleColumn]), colnames(cellMat))]
        })
        names(cellMatList) = unique(cellTypes)

        # Repackage the gene-sample matrices, the sampleData, transcript Genomic Ranges, and associated metadata (countInfo) into one SummarizedExperiment object. 
        rnaSE <- SummarizedExperiment::SummarizedExperiment(cellMatList, colData = sampleData, rowRanges = txList, metadata = countInfo)

    }else{

        #now iterate over each cell type, and generate a sample-gene matrix (while adding in empty columns for samples that didn't have any of a given cell type detected.)
        #Failing to add in empty samples will break the SummarizedExperiment object at the end. 

        cellTypes <- gsub("__.*", "", colnames(newSumMat))
        cellMatList <- lapply(unique(cellTypes), function(x){
    
            cellMat <- subsetMat[, x == cellTypes]
            colnames(cellMat) <- gsub(".*__","", colnames(cellMat))
            if(!all(sampleList %in% colnames(cellMat))){
                    emptyData = rep(0, dim(cellMat)[1])
                    emptyDF = do.call('cbind',lapply(1:sum(!sampleList %in% colnames(cellMat)), function(x) emptyData))
                    colnames(emptyDF) = sampleList[!sampleList %in% colnames(cellMat)]
                    cellMat <- cbind(cellMat, emptyDF)
            }
            cellMat[,match(unlist(sampleData[,sampleColumn]), colnames(cellMat))]
        })
        names(cellMatList) = cellTypes

        # Repackage the gene-sample matrices, the sampleData,  and associated metadata (countInfo) into one SummarizedExperiment object. 
        rnaSE <-  SummarizedExperiment::SummarizedExperiment(mat, colData = sampleData, metadata = countInfo)
    }
    
    return(rnaSE)
}