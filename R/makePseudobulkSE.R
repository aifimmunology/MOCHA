
## work flow


makePseudobulkSE <- function(SO, cellTypeColumn, sampleColumn, 
                                 filterByRowData = TRUE,
                                 TxDb = TxDb.Hsapiens.UCSC.hg38.refGene, 
                                 OrgDb = org.Hs.eg.db){


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



pbmc@meta.data$sample_celltype = paste(pbmc@meta.data[,'predicted.id'], pbmc@meta.data[,'Sample'], sep = '__')
mat <- Seurat::AverageExpression(pbmc, assays = 'RNA',  group.by = 'sample_celltype', return.seurat = FALSE)
metadata <- MOCHA:::sampleDataFromCellColData(pbmc@meta.data, 'Sample')

SummarizedExperiment(mat, colData = metadata )

sampleColumn = 'Sample'
cellTypeColumn = 'predicted.id'
percentMito = 'percent.mt'
fullMeta <- dplyr::group_by_at(as.data.frame(pbmc@meta.data),c(sampleColumn, cellTypeColumn))
countInfo <- dplyr::summarise(fullMeta, nCount = sum(nCount_RNA), nUMI = sum(nUMI), nGene = sum(nGene), 
                            percent.mito = mean(percent.mito), CellCount =  )

fullMeta <- dplyr::group_by(as.data.frame(pbmc@meta.data), Sample, predicted.id)
countInfo <- dplyr::summarise(fullMeta, nCount = sum(nCount_RNA), nUMI = sum(nUMI), nGene = sum(nGene), 
                            percent.mito = mean(percent.mito), CellCount =  )

countInfo <- dplyr::summarise(fullMeta, nCount = sum(nCount_RNA), nFeature = sum(nFeature_RNA),
                                     percent.mito = mean(vars(percentMito), na.rm = TRUE), CellCount = dplyr::n())

library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txList <- suppressWarnings(GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg38.knownGene, by = ("gene")))
names(txList) <- suppressWarnings(AnnotationDbi::mapIds(OrgDb, names(txList), "SYMBOL", "ENTREZID"))

sum(names(txList) %in% rownames(mat$RNA))

tmp <- rownames(mat$RNA)
tmp <- gsub("\\.*","", 

rownames(mat$RNA)[grepl("\\.",rownames(mat$RNA))]
            
mat <- Seurat::AverageExpression(pbmc, assays = 'RNA', slot = 'counts', group.by = 'sample_celltype', return.seurat = FALSE)
tmp_mat <- as.data.frame(mat$RNA)
colnames(tmp_mat) <- gsub(" ","_", colnames(tmp_mat))
tmp_mat$AllGenes <- rownames(tmp_mat)
tmp_mat$Genes <- gsub("\\.*","", tmp_mat$AllGenes) 


summarizeMat <- dplyr::group_by(tmp_mat, Genes) %>% dplyr::select(!AllGenes) %>% 
                    dplyr::summarise_all(mean)
newSumMat <- as.data.frame(summarizeMat[,-1])
rownames(newSumMat) <- summarizeMat$Genes 
            
txList <- genes(testTxDb)
            
txList$GeneSymbol <- suppressWarnings(AnnotationDbi::mapIds(OrgDb, txList$gene_id, "SYMBOL", "ENSEMBL"))
            
tmp <- stack(genes(TxDb, single.strand.genes.only = F))
tmp$GeneSymbol <- suppressWarnings(AnnotationDbi::mapIds(OrgDb, tmp$gene_id, "SYMBOL", "ENTREZID"))
     
GenomeInfoDb::keepStandardChromosomes(sort(tmp), species='Homo_sapiens',
                                      pruning.mode = 'coarse')

cellTypes <- gsub("__.*", "", colnames(subsetMat))
            
mat <- Seurat::AverageExpression(full_pbmc, assays = 'RNA', slot = 'counts', 
                                 group.by = 'sample_celltype', return.seurat = FALSE)
cellMatList <- lapply(unique(cellTypes), function(x){
    
    cellMat <- subsetMat[, x == cellTypes]
    colnames(cellMat) <- gsub(".*__","", colnames(cellMat))
    cellMat
    })
names(cellMatList) = unique(cellTypes)
        
rownames(subsetMat) == txList$GeneSymbol)

txList$GeneSymbol <- suppressWarnings(AnnotationDbi::mapIds(OrgDb, as.character(txList$name), "SYMBOL", "ENTREZID"))
    
full_pbmc@meta.data$sample_celltype = paste(full_pbmc@meta.data[,cellTypeColumn], full_pbmc@meta.data[,sampleColumn], sep = '__')
 mat <- Seurat::AverageExpression(full_pbmc, assays = 'RNA', slot = 'counts', group.by = 'sample_celltype', return.seurat = FALSE)

metadata <- MOCHA:::sampleDataFromCellColData(full_pbmc@meta.data, 'Sample')
                            
metadata <- MOCHA:::sampleDataFromCellColData(full_pbmc@meta.data, 'Sample')
                            
       cellColDataNoNA <- BiocGenerics::Filter(function(x) {
    !all(is.na(x))
  }, full_pbmc@meta.data)    
                            
cellColDT <- data.table::as.data.table(cellColDataNoNA)
BoolDT <- cellColDT[, lapply(.SD, function(x) {
    length(unique(x)) == 1
  }), by = c('Sample')]
                            
l)

se <- makePseudobulkSE(full_pbmc, 'predicted.id', 'Sample')