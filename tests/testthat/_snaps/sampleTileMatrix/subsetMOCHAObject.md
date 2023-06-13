# We can subset a sampleTileMatrix object by celltypes

    class: RangedSummarizedExperiment 
    dim: 22775 3 
    metadata(6): summarizedData Genome ... Directory History
    assays(1): C3
    rownames(22775): chr10:100083000-100083499 chr10:100206500-100206999
      ... chrY:23153000-23153499 chrY:2884500-2884999
    rowData names(2): C2 C3
    colnames(3): scATAC_BMMC_R1 scATAC_CD34_BMMC_R1 scATAC_PBMC_R1
    colData names(2): Sample PassQC

---

    class: SummarizedExperiment 
    dim: 1 3 
    metadata(0):
    assays(16): CellCounts FragmentCounts ... DoubletEnrichment
      BlacklistRatio
    rownames(1): C3
    rowData names(0):
    colnames(3): scATAC_BMMC_R1 scATAC_CD34_BMMC_R1 scATAC_PBMC_R1
    colData names(2): Sample PassQC

# We can subset a sampleTileMatrix object by Sample

    class: RangedSummarizedExperiment 
    dim: 22775 1 
    metadata(6): summarizedData Genome ... Directory History
    assays(2): C2 C3
    rownames(22775): chr10:100083000-100083499 chr10:100206500-100206999
      ... chrY:23153000-23153499 chrY:2884500-2884999
    rowData names(2): C2 C3
    colnames(1): scATAC_CD34_BMMC_R1
    colData names(2): Sample PassQC

---

    class: SummarizedExperiment 
    dim: 2 1 
    metadata(0):
    assays(16): CellCounts FragmentCounts ... DoubletEnrichment
      BlacklistRatio
    rownames(2): C2 C3
    rowData names(0):
    colnames(1): scATAC_CD34_BMMC_R1
    colData names(2): Sample PassQC

