# bulkDimReduction works on a 3 sample test dataset

    Code
      LSIObj
    Output
      class: SummarizedExperiment 
      dim: 2 6 
      metadata(6): CellCounts FragmentCounts ... OrgDb Directory
      assays(1): LSI
      rownames(2): LSI1 LSI2
      rowData names(0):
      colnames(6): C2__scATAC_BMMC_R1 C2__scATAC_CD34_BMMC_R1 ...
        C3__scATAC_CD34_BMMC_R1 C3__scATAC_PBMC_R1
      colData names(5): Sample PassQC CellType Freq FragNumber

---

    Code
      UMAPvalues
    Output
              UMAP1      UMAP2                  Sample PassQC CellType Freq
      1 -0.66386342  0.4887874      C2__scATAC_BMMC_R1      1       C2  379
      2 -0.38785958 -0.2135089 C2__scATAC_CD34_BMMC_R1      1       C2    0
      3  0.08295703 -0.8682348      C2__scATAC_PBMC_R1      1       C2    0
      4 -0.14599729  1.0024571      C3__scATAC_BMMC_R1      1       C3   22
      5  0.68430209 -0.6597805 C3__scATAC_CD34_BMMC_R1      1       C3  302
      6  0.43046117  0.2502796      C3__scATAC_PBMC_R1      1       C3    0
        FragNumber
      1    1040427
      2         NA
      3         NA
      4      44958
      5     966249
      6         NA

