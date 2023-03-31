# bulkLSI works on a 3 sample test dataset

    Code
      LSIObj
    Output
      class: SummarizedExperiment 
      dim: 2 6 
      metadata(6): CellCounts FragmentCounts ... OrgDb Directory
      assays(1): LSI
      rownames(2): PC1 PC2
      rowData names(0):
      colnames(6): C2__scATAC_BMMC_R1 C2__scATAC_CD34_BMMC_R1 ...
        C3__scATAC_CD34_BMMC_R1 C3__scATAC_PBMC_R1
      colData names(4): Sample PassQC CellType Freq

---

    Code
      UMAPvalues
    Output
             UMAP1       UMAP2                  Sample PassQC CellType Freq
      1  0.2267531  0.81612460      C2__scATAC_BMMC_R1      1       C2  379
      2  0.5380543 -0.98632940 C2__scATAC_CD34_BMMC_R1      1       C2    0
      3 -0.4860336  0.06429927      C2__scATAC_PBMC_R1      1       C2    0
      4 -0.4325334  1.03767268      C3__scATAC_BMMC_R1      1       C3   22
      5  0.3410281 -0.19679578 C3__scATAC_CD34_BMMC_R1      1       C3  302
      6 -0.1872684 -0.73497136      C3__scATAC_PBMC_R1      1       C3    0

