# getDifferentialAccessibleTiles works on a 3 sample test dataset

    GRanges object with 11316 ranges and 13 metadata columns:
              seqnames              ranges strand |                   Tile
                 <Rle>           <IRanges>  <Rle> |            <character>
          [1]    chr10 100083000-100083499      * | chr10:100083000-1000..
          [2]    chr10 100206500-100206999      * | chr10:100206500-1002..
          [3]    chr10 101090000-101090499      * | chr10:101090000-1010..
          [4]    chr10 101101500-101101999      * | chr10:101101500-1011..
          [5]    chr10 101279500-101279999      * | chr10:101279500-1012..
          ...      ...                 ...    ... .                    ...
      [11312]     chrY   16196500-16196999      * | chrY:16196500-16196999
      [11313]     chrY   22737500-22737999      * | chrY:22737500-22737999
      [11314]     chrY   22944000-22944499      * | chrY:22944000-22944499
      [11315]     chrY   23153000-23153499      * | chrY:23153000-23153499
      [11316]     chrY     2884500-2884999      * |   chrY:2884500-2884999
              CellPopulation     Foreground          Background   P_value
                 <character>    <character>         <character> <numeric>
          [1]             C3 scATAC_BMMC_R1 scATAC_CD34_BMMC_R1        NA
          [2]             C3 scATAC_BMMC_R1 scATAC_CD34_BMMC_R1        NA
          [3]             C3 scATAC_BMMC_R1 scATAC_CD34_BMMC_R1        NA
          [4]             C3 scATAC_BMMC_R1 scATAC_CD34_BMMC_R1        NA
          [5]             C3 scATAC_BMMC_R1 scATAC_CD34_BMMC_R1        NA
          ...            ...            ...                 ...       ...
      [11312]             C3 scATAC_BMMC_R1 scATAC_CD34_BMMC_R1        NA
      [11313]             C3 scATAC_BMMC_R1 scATAC_CD34_BMMC_R1        NA
      [11314]             C3 scATAC_BMMC_R1 scATAC_CD34_BMMC_R1        NA
      [11315]             C3 scATAC_BMMC_R1 scATAC_CD34_BMMC_R1        NA
      [11316]             C3 scATAC_BMMC_R1 scATAC_CD34_BMMC_R1        NA
              Test_Statistic       FDR  Log2FC_C  MeanDiff Avg_Intensity_Case
                   <numeric> <logical> <numeric> <numeric>          <numeric>
          [1]             NA      <NA>        NA        NA                 NA
          [2]             NA      <NA>        NA        NA                 NA
          [3]             NA      <NA>        NA        NA                 NA
          [4]             NA      <NA>        NA        NA                 NA
          [5]             NA      <NA>        NA        NA                 NA
          ...            ...       ...       ...       ...                ...
      [11312]             NA      <NA>        NA        NA                 NA
      [11313]             NA      <NA>        NA        NA                 NA
      [11314]             NA      <NA>        NA        NA                 NA
      [11315]             NA      <NA>        NA        NA                 NA
      [11316]             NA      <NA>        NA        NA                 NA
              Pct0_Case Avg_Intensity_Control Pct0_Control
              <numeric>             <numeric>    <numeric>
          [1]         1               11.3103            0
          [2]         1               14.3098            0
          [3]         1               11.8950            0
          [4]         1               11.3103            0
          [5]         1               11.3103            0
          ...       ...                   ...          ...
      [11312]         1               11.3103            0
      [11313]         1               14.6317            0
      [11314]         1               12.8948            0
      [11315]         1               11.3103            0
      [11316]         1               11.8950            0
      -------
      seqinfo: 24 sequences from an unspecified genome; no seqlengths

