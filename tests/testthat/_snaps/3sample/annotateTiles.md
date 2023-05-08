# annotateTiles works on a 3 sample test dataset

    Code
      SummarizedExperiment::rowRanges(STM)
    Output
      GRanges object with 22775 ranges and 4 metadata columns:
                                  seqnames              ranges strand |        C2
                                     <Rle>           <IRanges>  <Rle> | <logical>
        chr10:100083000-100083499    chr10 100083000-100083499      * |     FALSE
        chr10:100206500-100206999    chr10 100206500-100206999      * |     FALSE
        chr10:100228000-100228499    chr10 100228000-100228499      * |      TRUE
        chr10:101006500-101006999    chr10 101006500-101006999      * |      TRUE
        chr10:101090000-101090499    chr10 101090000-101090499      * |     FALSE
                              ...      ...                 ...    ... .       ...
           chrY:16196500-16196999     chrY   16196500-16196999      * |     FALSE
           chrY:22737500-22737999     chrY   22737500-22737999      * |     FALSE
           chrY:22944000-22944499     chrY   22944000-22944499      * |     FALSE
           chrY:23153000-23153499     chrY   23153000-23153499      * |     FALSE
             chrY:2884500-2884999     chrY     2884500-2884999      * |     FALSE
                                         C3    tileType        Gene
                                  <logical> <character> <character>
        chr10:100083000-100083499      TRUE    Promoter        CPN1
        chr10:100206500-100206999      TRUE  Intragenic        CHUK
        chr10:100228000-100228499     FALSE  Intragenic        CHUK
        chr10:101006500-101006999     FALSE  Intragenic       LZTS2
        chr10:101090000-101090499      TRUE  Intragenic      TLX1NB
                              ...       ...         ...         ...
           chrY:16196500-16196999      TRUE      Distal        <NA>
           chrY:22737500-22737999      TRUE      Distal        <NA>
           chrY:22944000-22944499      TRUE  Intragenic       TTTY4
           chrY:23153000-23153499      TRUE  Intragenic  DAZ1, DAZ4
             chrY:2884500-2884999      TRUE      Distal        <NA>
        -------
        seqinfo: 24 sequences from an unspecified genome; no seqlengths

# annotateTiles works on a GRanges

    Code
      ranges
    Output
      GRanges object with 22775 ranges and 4 metadata columns:
                                  seqnames              ranges strand |        C2
                                     <Rle>           <IRanges>  <Rle> | <logical>
        chr10:100083000-100083499    chr10 100083000-100083499      * |     FALSE
        chr10:100206500-100206999    chr10 100206500-100206999      * |     FALSE
        chr10:100228000-100228499    chr10 100228000-100228499      * |      TRUE
        chr10:101006500-101006999    chr10 101006500-101006999      * |      TRUE
        chr10:101090000-101090499    chr10 101090000-101090499      * |     FALSE
                              ...      ...                 ...    ... .       ...
           chrY:16196500-16196999     chrY   16196500-16196999      * |     FALSE
           chrY:22737500-22737999     chrY   22737500-22737999      * |     FALSE
           chrY:22944000-22944499     chrY   22944000-22944499      * |     FALSE
           chrY:23153000-23153499     chrY   23153000-23153499      * |     FALSE
             chrY:2884500-2884999     chrY     2884500-2884999      * |     FALSE
                                         C3    tileType        Gene
                                  <logical> <character> <character>
        chr10:100083000-100083499      TRUE    Promoter        CPN1
        chr10:100206500-100206999      TRUE  Intragenic        CHUK
        chr10:100228000-100228499     FALSE  Intragenic        CHUK
        chr10:101006500-101006999     FALSE  Intragenic       LZTS2
        chr10:101090000-101090499      TRUE  Intragenic      TLX1NB
                              ...       ...         ...         ...
           chrY:16196500-16196999      TRUE      Distal        <NA>
           chrY:22737500-22737999      TRUE      Distal        <NA>
           chrY:22944000-22944499      TRUE  Intragenic       TTTY4
           chrY:23153000-23153499      TRUE  Intragenic  DAZ1, DAZ4
             chrY:2884500-2884999      TRUE      Distal        <NA>
        -------
        seqinfo: 24 sequences from an unspecified genome; no seqlengths

