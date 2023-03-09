# annotateTiles works on a 3 sample test dataset

    Code
      rowRanges(STM)
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
      GRanges object with 25112 ranges and 4 metadata columns:
                                 seqnames              ranges strand |        C2
                                    <Rle>           <IRanges>  <Rle> | <logical>
            chr1:1000000-1000499     chr1     1000000-1000499      * |     FALSE
          chr1:10002000-10002499     chr1   10002000-10002499      * |     FALSE
          chr1:10002500-10002999     chr1   10002500-10002999      * |      TRUE
        chr1:100026000-100026499     chr1 100026000-100026499      * |     FALSE
          chr1:10003000-10003499     chr1   10003000-10003499      * |      TRUE
                             ...      ...                 ...    ... .       ...
          chr2:99952500-99952999     chr2   99952500-99952999      * |      TRUE
          chr2:99953000-99953499     chr2   99953000-99953499      * |      TRUE
          chr2:99953500-99953999     chr2   99953500-99953999      * |      TRUE
          chr2:99954000-99954499     chr2   99954000-99954499      * |      TRUE
          chr2:99954500-99954999     chr2   99954500-99954999      * |      TRUE
                                        C5    tileType        Gene
                                 <logical> <character> <character>
            chr1:1000000-1000499      TRUE    Promoter        HES4
          chr1:10002000-10002499      TRUE  Intragenic        RBP7
          chr1:10002500-10002999      TRUE  Intragenic        RBP7
        chr1:100026000-100026499      TRUE  Intragenic     SLC35A3
          chr1:10003000-10003499      TRUE  Intragenic        RBP7
                             ...       ...         ...         ...
          chr2:99952500-99952999      TRUE  Intragenic        AFF3
          chr2:99953000-99953499      TRUE  Intragenic        AFF3
          chr2:99953500-99953999      TRUE  Intragenic        AFF3
          chr2:99954000-99954499      TRUE  Intragenic        AFF3
          chr2:99954500-99954999     FALSE  Intragenic        AFF3
        -------
        seqinfo: 2 sequences from an unspecified genome; no seqlengths

