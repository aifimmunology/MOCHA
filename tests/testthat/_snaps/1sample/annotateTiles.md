# annotateTiles works on a 1 sample test dataset

    Code
      SummarizedExperiment::rowRanges(STM)
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

