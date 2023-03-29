# getCoverage works on a 1 sample test dataset

    Code
      covFiles
    Output
      $`C2#PBMCSmall__0.117146`
      GRanges object with 222676 ranges and 1 metadata column:
                 seqnames              ranges strand |     score
                    <Rle>           <IRanges>  <Rle> | <numeric>
             [1]     chr1            1-521571      * |         0
             [2]     chr1       521572-521598      * |         1
             [3]     chr1       521599-565288      * |         0
             [4]     chr1       565289-565304      * |         1
             [5]     chr1       565305-565335      * |         2
             ...      ...                 ...    ... .       ...
        [222672]     chr2 243040419-243040586      * |         1
        [222673]     chr2 243040587-243041454      * |         0
        [222674]     chr2 243041455-243041513      * |         1
        [222675]     chr2 243041514-243046389      * |         0
        [222676]     chr2 243046390-243046418      * |         1
        -------
        seqinfo: 2 sequences from hg38 genome
      
      $`C5#PBMCSmall__0.171033`
      GRanges object with 307923 ranges and 1 metadata column:
                 seqnames              ranges strand |     score
                    <Rle>           <IRanges>  <Rle> | <numeric>
             [1]     chr1            1-565290      * |         0
             [2]     chr1       565291-565293      * |         1
             [3]     chr1       565294-565329      * |         3
             [4]     chr1       565330-565345      * |         2
             [5]     chr1              565346      * |         1
             ...      ...                 ...    ... .       ...
        [307919]     chr2 243031330-243031590      * |         1
        [307920]     chr2 243031591-243034011      * |         0
        [307921]     chr2 243034012-243034339      * |         1
        [307922]     chr2 243034340-243073424      * |         0
        [307923]     chr2 243073425-243073624      * |         1
        -------
        seqinfo: 2 sequences from hg38 genome
      

# getCoverage works on a 1 sample test dataset with filterEmpty=TRUE

    Code
      covFiles
    Output
      $`C2#PBMCSmall__0.117146`
      GRanges object with 169007 ranges and 1 metadata column:
                 seqnames              ranges strand |     score
                    <Rle>           <IRanges>  <Rle> | <numeric>
             [1]     chr1       521572-521598      * |         1
             [2]     chr1       565289-565304      * |         1
             [3]     chr1       565305-565335      * |         2
             [4]     chr1       565336-565349      * |         1
             [5]     chr1       569375-569378      * |         2
             ...      ...                 ...    ... .       ...
        [169003]     chr2 243038802-243038972      * |         1
        [169004]     chr2 243039871-243040072      * |         1
        [169005]     chr2 243040419-243040586      * |         1
        [169006]     chr2 243041455-243041513      * |         1
        [169007]     chr2 243046390-243046418      * |         1
        -------
        seqinfo: 2 sequences from hg38 genome
      
      $`C5#PBMCSmall__0.171033`
      GRanges object with 258866 ranges and 1 metadata column:
                 seqnames              ranges strand |     score
                    <Rle>           <IRanges>  <Rle> | <numeric>
             [1]     chr1       565291-565293      * |         1
             [2]     chr1       565294-565329      * |         3
             [3]     chr1       565330-565345      * |         2
             [4]     chr1              565346      * |         1
             [5]     chr1       569387-569412      * |         1
             ...      ...                 ...    ... .       ...
        [258862]     chr2 243031253-243031269      * |         3
        [258863]     chr2 243031270-243031329      * |         2
        [258864]     chr2 243031330-243031590      * |         1
        [258865]     chr2 243034012-243034339      * |         1
        [258866]     chr2 243073425-243073624      * |         1
        -------
        seqinfo: 2 sequences from hg38 genome
      

# getSpecificCoverage works on a 1 sample test dataset

    Code
      counts
    Output
      $`C2#PBMCSmall__0.117146`
      GRanges object with 2 ranges and 1 metadata column:
            seqnames        ranges strand |     score
               <Rle>     <IRanges>  <Rle> | <numeric>
        [1]     chr1 565291-565349      * |      30.0
        [2]     chr1 569375-569412      * |      30.2
        -------
        seqinfo: 2 sequences from hg38 genome
      
      $`C5#PBMCSmall__0.171033`
      GRanges object with 4 ranges and 1 metadata column:
            seqnames              ranges strand |     score
               <Rle>           <IRanges>  <Rle> | <numeric>
        [1]     chr1       565291-565346      * |        36
        [2]     chr1       569387-569412      * |        26
        [3]     chr2 243031253-243031590      * |       144
        [4]     chr2 243034012-243034339      * |       328
        -------
        seqinfo: 2 sequences from hg38 genome
      

