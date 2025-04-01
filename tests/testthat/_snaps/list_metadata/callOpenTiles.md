# We can call peaks independent of ArchR

    Code
      tiles@metadata
    Output
      $summarizedData
      class: SummarizedExperiment 
      dim: 2 1 
      metadata(0):
      assays(3): CellCounts FragmentCounts nFrags
      rownames(2): C2 C5
      rowData names(0):
      colnames(1): PBMCSmall
      colData names(1): Sample
      
      $Genome
      [1] "BSgenome.Hsapiens.UCSC.hg19"
      
      $TxDb
      $TxDb$pkgname
      [1] "TxDb.Hsapiens.UCSC.hg38.refGene"
      
      $TxDb$metadata
                                       name
      1                             Db type
      2                  Supporting package
      3                         Data source
      4                              Genome
      5                            Organism
      6                         Taxonomy ID
      7                          UCSC Table
      8                          UCSC Track
      9                        Resource URL
      10                    Type of Gene ID
      11                       Full dataset
      12                   miRBase build ID
      13                  Nb of transcripts
      14                      Db created by
      15                      Creation time
      16 txdbmaker version at creation time
      17   RSQLite version at creation time
      18                    DBSCHEMAVERSION
                                                value
      1                                          TxDb
      2                               GenomicFeatures
      3                                          UCSC
      4                                          hg38
      5                                  Homo sapiens
      6                                          9606
      7                                       refGene
      8                                   UCSC RefSeq
      9                      https://genome.ucsc.edu/
      10                               Entrez Gene ID
      11                                          yes
      12                                         <NA>
      13                                        88819
      14          txdbmaker package from Bioconductor
      15 2024-03-23 22:26:48 +0000 (Sat, 23 Mar 2024)
      16                                       0.99.7
      17                                        2.3.5
      18                                          1.2
      
      
      $OrgDb
      $OrgDb$pkgname
      [1] "org.Hs.eg.db"
      
      $OrgDb$metadata
                       name                                                  value
      1     DBSCHEMAVERSION                                                    2.1
      2             Db type                                                  OrgDb
      3  Supporting package                                          AnnotationDbi
      4            DBSCHEMA                                               HUMAN_DB
      5            ORGANISM                                           Homo sapiens
      6             SPECIES                                                  Human
      7        EGSOURCEDATE                                             2024-Sep20
      8        EGSOURCENAME                                            Entrez Gene
      9         EGSOURCEURL                   ftp://ftp.ncbi.nlm.nih.gov/gene/DATA
      10          CENTRALID                                                     EG
      11              TAXID                                                   9606
      12       GOSOURCENAME                                                       
      13        GOSOURCEURL                                                       
      14       GOSOURCEDATE                                                       
      15     GOEGSOURCEDATE                                             2024-Sep20
      16     GOEGSOURCENAME                                            Entrez Gene
      17      GOEGSOURCEURL                   ftp://ftp.ncbi.nlm.nih.gov/gene/DATA
      18     KEGGSOURCENAME                                            KEGG GENOME
      19      KEGGSOURCEURL                   ftp://ftp.genome.jp/pub/kegg/genomes
      20     KEGGSOURCEDATE                                             2011-Mar15
      21       GPSOURCENAME              UCSC Genome Bioinformatics (Homo sapiens)
      22        GPSOURCEURL ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/database
      23       GPSOURCEDATE                                             2024-Sep22
      24       ENSOURCEDATE                                             2024-May14
      25       ENSOURCENAME                                                Ensembl
      26        ENSOURCEURL                ftp://ftp.ensembl.org/pub/current_fasta
      27       UPSOURCENAME                                                Uniprot
      28        UPSOURCEURL                                http://www.UniProt.org/
      29       UPSOURCEDATE                               Mon Sep 23 15:46:45 2024
      
      
      $History
      $History[[1]]
      [1] "callOpenTiles 1.2.0"
      
      

