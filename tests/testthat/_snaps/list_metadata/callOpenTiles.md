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
      [1] "hg19"
      
      $TxDb
      $TxDb$pkgname
      [1] "TxDb.Hsapiens.UCSC.hg38.refGene"
      
      $TxDb$metadata
                                             name
      1                                   Db type
      2                        Supporting package
      3                               Data source
      4                                    Genome
      5                                  Organism
      6                               Taxonomy ID
      7                                UCSC Table
      8                                UCSC Track
      9                              Resource URL
      10                          Type of Gene ID
      11                             Full dataset
      12                         miRBase build ID
      13                        Nb of transcripts
      14                            Db created by
      15                            Creation time
      16 GenomicFeatures version at creation time
      17         RSQLite version at creation time
      18                          DBSCHEMAVERSION
                                                value
      1                                          TxDb
      2                               GenomicFeatures
      3                                          UCSC
      4                                          hg38
      5                                  Homo sapiens
      6                                          9606
      7                                       refGene
      8                                   NCBI RefSeq
      9                       http://genome.ucsc.edu/
      10                               Entrez Gene ID
      11                                          yes
      12                                         <NA>
      13                                        88819
      14    GenomicFeatures package from Bioconductor
      15 2023-09-20 17:26:31 +0000 (Wed, 20 Sep 2023)
      16                                       1.53.2
      17                                        2.3.1
      18                                          1.2
      
      
      $OrgDb
      $OrgDb$pkgname
      [1] "org.Hs.eg.db"
      
      $OrgDb$metadata
                       name                                                 value
      1     DBSCHEMAVERSION                                                   2.1
      2             Db type                                                 OrgDb
      3  Supporting package                                         AnnotationDbi
      4            DBSCHEMA                                              HUMAN_DB
      5            ORGANISM                                          Homo sapiens
      6             SPECIES                                                 Human
      7        EGSOURCEDATE                                            2023-Sep11
      8        EGSOURCENAME                                           Entrez Gene
      9         EGSOURCEURL                  ftp://ftp.ncbi.nlm.nih.gov/gene/DATA
      10          CENTRALID                                                    EG
      11              TAXID                                                  9606
      12       GOSOURCENAME                                         Gene Ontology
      13        GOSOURCEURL http://current.geneontology.org/ontology/go-basic.obo
      14       GOSOURCEDATE                                            2023-07-27
      15     GOEGSOURCEDATE                                            2023-Sep11
      16     GOEGSOURCENAME                                           Entrez Gene
      17      GOEGSOURCEURL                  ftp://ftp.ncbi.nlm.nih.gov/gene/DATA
      18     KEGGSOURCENAME                                           KEGG GENOME
      19      KEGGSOURCEURL                  ftp://ftp.genome.jp/pub/kegg/genomes
      20     KEGGSOURCEDATE                                            2011-Mar15
      21       GPSOURCENAME             UCSC Genome Bioinformatics (Homo sapiens)
      22        GPSOURCEURL                                                      
      23       GPSOURCEDATE                                            2023-Aug20
      24       ENSOURCEDATE                                            2023-May10
      25       ENSOURCENAME                                               Ensembl
      26        ENSOURCEURL               ftp://ftp.ensembl.org/pub/current_fasta
      27       UPSOURCENAME                                               Uniprot
      28        UPSOURCEURL                               http://www.UniProt.org/
      29       UPSOURCEDATE                              Mon Sep 18 16:12:39 2023
      
      
      $History
      $History[[1]]
      [1] "callOpenTiles 1.1.0"
      
      

