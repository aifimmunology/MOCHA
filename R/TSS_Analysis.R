

##################################################################################################
### getStartWSites pulls out all peaks that fall in TSS sites and annotates them with the name of gene.
## @peakSet - GRanges object that contains all the peaks you want annotated
## @TxDb - Transcript Db object for organism. Default is Hg38.
## @Org - Organism Db for annotating gene names across databases. Must match TxDb.

getStartSites <- function(peakSet,
                          TxDb = TxDb.Hsapiens.UCSC.hg38.refGene,
                          Org = org.Hs.eg.db) {
  if (grepl("data.table|SummarizedExperiment", class(peakSet)[1])) {
    Granges <- ExtractGR(peakSet)
  } else if (class(peakSet)[1] == "GRanges") {
    GRanges <- peakSet
  }

  tss1 <- suppressWarnings(ensembldb::transcriptsBy(TxDb, by = ("gene")))

  names(tss1) <- suppressWarnings(mapIds(Org, names(tss1), "SYMBOL", "ENTREZID"))

  allT <- stack(tss1) %>%
    GenomicRanges::trim(.) %>%
    GenomicRanges::promoters(., upstream = 0, downstream = 0) %>%
    plyranges::mutate(exactTSS = start(.)) %>%
    plyranges::filter(!duplicated(exactTSS)) %>%
    plyranges::anchor_3p(.) %>%
    suppressWarnings(plyranges::stretch(., extend = 125)) %>%
    plyranges::unanchor() %>%
    GenomicRanges::trim()

  tpeaks <- plyranges::join_overlap_intersect(allT, GRanges)

  return(tpeaks)
}


##################################################################################################
### getTSS pulls out all peaks that fall in TSS sites and annotates them with the name of gene.
## @completeDAPs - GRanges object that contains the differential measurements across all peaks (unfiltered DAPs)
## @returnAllTSS - Flag to return all TSS sites with DAPs measurements, without filtering for alternative TSS ussage
## @nuancedTSS - True/False flag to determine if alternative TSS genes should be filtered out
##                if all their differential TSS usage falls within too small of a range. Default is TRUE
## @threshold -  Basepair distance Threshold for filtering out genes who's differential TSS sites all falls very close together,
## 		   and thus likely too close to be able to distinguish TSS behavior clearly.
## @effectSize - name of the metadata column with effect size measurements
##              (used to ensure differential TSS sites are not all differential in the same direction).
## @TxDb - Transcript Db object for organism. Default is Hg38.
## @Org - Organism Db for annotating gene names across databases. Must match TxDb.


getAltTSS <- function(completeDAPs,
                      returnAllTSS = FALSE,
                      nuancedTSS = TRUE,
                      nuancedTSSGap = 150,
                      threshold = 0.2,
                      effectSize = HL_EffectSize,
                      TxDb = TxDb.Hsapiens.UCSC.hg38.refGene,
                      Org = org.Hs.eg.db) {
  if (grepl("data.table|SummarizedExperiment", class(completeDAPs)[1])) {
    DAP_Granges <- ExtractGR(completeDAPs)
  } else if (class(completeDAPs)[1] == "GRanges") {
    DAP_GRanges <- completeDAPs
  }

  tss1 <- ensembldb::transcriptsBy(TxDb, by = ("gene"))

  names(tss1) <- mapIds(Org, names(tss1), "SYMBOL", "ENTREZID")

  allT <- stack(tss1) %>%
    trim(.) %>%
    promoters(., upstream = 0, downstream = 0) %>%
    plyranges::mutate(exactTSS = start(.)) %>%
    plyranges::filter(!duplicated(exactTSS)) %>%
    plyranges::anchor_3p(.) %>%
    plyranges::stretch(., extend = 125) %>%
    trim()

  tpeaks <- join_overlap_intersect(allT, DAP_GRanges)

  if (returnAllTSS) {
    return(tpeaks)
  }

  altTSS <- tpeaks %>%
    plyranges::filter(!duplicated(exactTSS)) %>%
    plyranges::group_by(name) %>%
    plyranges::filter(any(FDR < threshold) & any(duplicated(name))) %>%
    plyranges::filter(ifelse(all(!is.na(FDR)) & all(FDR < threshold, na.rm = TRUE),
      !all({{ effectSize }} > 0) & !all({{ effectSize }} < 0), TRUE
    )) %>%
    ungroup() %>%
    sort()

  if (!nuancedTSS) {
    nuancedGenes <- split(altTSS, as.character(altTSS$name)) %>%
      mclapply(., function(x) {
        tmp <- gaps(x) %>%
          filter(seqnames == seqnames(x) & strand == strand(x)) %>%
          width(.)
        ifelse(length(tmp) == 3 & tmp[2] < nuancedTSSGap, FALSE, TRUE)
      }) %>%
      unlist()

    altTSS <- altTSS %>% filter(name %in% nuancedGenes[nuancedGenes])
  }
  return(altTSS)
}
